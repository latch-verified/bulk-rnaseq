import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple, Union

from latch import large_task, message
from latch.types import LatchDir, LatchFile, file_glob
from wf.core.errors import TrimgaloreError
from wf.core.models import (PairedEndReads, Replicate, SingleEndReads,
                            TrimgaloreSalmonInput, TrimgaloreSalmonOutput)
from wf.core.utils import (concatenate_files, run, run_and_capture_output,
                           slugify)

# Each Salmon warning or error log starts with a timestamp surrounded in square
# brackets ('\d{4}' represents the first part of the timestamp - the year)
SALMON_ALERT_PATTERN = re.compile(r"\[(warning|error)\] (.+?)(?:\[\d{4}|$)", re.DOTALL)

DOCKERFILE_PATH = Path(__file__).parent / "Dockerfile"


def merge_replicates(
    replicates: List[Replicate], sample_name: str
) -> Union[Tuple[Path], Tuple[Path, Path]]:
    local_r1_path = f"{slugify(sample_name)}_r1_merged.fq"
    r1 = concatenate_files((str(x.r1.path) for x in replicates), local_r1_path)

    if isinstance(replicates[0], SingleEndReads):
        return (r1,)

    if not all(isinstance(x, PairedEndReads) for x in replicates):
        raise RuntimeError("Not all technical replicates were paired end")

    local_r2_path = f"{slugify(sample_name)}_r2_merged.fq"
    r2 = concatenate_files((str(x.r2.path) for x in replicates), local_r2_path)
    return (r1, r2)


@large_task(dockerfile=DOCKERFILE_PATH)
def trimgalore_salmon(input: TrimgaloreSalmonInput) -> Optional[TrimgaloreSalmonOutput]:

    SALMON_DIR = "/root/salmon_quant/"
    """Default location for salmon outputs."""

    REMOTE_PATH = f"latch:///{input.base_remote_output_dir}{input.run_name}/Quantification (salmon)/{input.sample_name}/"
    """Remote path prefix for LatchFiles + LatchDirs"""

    #
    # Trimgalore (separate!)
    try:
        outputs = [do_trimgalore(input, i, x) for i, x in enumerate(input.replicates)]
    except TrimgaloreError as e:
        print(f"Handling failure in trimming {input.sample_name} gracefully.")
        print(f"\tTrimming error ~ {e}")
        return None

    trimmed_replicates = [x[1] for x in outputs]
    trimgalore_reports = [y for x in outputs for y in x[0]]

    merged = merge_replicates(trimmed_replicates, input.sample_name)
    if len(merged) == 1:
        (r1,) = merged
        reads = ["-r", str(r1)]
    else:
        r1, r2 = merged
        reads = ["-1", str(r1), "-2", str(r2)]

    # Free space.
    for rep in trimmed_replicates:
        os.remove(rep.r1.path)
        if isinstance(rep, PairedEndReads):
            os.remove(rep.r2.path)

    #
    # Salmon (separate!)

    quant_cmd = [
        "salmon",
        "quant",
        "-i",
        input.salmon_index.local_path,
        "-l",
        "A",
        *reads,
        "--threads",
        str(96),
        "--validateMappings",
        "-o",
        SALMON_DIR,
    ]

    try:
        returncode, stdout = run_and_capture_output(quant_cmd)
    except subprocess.CalledProcessError as _:

        def parse_salmon_warning(
            alert_message: str, input: TrimgaloreSalmonInput
        ) -> str:
            if "of fragments were shorter than the k" in alert_message:
                return alert_message
            elif "Detected a *potential* strand bias" in alert_message:
                default = "salmon_quant/lib_format_counts.json"
                return alert_message.replace(
                    default, REMOTE_PATH + "lib_format_counts.json"
                )
            return alert_message

        identifier = f"sample {input.sample_name}"
        errors = []
        for alert_type, alert_message in re.findall(SALMON_ALERT_PATTERN, stdout):
            title = f"Salmon {alert_type} for {identifier}"
            if alert_type == "warning":
                message(
                    "warning",
                    {
                        "title": title,
                        "body": parse_salmon_warning(alert_message, input),
                    },
                )
            else:
                message("error", {"title": title, "body": alert_message})
                errors.append(alert_message)

        if returncode != 0:
            deets = "\n".join(["Error(s) occurred while running Salmon", *errors])
            return None

    # Also moves count files out of auxilliary directory.
    quant_path = f"/root/{slugify(input.sample_name)}_quant.sf"
    salmon_quant = Path(f"{SALMON_DIR}/quant.sf").rename(quant_path)

    junc_path = Path(f"/root/{input.sample_name}.bam.junc")
    if input.run_splicing:
        try:
            # TODO - gaw so bad
            print("\tDownloading STAR map.")
            run(
                [
                    "aws",
                    "s3",
                    "sync",
                    "s3://latch-genomes/Homo_sapiens/RefSeq/GRCh38.p14/STAR_index/",
                    "STAR_index",
                    "--no-progress",  # shh
                ]
            )

            print(f"\tMaking splice-aware alignment map {input.sample_name}")
            run(
                [
                    "STAR",
                    "--genomeDir",
                    "STAR_index",
                    "--twopassMode",
                    "Basic",
                    "--outSAMstrandField",
                    "intronMotif",
                    "--outSAMtype",
                    "BAM",
                    "SortedByCoordinate",
                    "--readFilesIn",
                    *[str(read) for read in merged],
                    "--runThreadN",
                    "96",
                ]
            )

            print(f"\tIndexing splice-aware alignment map {input.sample_name}")
            run(
                [
                    "samtools",
                    "index",
                    "Aligned.sortedByCoord.out.bam",
                    "-@",
                    "96",
                ]
            )

            print(f"\tIndexing splice-aware alignment map {input.sample_name}")
            run(
                [
                    "regtools",
                    "junctions",
                    "extract",
                    "Aligned.sortedByCoord.out.bam",
                    "-s",
                    "0",  # TODO strandedness
                    "-o",
                    str(junc_path),
                ]
            )

            print("\tDone with junctions")

        except subprocess.CalledProcessError as e:
            message(
                "error",
                {
                    "title": f"junction construction error for {input.sample_name}",
                    "body": str(e),
                },
            )
            print(f"Error building junction file: {e}")
            return None

    # Free space.
    for path in merged:
        os.remove(path)

    try:
        tximport_output_path = Path(f"{SALMON_DIR}/gene_abundance.sf")
        run(
            [
                "/root/wf/runTxImport.R",
                "--args",
                str(salmon_quant),
                input.gtf_path,
                tximport_output_path,
            ]
        )
    except subprocess.CalledProcessError as e:
        identifier = f"sample {input.sample_name}"
        message("error", {"title": f"tximport error for {identifier}", "body": str(e)})
        print(
            f"Unable to produce gene mapping from tximport. Error surfaced from tximport -> {e}"
        )
        print(f"Tximport error: {e}")
        return None

    Path(f"{SALMON_DIR}/cmd_info.json").resolve().unlink()
    shutil.rmtree(Path(f"{SALMON_DIR}/logs").resolve())

    if input.run_splicing:
        junction_file = LatchFile(junc_path, REMOTE_PATH + junc_path.name)
    else:
        # temp, Optional field in json data class no work
        junction_file = LatchFile(
            "/root/wf/run_tximport.R", REMOTE_PATH + junc_path.name
        )
    return TrimgaloreSalmonOutput(
        passed_salmon=True,
        passed_tximport=True,
        sample_name=input.sample_name,
        salmon_aux_output=LatchDir(SALMON_DIR, REMOTE_PATH),
        salmon_quant_file=LatchFile(salmon_quant, REMOTE_PATH + salmon_quant.name),
        junction_file=junction_file,
        gene_abundance_file=LatchFile(
            tximport_output_path, REMOTE_PATH + tximport_output_path.name
        ),
        trimgalore_reports=trimgalore_reports,
    )


def do_trimgalore(
    ts_input: TrimgaloreSalmonInput,
    replicate_index: int,
    reads: Replicate,
) -> Tuple[List[LatchFile], Replicate]:
    def _flag(name: str) -> List[str]:
        value = getattr(ts_input, name)
        return [f"--{name}", value] if value is not None else []

    flags = [*_flag("clip_r1"), *_flag("three_prime_clip_r1")]
    read_paths = [str(reads.r1.local_path)]
    if isinstance(reads, PairedEndReads):
        flags += ["--paired", *_flag("clip_r2"), *_flag("three_prime_clip_r2")]
        read_paths.append(str(reads.r2.local_path))

    local_output = f"{slugify(ts_input.sample_name)}_replicate_{replicate_index}"
    returncode, stdout = run_and_capture_output(
        [
            "trim_galore",
            "--cores",
            str(8),
            "--dont_gzip",
            "--output_dir",
            f"./{local_output}",
            *flags,
            *read_paths,
        ]
    )

    # todo(rohankan): examine trimgalore for useful warnings and add them here
    if returncode != 0:
        stdout = stdout.rstrip()
        stdout = stdout[stdout.rindex("\n") + 1 :]
        assert reads.r1.remote_path is not None
        path_name = reads.r1.remote_path.split("/")[-1]
        identifier = f"sample {ts_input.sample_name}, replicate {path_name}"
        message(
            "error",
            {"title": f"Trimgalore error for {identifier}", "body": stdout},
        )
        raise TrimgaloreError(stdout)

    def remote(middle: str) -> str:
        base = f"{ts_input.base_remote_output_dir}{ts_input.run_name}"
        tail = f"{ts_input.sample_name}/replicate_{replicate_index}/"
        return f"latch:///{base}/Quality Control Data/Trimming {middle} (TrimGalore)/{tail}"

    reads_directory = remote("Reads")
    if isinstance(reads, SingleEndReads):
        (r1,) = file_glob(f"{local_output}/*trimmed.fq*", reads_directory)
        trimmed_replicate = SingleEndReads(r1=r1)
    else:
        # File glob sorts files alphanumerically
        r1, r2 = file_glob(f"{local_output}/*val*.fq*", reads_directory)
        trimmed_replicate = PairedEndReads(r1=r1, r2=r2)

    # Delete unneeded files to free disk space
    os.remove(reads.r1.local_path)
    if isinstance(reads, PairedEndReads):
        os.remove(reads.r2.local_path)

    reports_directory = remote("Reports")
    reports = file_glob("*trimming_report.txt", reports_directory)

    return reports, trimmed_replicate
