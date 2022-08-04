"""latch/rnaseq"""

import csv
import os
import re
import shutil
import subprocess
import types
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Iterable, List, Optional, Tuple, Union
from urllib.parse import urlparse

import lgenome
from dataclasses_json import dataclass_json
from flytekit import LaunchPlan, task
from flytekitplugins.pod import Pod
from kubernetes.client.models import (V1Container, V1PodSpec,
                                      V1ResourceRequirements, V1Toleration)
from latch import map_task, message, small_task, workflow
from latch.types import LatchDir, LatchFile


def _capture_output(command: List[str]) -> Tuple[int, str]:
    captured_stdout = []

    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
        universal_newlines=True,
    ) as process:
        assert process.stdout is not None
        for line in process.stdout:
            print(line)
            captured_stdout.append(line)
        process.wait()
        returncode = process.returncode

    return returncode, "\n".join(captured_stdout)


def run(command: List[str], check: bool = True, capture_output: bool = False):
    return subprocess.run(command, check=check, capture_output=capture_output)


# TODO - patch latch with proper def __repr__ -> str
def ___repr__(self):
    return str(self.local_path)


LatchFile.__repr__ = types.MethodType(___repr__, LatchFile)


def file_glob(
    pattern: str, remote_directory: str, target_dir: Optional[Path] = None
) -> List[LatchFile]:
    """Constructs a list of LatchFiles from a glob pattern.
    Convenient utility for passing collections of files between tasks. See
    [nextflow's channels](https://www.nextflow.io/docs/latest/channel.html) or
    [snakemake's wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards).
    for similar functionality in other orchestration tools.
    The remote location of each constructed LatchFile will be consructed by
    appending the file name returned by the pattern to the directory
    represented by the `remote_directory`.
    Args:
        pattern: A glob pattern to match a set of files, eg. '*.py'. Will
            resolve paths with respect to the working directory of the caller.
        remote_directory: A valid latch URL pointing to a directory, eg.
            latch:///foo. This _must_ be a directory and not a file.
        target_dir: An optional Path object to define an alternate working
            directory for path resolution
    Returns:
        A list of instantiated LatchFile objects.
    Intended Use: ::
        @small_task
        def task():
            ...
            return file_glob("*.fastq.gz", "latch:///fastqc_outputs")
    """

    if not _is_valid_url(remote_directory):
        return []

    if target_dir is None:
        wd = Path.cwd()
    else:
        wd = target_dir
    matched = sorted(wd.glob(pattern))

    return [LatchFile(str(file), remote_directory + file.name) for file in matched]


def _is_valid_url(raw_url: str) -> bool:
    """A valid URL (as a source or destination of a LatchFile) must:
    * contain a latch or s3 scheme
    * contain an absolute path
    """
    try:
        parsed = urlparse(raw_url)
    except ValueError:
        return False
    if parsed.scheme not in ("latch", "s3"):
        return False
    if not parsed.path.startswith("/"):
        return False
    return True


def _get_96_spot_pod() -> Pod:
    """[ "c6i.24xlarge", "c5.24xlarge", "c5.metal", "c5d.24xlarge", "c5d.metal" ]"""

    primary_container = V1Container(name="primary")
    resources = V1ResourceRequirements(
        requests={"cpu": "90", "memory": "170Gi"},
        limits={"cpu": "96", "memory": "192Gi"},
    )
    primary_container.resources = resources

    return Pod(
        pod_spec=V1PodSpec(
            containers=[primary_container],
            tolerations=[
                V1Toleration(effect="NoSchedule", key="ng", value="cpu-96-spot")
            ],
        ),
        primary_container_name="primary",
    )


large_spot_task = task(task_config=_get_96_spot_pod(), retries=3)


@dataclass_json
@dataclass
class SingleEndReads:
    r1: LatchFile


@dataclass_json
@dataclass
class PairedEndReads:
    r1: LatchFile
    r2: LatchFile


class ReadType(Enum):
    single = "single"
    paired = "paired"


class Strandedness(Enum):
    auto = "auto"


Replicate = Union[SingleEndReads, PairedEndReads]


@dataclass_json
@dataclass
class Sample:
    name: str
    strandedness: Strandedness
    replicates: List[Replicate]


class LatchGenome(Enum):
    RefSeq_hg38_p14 = "Homo sapiens (RefSeq hg38.p14)"
    RefSeq_T2T_CHM13v2_0 = "Homo sapiens (RefSeq T2T-CHM13v2.0)"
    RefSeq_R64 = "Saccharomyces cerevisiae (RefSeq R64)"
    RefSeq_GRCm39 = "Mus musculus (RefSeq GRCm39)"


@dataclass_json
@dataclass
class GenomeData:
    gtf: LatchFile


@dataclass_json
@dataclass
class TrimgaloreInput:
    sample_name: str
    replicate: Replicate
    replicate_index: int
    run_name: str
    base_remote_output_dir: str
    clip_r1: Optional[int]
    clip_r2: Optional[int] = None
    three_prime_clip_r1: Optional[int] = None
    three_prime_clip_r2: Optional[int] = None


@dataclass_json
@dataclass
class TrimgaloreOutput:
    sample_name: str
    trimmed_replicate: Replicate
    reports: List[LatchFile]


@dataclass_json
@dataclass
class MergedSample:
    name: str
    reads: Replicate
    strandedness: Strandedness


@dataclass_json
@dataclass
class TrimgaloreSalmonInput:
    sample_name: str
    replicates: List[Replicate]
    run_name: str
    base_remote_output_dir: str
    latch_genome: str
    bams: List[List[LatchFile]]
    custom_names: List[str]
    custom_files: List[LatchFile]
    clip_r1: Optional[int] = None
    clip_r2: Optional[int] = None
    three_prime_clip_r1: Optional[int] = None
    three_prime_clip_r2: Optional[int] = None
    # custom_gtf: Optional[LatchFile] = None
    # custom_ref_genome: Optional[LatchFile] = None
    # custom_ref_trans: Optional[LatchFile] = None
    # star_index: Optional[LatchFile] = None
    # custom_salmon_index: Optional[LatchFile] = None
    save_indices: bool = False


@dataclass_json
@dataclass
class TrimgaloreSalmonOutput:
    passed_salmon: bool
    passed_tximport: bool
    sample_name: str
    sf_files: List[LatchFile]
    auxiliary_directory: List[LatchDir]
    # trimmed_reads: List[Replicate]
    trimgalore_reports: List[LatchFile]


@small_task
def prepare_trimgalore_salmon_inputs(
    samples: List[Sample],
    run_name: str,
    latch_genome: LatchGenome,
    bams: List[List[LatchFile]],
    save_indices: bool,
    clip_r1: Optional[int] = None,
    clip_r2: Optional[int] = None,
    three_prime_clip_r1: Optional[int] = None,
    three_prime_clip_r2: Optional[int] = None,
    custom_output_dir: Optional[LatchDir] = None,
    custom_gtf: Optional[LatchFile] = None,
    custom_ref_genome: Optional[LatchFile] = None,
    custom_ref_trans: Optional[LatchFile] = None,
    custom_salmon_index: Optional[LatchFile] = None,
) -> List[TrimgaloreSalmonInput]:
    custom_names = []
    custom_files = []
    if custom_ref_genome is not None:
        custom_names.append("genome")
        custom_files.append(custom_ref_genome)
    if custom_ref_trans is not None:
        custom_names.append("trans")
        custom_files.append(custom_ref_trans)
    if custom_salmon_index is not None:
        custom_names.append("index")
        custom_files.append(custom_salmon_index)
    if custom_gtf is not None:
        custom_names.append("gtf")
        custom_files.append(custom_gtf)

    return [
        TrimgaloreSalmonInput(
            sample_name=sample.name,
            replicates=sample.replicates,
            run_name=run_name,
            clip_r1=clip_r1,
            clip_r2=clip_r2,
            three_prime_clip_r1=three_prime_clip_r1,
            three_prime_clip_r2=three_prime_clip_r2,
            base_remote_output_dir=_remote_output_dir(custom_output_dir),
            latch_genome=latch_genome.name,
            bams=bams,
            custom_names=custom_names,
            custom_files=custom_files,
            save_indices=save_indices,
        )
        for sample in samples
    ]


def _merge_replicates(
    replicates: List[Replicate], sample_name: str
) -> Union[Tuple[Path], Tuple[Path, Path]]:
    local_r1_path = f"{sample_name}_r1_merged.fq"
    r1 = _concatenate_files((str(x.r1.path) for x in replicates), local_r1_path)

    if isinstance(replicates[0], SingleEndReads):
        return (r1,)

    if not all(isinstance(x, PairedEndReads) for x in replicates):
        raise RuntimeError("Not all technical replicates were paired end")

    local_r2_path = f"{sample_name}_r2_merged.fq"
    r2 = _concatenate_files((str(x.r2.path) for x in replicates), local_r2_path)
    return (r1, r2)


def _concatenate_files(filepaths: Iterable[str], output_path: str) -> Path:
    path = Path(output_path).resolve()
    with path.open("w") as output_file:
        for p in filepaths:
            p = p.removesuffix(".gz")
            with open(p, "r") as f:
                shutil.copyfileobj(f, output_file)
    return path


def _remote_output_dir(custom_output_dir: Optional[LatchDir]) -> str:
    if custom_output_dir is None:
        return "latch:///RNA-Seq Outputs/"
    remote_path = custom_output_dir.remote_path
    assert remote_path is not None
    if remote_path[-1] != "/":
        remote_path += "/"
    return remote_path


class TrimgaloreError(Exception):
    pass


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

    local_output = f"{ts_input.sample_name}_replicate_{replicate_index}"
    returncode, stdout = _capture_output(
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

    def _output_path(middle: str) -> str:
        base = f"{ts_input.base_remote_output_dir}{ts_input.run_name}"
        tail = f"{ts_input.sample_name}/replicate_{replicate_index}/"
        return f"{base}/Quality Control Data/Trimming {middle} (TrimGalore)/{tail}"

    reads_directory = _output_path("Reads")
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

    reports_directory = _output_path("Reports")
    reports = file_glob("*trimming_report.txt", reports_directory)

    return reports, trimmed_replicate


def _build_gentrome(genome: Path, transcript: Path) -> Path:
    run(["/root/gentrome.sh", str(genome), str(transcript)])
    return Path("/root/gentrome.fa")


def _unzip_if_needed(path: Path):
    is_gzipped = str(path).endswith(".gz")

    if is_gzipped:
        # Gunzip by default deletes .gz file; no need to manually delete
        run(["gunzip", str(path)])
        unzipped_path_string = str(path).removesuffix(".gz")
        return Path(unzipped_path_string)

    return path


def _build_transcriptome(genome: Path, gtf: Path) -> Path:
    run(
        [
            "/root/RSEM-1.3.3/rsem-prepare-reference",
            "--gtf",
            str(gtf),
            "--num-threads",
            "96",
            str(genome),
            "genome",
        ]
    )
    return Path("/root/genome.transcripts.fa")


def _build_index(gentrome: Path) -> Path:
    run(
        [
            "salmon",
            "index",
            "-t",
            str(gentrome),
            "-i",
            "salmon_index",
            "--decoys",
            "decoys.txt",  # Comes from gentrome.sh
            "-k",
            "31",
            "--threads",
            "96",
        ]
    )
    return Path("/root/salmon_index")


@large_spot_task
def trimgalore_salmon(input: TrimgaloreSalmonInput) -> TrimgaloreSalmonOutput:
    outputs = [do_trimgalore(input, i, x) for i, x in enumerate(input.replicates)]
    trimmed_replicates = [x[1] for x in outputs]
    trimgalore_reports = [y for x in outputs for y in x[0]]

    gm = lgenome.GenomeManager(input.latch_genome)

    custom_salmon_index = None
    custom_ref_genome = None
    custom_ref_transcriptome = None
    custom_gtf = None
    for name, file in zip(input.custom_names, input.custom_files):
        if name == "index":
            custom_salmon_index = file
        elif name == "trans":
            custom_ref_transcriptome = file
        elif name == "genome":
            custom_ref_genome = file
        elif name == "gtf":
            custom_gtf = file

    gtf_path = None
    if custom_gtf is not None:
        gtf_path = _unzip_if_needed(custom_gtf.local_path)

    if custom_salmon_index is not None:
        # TODO: validate provided index...
        run(["tar", "-xzvf", custom_salmon_index.local_path])

        if not Path("salmon_index").is_dir():
            body = "The custom Salmon index provided must be a directory named 'salmon_index'"
            message("error", {"title": "Invalid custom Salmon index", "body": body})
            raise MalformedSalmonIndex(body)

        local_index = Path("salmon_index")

    elif custom_ref_genome is not None:
        if custom_ref_transcriptome is None and custom_gtf is None:
            body = (
                "Both a custom reference genome and GTF file need to be provided "
                "to build a local index for Salmon"
            )
            message("error", {"title": "Unable to build local index", "body": body})
            raise InsufficientCustomGenomeResources(body)

        ref_genome_path = _unzip_if_needed(custom_ref_genome.local_path)

        if custom_ref_transcriptome is None:
            if gtf_path is None:
                gtf_path = _unzip_if_needed(custom_gtf.local_path)
            ref_transcriptome = _build_transcriptome(ref_genome_path, gtf_path)
        else:
            ref_transcriptome = _unzip_if_needed(custom_ref_transcriptome.local_path)

        gentrome = _build_gentrome(ref_genome_path, ref_transcriptome)
        local_index = _build_index(gentrome)
    else:
        local_index = gm.download_salmon_index(show_progress=False)

    merged = _merge_replicates(trimmed_replicates, input.sample_name)
    if len(merged) == 1:
        (r1,) = merged
        reads = ["-r", str(r1)]
    else:
        r1, r2 = merged
        reads = ["-1", str(r1), "-2", str(r2)]

    # Delete unneeded files to free space
    for rep in trimmed_replicates:
        os.remove(rep.r1.path)
        if isinstance(rep, PairedEndReads):
            os.remove(rep.r2.path)

    returncode, stdout = _capture_output(
        [
            "salmon",
            "quant",
            "-i",
            str(local_index),
            "-l",
            "A",
            *reads,
            "--threads",
            str(96),
            "--validateMappings",
            "-o",
            "salmon_quant",
        ]
    )

    # Delete unneeded files to free space
    for path in merged:
        os.remove(path)

    # todo(rohankan): Add info messages too? Decide which Salmon logs
    # are interesting or important enough to show on console
    identifier = f"sample {input.sample_name}"
    errors = []
    for alert_type, alert_message in re.findall(_SALMON_ALERT_PATTERN, stdout):
        data = {"title": f"Salmon {alert_type} for {identifier}", "body": alert_message}
        if alert_type == "warning":
            message("warning", data)
        else:
            message("error", data)
            errors.append(alert_message)

    if returncode != 0:
        deets = "\n".join(["Error(s) occurred while running Salmon", *errors])
        raise RuntimeError(deets)

    salmon_quant_file = LatchFile(
        "/root/salmon_quant/quant.sf",
        _salmon_output_path(input, f"{input.sample_name}_quant.sf"),
    )
    sf_files = [salmon_quant_file]

    try:
        gtf_path = (
            custom_gtf.local_path
            if custom_gtf is not None
            else gm.download_gtf(show_progress=False)
        )

        tximport_output_path = "/root/salmon_quant/genome_abundance.sf"
        subprocess.run(
            [
                "/root/wf/run_tximport.R",
                "--args",
                str(salmon_quant_file.path),
                gtf_path,
                tximport_output_path,
            ],
            capture_output=True,
        )

        sf_files.append(
            LatchFile(
                tximport_output_path,
                _salmon_output_path(input, f"{input.sample_name}_genome_abundance.sf"),
            )
        )
    except subprocess.CalledProcessError as e:
        message("error", {"title": f"tximport error for {identifier}", "body": str(e)})
        print(
            f"Unable to produce gene mapping from tximport. Error surfaced from tximport -> {e}"
        )
        raise RuntimeError(f"Tximport error: {e}")

    auxiliary_directory = LatchDir(
        "/root/salmon_quant/aux_info",
        _salmon_output_path(input, "Auxiliary Info"),
    )

    return TrimgaloreSalmonOutput(
        passed_salmon=True,
        passed_tximport=True,
        sample_name=input.sample_name,
        sf_files=sf_files,
        auxiliary_directory=[auxiliary_directory],
        # trimmed_reads=trimmed_replicates,
        trimgalore_reports=trimgalore_reports,
    )


class InsufficientCustomGenomeResources(Exception):
    pass


class MalformedSalmonIndex(Exception):
    pass


class SalmonError(Exception):
    pass


# Each Salmon warning or error log starts with a timestamp surrounded in square
# brackets ('\d{4}' represents the first part of the timestamp - the year)
_SALMON_ALERT_PATTERN = re.compile(r"\[(warning|error)\] (.+?)(?:\[\d{4}|$)", re.DOTALL)


def _salmon_output_path(inp: TrimgaloreSalmonInput, suffix: str) -> str:
    return f"{inp.base_remote_output_dir}{inp.run_name}/Quantification (salmon)/{inp.sample_name}/{suffix}"


_COUNT_TABLE_GENE_ID_COLUMN = "gene_id"


@small_task
def count_matrix_and_multiqc(
    run_name: str,
    outputs: List[TrimgaloreSalmonOutput],
    output_directory: Optional[LatchDir],
    latch_genome: LatchGenome,
    custom_gtf: Optional[LatchFile] = None,
) -> List[LatchFile]:
    output_files = []

    # Beginning creation of combined count matrix
    tximport_outputs = [x for x in outputs if x.passed_tximport]
    if len(tximport_outputs) == len(outputs):
        message(
            "info",
            {
                "title": "Generating count matrix from all samples",
                "body": "\n".join(f"- {x.sample_name}" for x in tximport_outputs),
            },
        )

        combined_counts = defaultdict(dict)
        for output in tximport_outputs:
            genome_abundance_file = next(
                x for x in output.sf_files if x.remote_path.endswith("dance.sf")
            )
            with open(genome_abundance_file.local_path, "r") as f:
                for row in csv.DictReader(f, dialect=csv.excel_tab):
                    gene_name = row["Name"]
                    combined_counts[gene_name][output.sample_name] = row["NumReads"]

        raw_count_table_path = Path("./counts.tsv").resolve()
        with raw_count_table_path.open("w") as file:
            sample_names = (x.sample_name for x in tximport_outputs)
            writer = csv.DictWriter(
                file,
                [_COUNT_TABLE_GENE_ID_COLUMN, *sample_names],
                delimiter="\t",
            )
            writer.writeheader()
            for gene_id, data in combined_counts.items():
                data[_COUNT_TABLE_GENE_ID_COLUMN] = gene_id
                writer.writerow(data)

        base = _remote_output_dir(output_directory) + run_name
        output_path = f"{base}/Quantification (salmon)/counts.tsv"
        count_matrix_file = LatchFile(str(raw_count_table_path), output_path)
        output_files.append(count_matrix_file)
    else:
        message(
            "warning",
            {
                "title": "Unable to create combined count matrix",
                "body": "Some samples failed in the earlier step",
            },
        )

    # Beginning creation of MultiQC report
    salmon_outputs = [x for x in outputs if x.passed_salmon]
    paths = [Path(x.auxiliary_directory[0].local_path) for x in salmon_outputs]
    for output in salmon_outputs:
        paths += [Path(x.local_path) for x in output.sf_files]

    try:
        subprocess.run(["multiqc", *paths], check=True)
        base = _remote_output_dir(output_directory) + run_name
        multiqc_report_file = LatchFile(
            "/root/multiqc_report.html",
            f"{base}/multiqc_report.html",
        )
        output_files.append(multiqc_report_file)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while generating MultiQC report -> {e}")
        message(
            "error",
            {
                "title": "Unable to generate MultiQC report",
                "body": "See logs for more information",
            },
        )

    return output_files


class AlignmentTools(Enum):
    star_salmon = "Traditional Alignment + Quantification"
    salmon = "Selective Alignment + Quantification"


@workflow
def rnaseq(
    samples: List[Sample],
    alignment_quantification_tools: AlignmentTools,
    ta_ref_genome_fork: str,
    sa_ref_genome_fork: str,
    output_location_fork: str,
    run_name: str,
    latch_genome: LatchGenome,
    bams: List[List[LatchFile]],
    custom_gtf: Optional[LatchFile] = None,
    custom_ref_genome: Optional[LatchFile] = None,
    custom_ref_trans: Optional[LatchFile] = None,
    star_index: Optional[LatchFile] = None,
    salmon_index: Optional[LatchFile] = None,
    save_indices: bool = False,
    custom_output_dir: Optional[LatchDir] = None,
) -> List[LatchFile]:
    """Perform alignment and quantification on Bulk RNA-Sequencing reads

    Bulk RNA-Seq (Alignment and Quantification)
    ----

    This workflow will produce gene and transcript counts from bulk RNA-seq
    sample reads.

    ## Workflow Anatomy

    #### Disclaimer

    This workflow assumes that your sequencing reads were derived from *short-read
    cDNA seqeuncing* ( as opposed to long-read cDNA/direct RNA sequencing). If in
    doubt, you can likely make the same assumption, as it is by far the most common
    form of "RNA-sequencing".

    ### Brief Summary of RNA-seq

    This workflow ingests short-read sequencing files (in FastQ format) that came
    from the following sequence of steps[^1]:

      - RNA extraction from sample
      - cDNA synthesis from extracted RNA
      - adaptor ligation / library prep
      - (likely) PCR amplification of library
      - sequencing of library

    You will likely end up with one or more FastQ files from this process that hold
    the sequencing reads in raw text form. This will be the starting point of our
    workflow.

    (If you have a `.bcl` file, this holds the raw output of a sequencing machine.
    There are there are [external
    tools](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)
    that can convert these files to FastQ format, which you will need before you can
    proceed).

    ### Quality Control

    As a pre-processing step, its important to check the quality of your sequencing
    files. FastQC is the industry staple for generating a report of useful summary
    statistics[^2] and is available if you double-click on a file on the [LatchBio
    platform](https://console.latch.bio).

    The following are the most useful of these statistics:

      - *Per base sequence quality* gives the per-site distribution over the length
    of the read
      - *Sequence duplication levels* reveals duplicated reads, indicating degraded
    RNA samples or aggressive PCR cycling[^1]

    For a full breakdown of the values and their interpretation, we refer the
    reader to this
    [tutorial](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html).

    #### Trimming

    Short-read sequencing introduces adapters, small sequences attached to the 5'
    and 3' end of cDNA fragments, that are present as artifacts in our FastQ files
    and must be removed.

    We have yet to identify a comprehensive review of the various trimming tools, so
    we have selected [TrimGalore](https://github.com/FelixKrueger/TrimGalore)
    trusted by researchers we work with out of UCSF and Stanford, until we are able
    to do so ourself.

    ### Alignment

    Alignment is the process of assigning a sequencing read a location on a
    reference genome or transcriptome. It is the most computationally expensive step
    of the workflow, requiring a comparison against the entire reference sequence
    for each of millions of reads.

    Transcript alignment was initially conducted similarly to genomic alignment,
    using tools like
    [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to rigorously
    recover reference coordinates for each read. This was eschewed for a lighter
    "pseudo-alignment" in the years that followed that assigned each read to a
    transcript rather than an exact location, saving time and resources. However,
    while these methods are faster, they have proven to be less accurate.[^3]

    In 2020, the **Selective Alignment** algorithm was introduced that performed a
    similar lightweight read assignment while simultaneously *outperforming*
    traditional alignment methods in accuracy.[^3] We utilize
    [salmon](https://github.com/COMBINE-lab/salmon) to implement selective
    alignment.

    ### Gene Count Quantification

    Selective Alignment produces estimations of transcript abundances. Recall that
    that there can be multiple transcripts for any single gene. It is desirable to
    have estimated gene counts for two reasons:

      1. gene counts are a more stable measure of transcription.\*
      2. gene counts are more interpretable

    \* Stability is loosely defined as consistent correlation with ground truth
    counts as the available (transcript) annotations begin to drop out. [^4]

    We utilize [tximport](ps://bioconductor.org/packages/release/bioc/html/tximport.html) to
    perform the conversion of transcripts to read counts.

    [^1]: Stark, Rory; Grzelak, Marta; Hadfield, James (2019). RNA sequencing: the teenage years. Nature Reviews Genetics, (), –. doi:10.1038/s41576-019-0150-2
    [^2]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    [^3]: Srivastava, A., Malik, L., Sarkar, H. et al. Alignment and mapping methodology influence transcript abundance estimation. Genome Biol 21, 239 (2020). https://doi.org/10.1186/s13059-020-02151-8
    [^4]: Soneson C, Love MI and Robinson MD. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences [version 1; peer review: 2 approved]. F1000Research 2015, 4:1521


    __metadata__:
        display_name: Bulk RNAseq
        author:
            name: LatchBio
            email:
            github:
        repository: github.com/latchbio/rnaseq
        license:
            id: MIT
        flow:
        - section: Samples
          flow:
            - text: >-
                  Sample files can be provided and their read type can be
                  inferred from their name or this information can be specified manually.
                  Sample strandedness is inferred automatically (learn more).

            - params:
                - samples
        - section: Alignment and Quantification
          flow:
            - text: >-
                  This workflow uses Salmon's selective alignment described in this
                  [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8),
                  which achieves greater accuracy than traditional alignment methods while
                  using less computational resources.

            - fork: alignment_quantification_tools
              flows:
                traditional:
                    display_name: Selective Alignment
                    flow:
                        - fork: ta_ref_genome_fork
                          flows:
                            database:
                                display_name: Select from Latch Genome Database
                                flow:
                                    - text: >-
                                        We have curated a set of reference
                                        genome data for ease and
                                        reproducibility. More information about
                                        these managed files can be found
                                        [here](https://github.com/latchbio/latch-genomes).
                                    - params:
                                        - latch_genome
                            custom:
                                display_name: Provide Custom Genome
                                _tmp_unwrap_optionals:
                                    - custom_gtf
                                    - custom_ref_genome
                                flow:
                                    - params:
                                        - custom_ref_genome
                                        - custom_gtf
                                    - spoiler: Optional Params
                                      flow:
                                        - text: >-
                                            These files will be generated from the
                                            GTF/Genome files if not provided.
                                        - params:
                                            - salmon_index
                                            - custom_ref_trans
        - section: Output Location
          flow:
          - params:
              - run_name
          - fork: output_location_fork
            flows:
                default:
                    display_name: Default
                    flow:
                    - text:
                        Output will be at default location in the data
                        viewer - RNA-Seq Outputs/"Run Name"
                custom:
                    display_name: Specify Custom Path
                    _tmp_unwrap_optionals:
                        - custom_output_dir
                    flow:
                    - params:
                        - custom_output_dir
    Args:

        samples:
            Here you can organize your FastQ files by sample and add technical
            replicates for each sample.  Biological replicates should be
            organized as separate samples.

          __metadata__:
            display_name: Sample Sheet
            batch_table_column: true
            _tmp:
                custom_ingestion: auto

        alignment_quantification_tools:

          __metadata__:
            display_name: Alignment and Quantification Method

        latch_genome:
          Curated reference files for specific genome sources and builds.

          __metadata__:
            batch_table_column: true
            display_name: Genome Database Option

        sa_ref_genome_fork:
          Select a reference genome from our curated database or provide your own.

          __metadata__:
            display_name: Reference Genome Source

        ta_ref_genome_fork:
          Select a reference genome from our curated database or provide your own.

          __metadata__:
            display_name: Reference Genome Source

        custom_ref_genome:
          The reference genome you want to align you samples to.

          __metadata__:
            display_name: Reference Genome File
            appearance:
                detail: (.fasta, .fasta.gz, .fa, .fa.gz, .fna, .fna.gz)

        custom_gtf:
          The gene annonation file that corresponds to the reference genome
          provided.

          __metadata__:
            display_name: Annotation File
            appearance:
                detail: (.gtf)

        bams:
          foobar

          __metadata__:
            display_name: bams

        custom_ref_trans:
          If not provided the workflow will generate from the Annotation File
          and Reference Genome File.

          __metadata__:
            display_name: Reference Transcript File (optional)
            appearance:
                detail: (.fasta, .fasta.gz, .fa, .fa.gz, .fna, .fna.gz)

        star_index:
          You are able to provide a zipped prebuilt STAR alignment index for
          your genome. This will speed up run time as the index is generated if
          none is provided. In output settings you are able to save indices
          from a run to be used in future runs.

          __metadata__:
            display_name: Provide Prebuilt STAR Index

        salmon_index:
            You are able to provide a zipped prebuilt Salmon selective alignment
            index for your genome. This will speed up run time as the index is
            generated if none is provided.

          __metadata__:
            display_name: salmon Index
            appearance:
                detail: (.tar.gz is only accepted extension)

        save_indices:
            If you provided a custom genome you can output the alignment
            indexes generated from this run for use in future runs. This will
            speed up runtime since the workflow doesn't have to then regenerate
            the indexes.

          __metadata__:
            display_name: Save Generated Reference Indexes

        run_name:
          A name for this analysis run, this will be used to name outputs from
          this run.

          __metadata__:
            batch_table_column: true
            display_name: Run Name

        output_location_fork:

        custom_output_dir:
          You can provide a custom location where this run's analysis outputs
          will be located.

          __metadata__:
            display_name: Custom Output Location
    """
    inputs = prepare_trimgalore_salmon_inputs(
        samples=samples,
        run_name=run_name,
        clip_r1=None,
        clip_r2=None,
        three_prime_clip_r1=None,
        three_prime_clip_r2=None,
        custom_output_dir=custom_output_dir,
        latch_genome=latch_genome,
        bams=bams,
        custom_gtf=custom_gtf,
        custom_ref_genome=custom_ref_genome,
        custom_ref_trans=custom_ref_trans,
        custom_salmon_index=salmon_index,
        save_indices=save_indices,
    )
    outputs = map_task(trimgalore_salmon)(input=inputs)
    return count_matrix_and_multiqc(
        run_name=run_name,
        outputs=outputs,
        output_directory=custom_output_dir,
        latch_genome=latch_genome,
        custom_gtf=custom_gtf,
    )


if __name__ == "wf":
    LaunchPlan.create(
        "wf.__init__.rnaseq.Test Data",
        rnaseq,
        default_inputs={
            "samples": [
                Sample(
                    name="test_sample",
                    strandedness=Strandedness.auto,
                    replicates=[
                        PairedEndReads(
                            r1=LatchFile(
                                "s3://latch-public/welcome/rnaseq/r1.fastq",
                            ),
                            r2=LatchFile(
                                "s3://latch-public/welcome/rnaseq/r2.fastq",
                            ),
                        )
                    ],
                )
            ],
            "run_name": "Test Run",
        },
    )
