from pathlib import Path
from typing import List, Optional

import lgenome
from latch import large_task, message
from latch.types import LatchDir, LatchFile

from wf.errors import InsufficientCustomGenomeResources, MalformedSalmonIndex
from wf.models import LatchGenome, Sample, TrimgaloreSalmonInput
from wf.types import FastaFile, GtfFile
from wf.utils import remote_output_dir, run, unzip_if_needed

prepare_inputs_dockerfile = (
    Path(__file__).parent.parent / "dockerfiles/Dockerfile.prepare_inputs"
)


@large_task(dockerfile=prepare_inputs_dockerfile)
def prepare_inputs(
    samples: List[Sample],
    run_name: str,
    latch_genome: LatchGenome,
    custom_output_dir: Optional[LatchDir] = None,
    custom_gtf: Optional[LatchFile] = None,
    custom_ref_genome: Optional[LatchFile] = None,
    custom_ref_trans: Optional[LatchFile] = None,
    custom_salmon_index: Optional[LatchFile] = None,
    run_splicing: bool = False,
) -> List[TrimgaloreSalmonInput]:
    """TODO"""

    gm = lgenome.GenomeManager(latch_genome.name)

    gtf_path = (
        unzip_if_needed(custom_gtf.local_path)
        if custom_gtf is not None
        else gm.download_gtf(show_progress=False)
    )

    if custom_salmon_index is not None:

        run(["tar", "-xzvf", custom_salmon_index.local_path])

        if not Path("salmon_index").is_dir():
            body = "The custom Salmon index provided must be a directory named 'salmon_index'"
            message("error", {"title": "Invalid custom Salmon index", "body": body})
            raise MalformedSalmonIndex(body)

        salmon_index_path = Path("salmon_index")

    elif custom_ref_genome is not None:
        if custom_ref_trans is None and custom_gtf is None:
            body = (
                "Both a custom reference genome and GTF file need to be provided "
                "to build a local index for Salmon"
            )
            message("error", {"title": "Unable to build local index", "body": body})
            raise InsufficientCustomGenomeResources(body)

        ref_genome = unzip_if_needed(custom_ref_genome.local_path)
        if custom_ref_trans is None:
            ref_trans = rsem_prepare_reference(ref_genome, gtf_path)
        else:
            ref_trans = unzip_if_needed(custom_ref_trans.local_path)

        gentrome = build_gentrome(ref_genome, ref_trans)
        salmon_index_path = build_salmon_index(gentrome)
    else:
        salmon_index_path = gm.download_salmon_index(show_progress=False)

    return [
        TrimgaloreSalmonInput(
            sample_name=sample.name,
            replicates=sample.replicates,
            run_name=run_name,
            base_remote_output_dir=remote_output_dir(custom_output_dir),
            gtf=LatchFile(gtf_path),
            salmon_index=LatchDir(salmon_index_path),
            run_splicing=run_splicing,
        )
        for sample in samples
    ]


def rsem_prepare_reference(genome: FastaFile, gtf: GtfFile) -> Path:
    """TODO"""

    genome = unzip_if_needed(genome)
    unzip_if_needed(gtf)

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


def build_gentrome(genome: FastaFile, transcriptome: FastaFile) -> Path:
    """TODO"""

    run(["/root/gentrome.sh", str(genome), str(transcriptome)])

    return Path("/root/gentrome.fa")


def build_salmon_index(gentrome: Path) -> Path:
    """TODO"""

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
