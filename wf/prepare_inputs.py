from pathlib import Path
from typing import List, Optional

import lgenome
from latch import large_task, message
from latch.types import LatchDir, LatchFile

from .errors import InsufficientCustomGenomeResources, MalformedSalmonIndex
from .models import LatchGenome, Sample, TrimgaloreSalmonInput
from .utils import run


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


def _build_gentrome(genome: Path, transcript: Path) -> Path:
    run(["/root/gentrome.sh", str(genome), str(transcript)])
    return Path("/root/gentrome.fa")


@large_task
def prepare_inputs(
    samples: List[Sample],
    run_name: str,
    latch_genome: LatchGenome,
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
    run_splicing: bool = False,
) -> List[TrimgaloreSalmonInput]:
    """Prepare all reference files one time.

    This includes building custom indices or downloading managed files.
    """

    gm = lgenome.GenomeManager(latch_genome)

    gtf_path = (
        _unzip_if_needed(custom_gtf.local_path)
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

        ref_genome_path = _unzip_if_needed(custom_ref_genome.local_path)

        if custom_ref_trans is None:
            ref_trans = _build_transcriptome(ref_genome_path, gtf_path)
        else:
            ref_trans = _unzip_if_needed(custom_ref_trans.local_path)

        gentrome = _build_gentrome(ref_genome_path, ref_trans)
        salmon_index_path = _build_index(gentrome)
    else:
        salmon_index_path = gm.download_salmon_index(show_progress=False)

    return [
        TrimgaloreSalmonInput(
            sample_name=sample.name,
            replicates=sample.replicates,
            run_name=run_name,
            base_remote_output_dir=_remote_output_dir(custom_output_dir),
            clip_r1=clip_r1,
            clip_r2=clip_r2,
            three_prime_clip_r1=three_prime_clip_r1,
            three_prime_clip_r2=three_prime_clip_r2,
            gtf=LatchFile(gtf_path),
            salmon_index=LatchDir(salmon_index_path),
            save_indices=save_indices,
            run_splicing=run_splicing,
        )
        for sample in samples
    ]
