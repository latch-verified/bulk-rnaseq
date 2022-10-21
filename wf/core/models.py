from dataclasses import dataclass
from enum import Enum
from typing import List, Union

from dataclasses_json import dataclass_json
from latch.types import LatchDir, LatchFile


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
    gtf: LatchFile
    salmon_index: LatchDir
    run_splicing: bool = False


@dataclass_json
@dataclass
class TrimgaloreSalmonOutput:
    passed_salmon: bool
    passed_tximport: bool
    sample_name: str
    salmon_aux_output: LatchDir
    salmon_quant_file: LatchFile
    """Tab-separated transcript quantification file."""
    gene_abundance_file: LatchFile
    """Gene abundance file from tximport."""
    trimgalore_reports: List[LatchFile]
    """Temporary field for junction file."""
    junction_file: LatchFile


@dataclass_json
@dataclass
class TrimgaloreOutput:
    sample_name: str
    trimmed_replicate: Replicate
    reports: List[LatchFile]


class AlignmentTools(Enum):
    star_salmon = "Traditional Alignment + Quantification"
    salmon = "Selective Alignment + Quantification"
