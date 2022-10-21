from typing import Annotated

from flytekit.core.annotation import FlyteAnnotation
from latch.types import LatchFile, LatchRule
from typing_extensions import TypeAlias

FastaFile: TypeAlias = Annotated[
    LatchFile,
    FlyteAnnotation(
        {
            "rules": [
                LatchRule(
                    message="Provided path is not a valid .fasta file",
                    regex="^.*\.(fa|fa.gz|fasta|fasta.gz|fq|fq.gz)$",
                ).dict
            ]
        },
    ),
]

GtfFile: TypeAlias = Annotated[
    LatchFile,
    FlyteAnnotation(
        {
            "rules": [
                LatchRule(
                    message="Provided path is not a valid .gtf file",
                    regex="^.*\.(gtf|gtf.gz)$",
                ).dict
            ]
        },
    ),
]
