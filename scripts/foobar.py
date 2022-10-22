import subprocess

from latch.types import LatchFile
from wf.core.models import (
    LatchGenome,
    PairedEndReads,
    Sample,
    SingleEndReads,
    Strandedness,
)
from wf.tasks.prepare_inputs.prepare_inputs import prepare_inputs

zipped_paired_sample = Sample(
    name="gene",
    strandedness=Strandedness.auto,
    replicates=[
        PairedEndReads(
            r1=LatchFile("s3://latch-public/test-data/4107/gene.r1.fq.gz"),
            r2=LatchFile("s3://latch-public/test-data/4107/gene.r2.fq.gz"),
        ),
    ],
)

zipped_single_sample = Sample(
    name="gene",
    strandedness=Strandedness.auto,
    replicates=[
        SingleEndReads(
            r1=LatchFile("s3://latch-public/test-data/4107/gene.r2.fq.gz"),
        ),
    ],
)

# import time
# time.sleep(10000)

# for s in (zipped_paired_sample, zipped_single_sample):
#     for gtf in (
#         LatchFile("s3://latch-public/test-data/1/gene.gtf.gz"),
#         LatchFile("s3://latch-public/test-data/1/gene.gtf"),
#     ):
#         for fa in (
#             LatchFile("s3://latch-public/test-data/1/gene.fa"),
#             LatchFile("s3://latch-public/test-data/1/gene.fa.gz"),
#         ):
#
#             import time
#
#             time.sleep(1000000)
#             prepare_inputs(
#                 samples=[s],
#                 run_name="test",
#                 latch_genome=LatchGenome.RefSeq_hg38_p14,
#                 custom_output_dir=None,
#                 custom_gtf=gtf,
#                 custom_ref_genome=fa,
#                 custom_ref_trans=None,
#                 custom_salmon_index=None,
#                 run_splicing=False,
#             )

print("before")

print(
    prepare_inputs(
        samples=[zipped_paired_sample],
        run_name="test",
        latch_genome=LatchGenome.RefSeq_hg38_p14,
        custom_output_dir=None,
        custom_gtf=LatchFile("s3://latch-public/test-data/4107/gene.gtf.gz"),
        custom_ref_genome=LatchFile("s3://latch-public/test-data/4107/gene.fa.gz"),
        custom_ref_trans=None,
        custom_salmon_index=None,
        run_splicing=False,
    )
)

print("after")

# prepare_inputs(
#     samples=[],
#     run_name="test",
#     latch_genome=LatchGenome.RefSeq_hg38_p14,
#     custom_output_dir=None,
#     custom_gtf=LatchFile("s3://latch-public/test-data/1/gene.gtf"),
#     custom_ref_genome=LatchFile("s3://latch-public/test-data/1/gene.fa.gz"),
#     custom_ref_trans=None,
#     custom_salmon_index=None,
#     run_splicing=False,
# )
#
# prepare_inputs(
#     samples=[],
#     run_name="test",
#     latch_genome=LatchGenome.RefSeq_hg38_p14,
#     custom_output_dir=None,
#     custom_gtf=LatchFile("s3://latch-public/test-data/1/gene.gtf.gz"),
#     custom_ref_genome=LatchFile("s3://latch-public/test-data/1/gene.fa"),
#     custom_ref_trans=None,
#     custom_salmon_index=None,
#     run_splicing=False,
# )
#
# prepare_inputs(
#     samples=[],
#     run_name="test",
#     latch_genome=LatchGenome.RefSeq_hg38_p14,
#     custom_output_dir=None,
#     custom_gtf=LatchFile("s3://latch-public/test-data/1/gene.gtf"),
#     custom_ref_genome=LatchFile("s3://latch-public/test-data/1/gene.fa"),
#     custom_ref_trans=None,
#     custom_salmon_index=None,
#     run_splicing=False,
# )
#
#
# prepare_inputs(
#     samples=[],
#     run_name="test",
#     latch_genome=LatchGenome.RefSeq_hg38_p14,
#     custom_output_dir=None,
#     custom_gtf=LatchFile("s3://latch-public/test-data/1/gene.gtf"),
#     custom_ref_genome=LatchFile("s3://latch-public/test-data/1/gene.fa"),
#     custom_ref_trans=None,
#     custom_salmon_index=None,
#     run_splicing=False,
# )
