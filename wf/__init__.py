"""latch/rnaseq"""

import functools
from typing import Annotated, List, Optional

from flytekit.core.annotation import FlyteAnnotation
from latch import map_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile

from wf.core.models import (AlignmentTools, LatchGenome, PairedEndReads,
                            Sample, SingleEndReads, Strandedness)
from wf.subworkflows.deseq2 import deseq2_wf
from wf.tasks.count_matrix_and_multiqc import count_matrix_and_multiqc
from wf.tasks.leafcutter import leafcutter
from wf.tasks.prepare_inputs import prepare_inputs
from wf.tasks.trimgalore_salmon import trimgalore_salmon

print = functools.partial(print, flush=True)


@workflow
def rnaseq(
    samples: List[Sample],
    alignment_quantification_tools: AlignmentTools,
    ta_ref_genome_fork: str,
    sa_ref_genome_fork: str,
    output_location_fork: str,
    run_name: str,
    latch_genome: LatchGenome,
    conditions_source: str = "none",
    manual_conditions: Annotated[
        List[List[str]],
        FlyteAnnotation({"_tmp_hack_deseq2": "manual_design_matrix"}),
    ] = [],
    conditions_table: Optional[
        Annotated[
            LatchFile,
            FlyteAnnotation(
                {
                    "_tmp_hack_deseq2": "design_matrix",
                    "rules": [
                        {
                            "regex": r".*\.(csv|tsv|xlsx)$",
                            "message": "Expected a CSV, TSV, or XLSX file",
                        }
                    ],
                }
            ),
        ]
    ] = None,
    design_matrix_sample_id_column: Optional[
        Annotated[str, FlyteAnnotation({"_tmp_hack_deseq2": "design_id_column"})]
    ] = None,
    design_formula: Annotated[
        List[List[str]],
        FlyteAnnotation(
            {
                "_tmp_hack_deseq2": "design_formula",
                "_tmp_hack_deseq2_allow_clustering": True,
            }
        ),
    ] = [],
    custom_gtf: Optional[LatchFile] = None,
    custom_ref_genome: Optional[LatchFile] = None,
    custom_ref_trans: Optional[LatchFile] = None,
    star_index: Optional[LatchFile] = None,
    salmon_index: Optional[LatchFile] = None,
    run_splicing: bool = False,
    custom_output_dir: Optional[LatchDir] = None,
):
    ts_inputs = prepare_inputs(
        samples=samples,
        run_name=run_name,
        custom_output_dir=custom_output_dir,
        latch_genome=latch_genome,
        custom_gtf=custom_gtf,
        custom_ref_genome=custom_ref_genome,
        custom_ref_trans=custom_ref_trans,
        custom_salmon_index=salmon_index,
        run_splicing=run_splicing,
    )
    outputs = map_task(trimgalore_salmon)(ts_input=ts_inputs)
    count_matrix_file, multiqc_report_file = count_matrix_and_multiqc(
        run_name=run_name,
        ts_outputs=outputs,
        output_directory=custom_output_dir,
    )
    a, b, c = leafcutter(
        run_splicing=run_splicing,
        run_name=run_name,
        ts_outputs=outputs,
        output_directory=custom_output_dir,
        manual_conditions=manual_conditions,
    )
    deseq2_wf(
        report_name=run_name,
        count_table_source="single",
        raw_count_table=count_matrix_file,
        raw_count_tables=[],
        count_table_gene_id_column="gene_id",
        output_location_type="default",
        output_location=custom_output_dir,
        conditions_source=conditions_source,
        manual_conditions=manual_conditions,
        conditions_table=conditions_table,
        design_matrix_sample_id_column=design_matrix_sample_id_column,
        design_formula=design_formula,
        number_of_genes_to_plot=30,
    )


LaunchPlan(
    rnaseq,
    "Small Data - Human 10K Reads",
    {
        "samples": [
            Sample(
                name="Test",
                strandedness=Strandedness.auto,
                replicates=[
                    SingleEndReads(
                        r1=LatchFile(
                            "s3://latch-public/verified/bulk-rnaseq/small-test/10k_reads_human.fastq.gz",
                        ),
                    ),
                ],
            ),
        ],
        "run_name": "Small Test",
    },
)

LaunchPlan(
    rnaseq,
    "Small Data - Yeast Chr I",
    {
        "samples": [
            Sample(
                name="chrI",
                strandedness=Strandedness.auto,
                replicates=[
                    PairedEndReads(
                        r1=LatchFile("s3://latch-public/test-data/1/gene.r1.fq.gz"),
                        r2=LatchFile("s3://latch-public/test-data/1/gene.r2.fq.gz"),
                    ),
                ],
            ),
        ],
        "run_name": "yeast chrI",
    },
)

LaunchPlan(
    rnaseq,
    "Test Data - CoCl2 vs Control (Knyazev, 2021)",
    {
        "samples": [
            Sample(
                name="Control rep 1",
                strandedness=Strandedness.auto,
                replicates=[
                    SingleEndReads(
                        r1=LatchFile(
                            "s3://latch-public/verified/bulk-rnaseq/IBD/Control_rep1.fastq.gz",
                        ),
                    ),
                ],
            ),
            Sample(
                name="Control rep 2",
                strandedness=Strandedness.auto,
                replicates=[
                    SingleEndReads(
                        r1=LatchFile(
                            "s3://latch-public/verified/bulk-rnaseq/IBD/Control_rep2.fastq.gz",
                        ),
                    ),
                ],
            ),
            Sample(
                name="Control rep 3",
                strandedness=Strandedness.auto,
                replicates=[
                    SingleEndReads(
                        r1=LatchFile(
                            "s3://latch-public/verified/bulk-rnaseq/IBD/Control_rep3.fastq.gz",
                        ),
                    ),
                ],
            ),
            Sample(
                name="CoCl2 rep 1",
                strandedness=Strandedness.auto,
                replicates=[
                    SingleEndReads(
                        r1=LatchFile(
                            "s3://latch-public/verified/bulk-rnaseq/IBD/CoCl2_rep1.fastq.gz",
                        ),
                    ),
                ],
            ),
            Sample(
                name="CoCl2 rep 2",
                strandedness=Strandedness.auto,
                replicates=[
                    SingleEndReads(
                        r1=LatchFile(
                            "s3://latch-public/verified/bulk-rnaseq/IBD/CoCl2_rep2.fastq.gz",
                        ),
                    ),
                ],
            ),
            Sample(
                name="Oxy rep 1",
                strandedness=Strandedness.auto,
                replicates=[
                    SingleEndReads(
                        r1=LatchFile(
                            "s3://latch-public/verified/bulk-rnaseq/IBD/Oxy_rep1.fastq.gz",
                        ),
                    ),
                ],
            ),
            Sample(
                name="Oxy rep 2",
                strandedness=Strandedness.auto,
                replicates=[
                    SingleEndReads(
                        r1=LatchFile(
                            "s3://latch-public/verified/bulk-rnaseq/IBD/Oxy_rep2.fastq.gz",
                        ),
                    ),
                ],
            ),
            Sample(
                name="Oxy rep 3",
                strandedness=Strandedness.auto,
                replicates=[
                    SingleEndReads(
                        r1=LatchFile(
                            "s3://latch-public/verified/bulk-rnaseq/IBD/Oxy_rep3.fastq.gz",
                        ),
                    ),
                ],
            ),
        ],
        "run_name": "Inflammatory Bowel Disease",
    },
)
