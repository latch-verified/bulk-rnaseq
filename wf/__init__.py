"""latch/rnaseq"""

import functools
from typing import Annotated, List, Optional

from flytekit.core.annotation import FlyteAnnotation
from latch import map_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile

from wf.core.metadata import metadata
from wf.core.models import (AlignmentTools, LatchGenome, PairedEndReads,
                            Sample, SingleEndReads, Strandedness)
from wf.subworkflows.deseq2 import deseq2_wf
from wf.tasks.count_matrix_and_multiqc import count_matrix_and_multiqc
from wf.tasks.leafcutter import leafcutter
from wf.tasks.prepare_inputs import prepare_inputs
from wf.tasks.trimgalore_salmon import trimgalore_salmon

print = functools.partial(print, flush=True)


@workflow(metadata)
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
    """Perform alignment and quantification on Bulk RNA-Sequencing reads

        Bulk RNA-Seq (Alignment and Quantification)
        ----

        This workflow will produce gene and transcript counts from bulk RNA-seq
        sample reads.

        # Workflow Anatomy

        # Disclaimer

        This workflow assumes that your sequencing reads were derived from *short-read
        cDNA seqeuncing* ( as opposed to long-read cDNA/direct RNA sequencing). If in
        doubt, you can likely make the same assumption, as it is by far the most common
        form of "RNA-sequencing".

        # Brief Summary of RNA-seq

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

        # Quality Control

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

        # Trimming

        Short-read sequencing introduces adapters, small sequences attached to the 5'
        and 3' end of cDNA fragments, that are present as artifacts in our FastQ files
        and must be removed.

        We have yet to identify a comprehensive review of the various trimming tools, so
        we have selected [TrimGalore](https://github.com/FelixKrueger/TrimGalore)
        trusted by researchers we work with out of UCSF and Stanford, until we are able
        to do so ourself.

        # Alignment

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

        # Gene Count Quantification

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
    """


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
