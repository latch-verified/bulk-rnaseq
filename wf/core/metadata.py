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


    __metadata__:
        display_name: Bulk RNAseq
        wiki_url: https://www.latch.wiki/bulk-rna-seq-end-to-end
        video_tutorial: https://www.loom.com/share/dfba09ba6f524722b5d829f2424a3a3f
        author:
            name: LatchBio
            email: help@latch.bio
            github: github.com/latchbio
        repository: github.com/latch-verified/bulk-rnaseq
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
        - section: Sample Conditions for Differential Expression Analysis (Control vs Treatment, etc.)
          flow:
            - text: >-
                  You can (optionally) group samples into condition groups so
                  that downstream differential expression can be run on the
                  resulting sample transcript counts. For example, labeling a
                  subset of your samples as "Treatment" and another subset as
                  "Control" will yield a list of transcripts/genes that are
                  statistically different between the two groups.
            - fork: conditions_source
              flows:
                none:
                    display_name: No Differential Expression
                    flow:
                    - text: >-
                          Select "Manual Input" or "File" to construct your
                              condition groups.
                manual:
                    display_name: Manual Input
                    flow:
                    - params:
                        - manual_conditions
                    - params:
                        - run_splicing
                table:
                    display_name: File
                    _tmp_unwrap_optionals:
                        - conditions_table
                        - design_matrix_sample_id_column
                        - design_formula
                    flow:
                    - text: >-
                        Table with sample IDs and experimental conditions
                    - params:
                        - conditions_table
                        - design_matrix_sample_id_column
                        - design_formula
                    - params:
                        - run_splicing
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

        manual_conditions:

          __metadata__:
            display_name: Apply conditions to your samples:

        conditions_table:

          __metadata__:
            display_name: Design Matrix
            appearance:
                batch_table_column: true

        design_matrix_sample_id_column:

            __metadata__:
              display_name: Sample ID Column

        design_formula:

            __metadata__:
              display_name: Design Formula


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

        run_splicing:

          __metadata__:
            display_name: Run Differential Splicing Analysis

        custom_output_dir:
          You can provide a custom location where this run's analysis outputs
          will be located.

          __metadata__:
            display_name: Custom Output Location
    """
