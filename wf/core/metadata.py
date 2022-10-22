from latch.types import (Fork, ForkBranch, LatchAuthor, LatchMetadata,
                         LatchParameter, Params, Section, Spoiler, Text)

metadata = LatchMetadata(
    display_name="Bulk RNASeq",
    wiki_url="https://www.latch.wiki/bulk-rna-seq-end-to-end",
    video_tutorial="https://www.loom.com/share/dfba09ba6f524722b5d829f2424a3a3f",
    author=LatchAuthor(
        name="LatchBio",
        email="help@latch.bio",
        github="github.com/latchbio",
    ),
    repository="github.com/latch-verified/bulk-rnaseq",
    license="MIT",
    parameters={
        "samples": LatchParameter(
            display_name="Sample Sheet",
            description=(
                "Here you can organize your FastQ files by sample and add technical replicates for each "
                "sample. Biological replicates should be organized as separate samples."
            ),
            batch_table_column=True,
            _custom_ingestion="auto",
        ),
        "manual_conditions": LatchParameter(
            display_name="Apply conditions to your samples:"
        ),
        "conditions_table": LatchParameter(
            display_name="Design Matrix",
            batch_table_column=True,
        ),
        "design_matrix_sample_id_column": LatchParameter(
            display_name="Sample ID Column"
        ),
        "design_formula": LatchParameter(display_name="Design Formula"),
        "alignment_quantification_tools": LatchParameter(
            display_name="Alignment and Quantification Method"
        ),
        "latch_genome": LatchParameter(
            display_name="Genome Database Option",
            description="Curated reference files for specific genome sources and builds.",
            batch_table_column=True,
        ),
        "sa_ref_genome_fork": LatchParameter(
            display_name="Reference Genome Source",
            description="Select a reference genome from our curated database or provide your own.",
        ),
        "ta_ref_genome_fork": LatchParameter(
            display_name="Reference Genome Source",
            description="Select a reference genome from our curated database or provide your own.",
        ),
        "custom_ref_genome": LatchParameter(
            display_name="Reference Genome File",
            description="The reference genome you want to align you samples to.",
            detail="(.fasta, .fasta.gz, .fa, .fa.gz, .fna, .fna.gz)",
        ),
        "custom_gtf": LatchParameter(
            display_name="Annotation File",
            description="The gene annonation file that corresponds to the reference genome provided.",
            detail="(.gtf)",
        ),
        "custom_ref_trans": LatchParameter(
            display_name="Reference Transcript File (optional)",
            description="If not provided the workflow will generate from the Annotation File and Reference Genome File.",
            detail="(.fasta, .fasta.gz, .fa, .fa.gz, .fna, .fna.gz)",
        ),
        "star_index": LatchParameter(
            display_name="Provide Prebuilt STAR Index",
            description=(
                "You are able to provide a zipped prebuilt STAR alignment index for your genome. This will "
                "speed up run time as the index is generated if none is provided. In output settings you are able to save "
                "indices from a run to be used in future runs."
            ),
        ),
        "salmon_index": LatchParameter(
            display_name="Salmon Index",
            description=(
                "You are able to provide a zipped prebuilt Salmon selective alignment index for your genome. "
                "This will speed up run time as the index is generated if none is provided."
            ),
            detail="(.tar.gz is only accepted extension)",
        ),
        "save_indices": LatchParameter(
            display_name="Save Generated Reference Indexes",
            description=(
                "If you provided a custom genome you can output the alignment indexes generated from this run "
                "for use in future runs. This will speed up runtime since the workflow doesn't have to then regenerate the "
                "indexes."
            ),
        ),
        "run_name": LatchParameter(
            display_name="Run Name",
            description="A name for this analysis run, this will be used to name outputs from this run.",
            batch_table_column=True,
        ),
        "output_location_fork": LatchParameter(),
        "run_splicing": LatchParameter(
            display_name="Run Differential Splicing Analysis"
        ),
        "custom_output_dir": LatchParameter(
            display_name="Custom Output Location",
            description="You can provide a custom location where this run's analysis outputs will be located.",
        ),
    },
    flow=[
        Section(
            "Samples",
            Text(
                "Sample files can be provided and their read type can be inferred from their name or this information "
                "can be specified manually. Sample strandedness is inferred automatically (learn more)."
            ),
            Params("samples"),
        ),
        Section(
            "Sample Conditions for Differential Expression Analysis (Control vs Treatment, etc.)",
            Text(
                "You can (optionally) group samples into condition groups so that downstream differential expression "
                "can be run on the resulting sample transcript counts. For example, labeling a subset of your samples "
                'as "Treatment" and another subset as "Control" will yield a list of transcripts/genes that are '
                "statistically different between the two groups.",
            ),
            Fork(
                "conditions_source",
                "Conditions Source",
                none=ForkBranch(
                    "No Differential Expression",
                    Text(
                        'Select "Manual Input" or "File" to construct your condition groups.'
                    ),
                ),
                manual=ForkBranch(
                    "Manual Input",
                    Params("manual_conditions"),
                    Params("run_splicing"),
                ),
                table=ForkBranch(
                    "File",
                    Text("Table with sample IDs and experimental conditions"),
                    Params(
                        "conditions_table",
                        "design_matrix_sample_id_column",
                        "design_formula",
                    ),
                    Params("run_splicing"),
                ),
            ),
        ),
        Section(
            "Alignment and Quantification",
            Text(
                "This workflow uses Salmon's selective alignment described in this "
                "[paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8), which achieves "
                "greater accuracy than traditional alignment methods while using less computational resources."
            ),
            Fork(
                "alignment_quantification_tools",
                "Alignment Quantification Tools",
                traditional=ForkBranch(
                    "Selective Alignment",
                    Fork(
                        "ta_ref_genome_fork",
                        "Reference Genome Source",
                        database=ForkBranch(
                            "Select from Latch Genome Database",
                            Text(
                                "We have curated a set of reference genome data for ease and reproducibility. More "
                                "information about these managed files can be found "
                                "[here](https://github.com/latchbio/latch-genomes)."
                            ),
                            Params("latch_genome"),
                        ),
                        custom=ForkBranch(
                            "Provide Custom Genome",
                            Params(
                                "custom_ref_genome",
                                "custom_gtf",
                            ),
                            Spoiler(
                                "Optional Params",
                                Text(
                                    "These files will be generated from the GTF/Genome files if not provided."
                                ),
                                Params("salmon_index", "custom_ref_trans"),
                            ),
                        ),
                    ),
                ),
            ),
        ),
        Section(
            "Output Location",
            Params("run_name"),
            Fork(
                "output_location_fork",
                "Output Location",
                default=ForkBranch(
                    "Default",
                    Text(
                        'Output will be at default location in the data viewer - RNA-Seq Outputs/"Run Name"'
                    ),
                ),
                custom=ForkBranch("Specify Custom Path", Params("custom_output_dir")),
            ),
        ),
    ],
)
