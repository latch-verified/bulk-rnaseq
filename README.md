<html>
<p align="center">
  <img src="https://user-images.githubusercontent.com/31255434/182674961-b5b9ce91-ef56-48e7-80d1-cca029d25f78.jpg" alt="Latch Verified" width="100">
</p>

<h1 align="center">
  Bulk RNA-seq
</h1>

<p align="center">
<strong>
Latch Verified
</strong>
</p>

<p align="center">
  Produce transcript/count matrices from sequencing reads.
</p>

<p align="center">
  <a href="https://github.com/latch-verified/bulk-rnaseq/releases/latest">
    <img src="https://img.shields.io/github/release/latch-verified/bulk-rnaseq.svg" alt="Current Release" />
  </a>
  <a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/LICENSE-MIT-brightgreen.svg" alt="License" />
  </a>
  <img src="https://img.shields.io/github/commit-activity/w/latch-verified/bulk-rnaseq.svg?style=plastic" alt="Commit Activity" />
  <img src="https://img.shields.io/github/commits-since/latch-verified/bulk-rnaseq/latest.svg?style=plastic" alt="Commits since Last Release" />
</p>

<h3 align="center">
  <a href="https://console.latch.bio/explore/65992/info">Hosted Interface</a>
  <span> · </span>
  <a href="https://docs.latch.bio">SDK Documentation</a>
  <span> · </span>
  <a href="https://join.slack.com/t/latchbiosdk/shared_invite/zt-193ibmedi-WB6mBu2GJ2WejUHhxMOuwg">Slack Community</a>
</h3>

</html>


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

We have yet to identify a comprehensive review of the various trimming tools to
benchmark both accuracy and speed, so
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
