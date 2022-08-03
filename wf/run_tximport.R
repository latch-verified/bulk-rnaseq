#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop("One argument must be supplied (GTF file).n", call. = FALSE)
}
library(tximport)
library(GenomicFeatures)

quant_path <- args[2]
gtf_path <- args[3]
out_path <- args[4]

txdb <- makeTxDbFromGFF(gtf_path, format = "gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

txi <- tximport(quant_path, type = "salmon", tx2gene = tx2gene)

counts <- txi$counts

df <- data.frame(Name = row.names(counts), NumReads = counts)

write.table(df, file = out_path, quote = FALSE, row.names = FALSE, sep = "\t")
