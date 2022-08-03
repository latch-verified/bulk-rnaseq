#!/usr/bin/env Rscript
if (!require("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite")
}
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("tximport")
BiocManager::install("GenomicFeatures")

library(tximport)
library(GenomicFeatures)
