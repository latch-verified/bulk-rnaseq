#!/usr/bin/env Rscript
if (!require("jsonlite", quietly = TRUE)) {
  pak::pak("jsonlite")
}
if (!require("BiocManager", quietly = TRUE)) {
  pak::pak("BiocManager")
}

pak::pak(c(
  "tximport",
  "RCurl",
  "GenomicFeatures"
))

library("tximport")
library("GenomicFeatures")
