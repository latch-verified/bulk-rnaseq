#!/usr/bin/env Rscript

install.packages('remotes', repos='http://cran.us.r-project.org')

install.packages("rstan")
remotes::install_github("stan-dev/rstantools")
BiocManager::install("DirichletMultinomial")
remotes::install_github("davidaknowles/leafcutter/leafcutter")

library("leafcutter")
