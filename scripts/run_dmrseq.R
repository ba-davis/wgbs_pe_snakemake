#!/usr/bin/Rscript

source("scripts/dmrseq_functions.R")

#---------------#
# Load Libraries

library(bsseq)
library(tidyverse)
library(dmrseq)

#---------------#
# Get Command Line Arguments

args <- commandArgs(trailingOnly = TRUE)

# path to bsseq object created during bsseq ide
bs_rds <- args[1]

# path to design txt file
design_file <- args[2]

# desired outdir
outdir <- args[3]

# minimum CpG coverage when filtering
min_cov <- as.numeric(args[4])

# minimum percentage of samples CpG must have min_cov in
# using only the samples involved in the specific comparison
perc_samp <- as.numeric(args[5])

# pval cutoff
cutoff <- as.numeric(args[6])

# test covariate
# column in design file
test_cov <- args[7]

# match covariate
# column in design file
match_cov <- args[8]

# check for NULL value
if (match_cov == "NULL") {
  match_cov <- NULL
}

# number of cores to use for dmrseq
num_cores <- args[9]

#---------------#
# Set up Directory Structure

workdir <- getwd()

# set working directory to output folder
dir.create(file.path(outdir), showWarnings = FALSE)
setwd(file.path(outdir))
# create output directory for RDS objects
dir.create(file.path("RData"), showWarnings = FALSE)
# set workdir back to main folder
setwd(file.path(workdir))

#---------------#
# Read in Data

# read in BSseq object
bs <- readRDS(bs_rds)

# read in design file
design_df <- read.delim(design_file, header=T)
# ensure correct colnames of required columns
colnames(design_df)[1:5] <- c("sample_name",
  "group", "comparison", "comparison_designation",
  "comparison_name")

#---------------#
# Process Data

# Set workdir to outdir
setwd(file.path(outdir))

# Run each comparison in the design file
run_design_comps(bs, design_df,
  min_cov = min_cov,
  perc_samp = perc_samp,
  cutoff = cutoff,
  test_cov = test_cov,
  match_cov = match_cov,
  num_cores = num_cores)

sink("dmrseq_diff_complete.txt")
print("dmrseq diff analysis is complete.")
sink()
