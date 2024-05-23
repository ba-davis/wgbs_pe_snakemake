#!/usr/bin/Rscript

# Execute MethylKit DMR Analysis

source("scripts/dmr_functions.R")

#---------------#
# Get Command Line Arguments

args <- commandArgs(trailingOnly = TRUE)

# methylRawList object RDS file
my_obj_rds <- args[1]

lo_count <- args[2]
hi_percent <- args[3]

cov_bases <- args[4]

tile_mpg <- args[5]

design_file <- args[6]

outdir <- args[7]

# check for NULL value
if (hi_percent == "NULL") {
  hi_percent <- NULL
} else {
  # Convert the value to numeric
  hi_percent <- as.numeric(hi_percent)
}
# check tile_mpg param
if (tile_mpg == "NULL") {
  tile_mpg <- NULL
} else {
  # Convert the value to integer
  tile_mpg <- as.integer(tile_mpg)
}

#---------------#
# Set up Directory Structure

# set working directory to output folder
dir.create(file.path(outdir), showWarnings = FALSE)
setwd(file.path(outdir))
# create output directory for RDS objects
dir.create(file.path("RData"), showWarnings = FALSE)
# set workdir back to main folder
setwd("../../..")

#---------------#
# Read in Data

print("Reading in RDS file.")

# read in methylrawlist file
my_obj <- readRDS(my_obj_rds)

# read in design file
design_df <- read.delim(design_file, header=T)
# ensure correct colnames of required columns
colnames(design_df)[1:5] <- c("sample_name",
  "group", "comparison", "comparison_designation",
  "comparison_name")

print("Coverage Filtering")

# Coverage Filter CpGs
myObj.filtered <- filterByCoverage(my_obj,
  lo.count = lo_count,
  lo.perc = NULL,
  hi.count = NULL,
  hi.perc = hi_percent)

print("Normalizing median.")

# Normalize by Coverage
myObj.filtered.normalized <- normalizeCoverage(myObj.filtered, method = "median")

#---------------#
# Tile

print("Creating tiles per sample now.")

my.tiles <- tileMethylCounts(myObj.filtered.normalized,
  win.size = 1000,
  step.size = 1000,
  cov.bases = cov_bases,
  mc.cores = 8)

#---------------#
# Merge

print("Merging tiles now.")

meth <- methylKit::unite(my.tiles,
  min.per.group = tile_mpg,
  mc.cores = 8)

print(paste0("Merging complete. Number of tiles: ", nrow(meth)))

# export merged object
saveRDS(meth, paste0(outdir, "/RData/dmr_meth.RDS"))

#---------------#
# Reorganize and Calculate DMRs

# For each comparison, need:
# 1. control group sample names
# 2. treatment group sample names
# 3. treatment vector (define based off length of above 2 variables)
# 4. covariate data frame, same order as above sample names
# 5. name of comparison for exported files
# 6. meth.diff cutoff (10)
# 7. qvalue cutoff (0.1 or 0.05)
# 8. DMR or DMC

# Set workdir to outdir
setwd(file.path(outdir))

# Run DMR function on each comparison in the design matrix
# Covariates will be included if they exist as additional columns in the design file

run_design_comps(meth = meth,
  design_df = design_df,
  meth_diff = 10,
  qval_cutoff = 0.1,
  dm_type = "DMR")

sink("diff_complete.txt")
print("Diff Analysis is complete.")
sink()

# create seq_run as covariate
#cov <- data.frame(seq_run=c("Run1","Run2","Run1","Run1","Run2","Run2","Run1","Run1","Run2","Run1","Run2","Run2"))
# cov not used below

# GR3 vs GR1
#c1 <- reorganize(meth,
#  sample.ids=c("2_A1","2_A2","4_A1","4_A2","6_A1","6_A2","8_A1","8_A2","10_A0","10_A1","12_A0",
#               "2_E1","2_E2","4_E1","4_E2","6_E1","6_E2","8_E0","8_E1","8_E2","10_E0","10_E1","12_E0"),
#  treatment=c(0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1)
#)

# Calculate DMRs
# exports:
#   1. Sig DMR results txt file
#   2. Sig DMR results bed file (ensembl contigs)
#   3. Hyper vs Hypo DMR pie chart
#   4. Both GREAT bed file
#   5. Hyper GREAT bed file
#   6. Hypo GREAT bed file
#res1 <- calc.DMRs(c1,
#                  covariate=NULL,
#                  overdispersion="MN",
#                  test="Chisq",
#                  comparison="GR3_v_GR1",
#                  meth.diff=10,
#                  qval=0.1,
#                  type="DMR",
#                  mc=1
#)
# [1] "nrow of unfiltered diff results is: 95406"
# [1] "Number of sig DMRs: 339"
# [1] "nrow of myDiff.sig.df is: 339"
# [1] "number of hyper DMRs: 173"
# [1] "number of hypo DMRs: 166"
#
# Make bed file track
#makeBED(res1, "GR3_v_GR1", "DMR")

