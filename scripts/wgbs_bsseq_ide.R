
# conda activate bsseq

# Exploratory Data Analysis of WGBS Data

# Reads in bismark cytosine reports (deduped for wgbs) via
#   bsseq::read.bismark()
#   By default, does not remove zero cov cpgs, and collapses strands

# Inputs:
#   - directory pointing to .CpG_report.txt.gz files from bismakr meth extract
#   - metadata excel file
#       - first column must be "sample_name" and match sample names prefixing the input files
#   - cores: number of cores to use when reading in input files and creating bsseq obj

# Load Libraries
library(bsseq)
library(tidyverse)
library(readxl)

# source my custom functions for pca
source("general_functions.R")
source("bsseq_ide_functions.R")

# Input Variables
cores <- 8
indir <- "data/meth_extract"
files <- list.files(indir,
    pattern = "*.CpG_report.txt.gz",
    full = TRUE)
coldata <- as.data.frame(read_xlsx("data/metadata/metadata.xlsx"))
pca_cov <- 1
pca_color_var <- "Group"
pca_shape_var <- "NULL"
pca_label_var <- "sample_name"

# read in cytosine reports (deduped)
# collapse strands (this halves the number of loci)
bs <- bsseq::read.bismark(files = files,
    rmZeroCov = FALSE,
    strandCollapse = TRUE,
    verbose = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(workers = cores,
        progressbar = FALSE),
    nThread = 1)

# clean sample names obtained from filenames
sampleNames(bs) <- gsub(paste0(indir, "/"), "", sampleNames(bs))
sampleNames(bs) <- gsub("_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.gz",
    "", sampleNames(bs))

# add metadata to bsseq obj
# check that the sample names obtained from the cytosine report file names
# match the sample names from the first column of the metadata file
if (!identical(sampleNames(bs), coldata$sample_name)) {
  print("CAUTION. Sample Names in metadata file do not match file names.")
  print("Check the sample_name column and change to match input file name prefixes.")
} else {
  pData(bs) <- cbind(pData(bs), coldata)
}

#---------------------------------------------------#
# Useful IDE things
# Coverage Plots
# Meth Stats Plots

# Filter by cov (low cov and extremely high cov)

# PCA
# (on smoothed methylation values?)
# no missing values
# if you have too many PC's, maybe try NMF or tSNE

# seems like general bssmooth pipeline is:
#   - I guess keep all cpgs, even if they have 0 coverage
#   - smoothing is done separately on each sample:
#     - only uses data where the coverage for that sample is nonzero
#     - this estimates a genome wide methylation profile, which is
#       then evaluated for all cpgs in the bsseq object
#     - thus, after smoothing, every cpg in the bsseq object as a
#       methylation value
#   - note that smoothed methylation values make less sense if they are
#     in an area with no covered cpgs nearby. Remove cpgs after smoothing.
#   - example: keep cpgs where majority of samples per group have
#     at least 2x cov (DMRichr has 1x as default)

#-----------------------#

# Coverage plots

# Plot a cov distribution barplot for every sample in a bsseq obj
plot_cov_dist_bs(bs)

#-----------------------#

# Meth Stats Plots

# Plot a methylation percentage distribution barplot for every
# sample in a bsseq obj
plot_meth_stats_bs(bs)

#-----------------------#

# PCA

# subset to cpgs with at least 1X cov in all samples

# Get read coverage for each CpG
cov <- getCoverage(bs)
# Compute the number of samples
n_samples <- dim(cov)[2]
# Subset the BSseq object to only include CpGs with at least 1 read coverage in all samples
keep <- apply(cov, 1, function(x) sum(x >= pca_cov) == n_samples)
bs_sub <- bs[keep, ]

num_cpg <- length(bs_sub)
print(paste0("For PCA, using ", num_cpg, " CpGs with at least ", pca_cov, "X coverage in all samples."))

# get percent meth df of filtered bsseq obj
meth_df <- as.data.frame(bsseq::getMeth(bs_sub, type = "raw"))

# Perform PCA on the percent methylation matrix
# prcomp expects samples to be rows and features to be columns, so transpose
# remove cpgs (columns) with zero variance
mypca <- meth_df %>%
    t() %>%
    .[, which(apply(., 2, var) != 0)] %>%
    prcomp(., center = TRUE, scale = TRUE)

describe_pca(mypca)

# grab the principle components as a new df
my_pca_df <- as.data.frame(mypca$x)

# add the relevant variable columns from coldata to the pca df
my_pca_df2 <- cbind(my_pca_df, coldata)

# get eigen values (standard deviation squared)
eigs <- mypca$sdev^2

plot_pca(my_pca_df2, "PC1", "PC2",
    color_var=pca_color_var, shape_var=pca_shape_var, label_var=pca_label_var,
    eigs=eigs, num_cpg = num_cpg)

plot_pca(my_pca_df2, "PC1", "PC3",
    color_var=pca_color_var, shape_var=pca_shape_var, label_var=pca_label_var,
    eigs=eigs, num_cpg = num_cpg)

plot_pca(my_pca_df2, "PC2", "PC3",
    color_var=pca_color_var, shape_var=pca_shape_var, label_var=pca_label_var,
    eigs=eigs, num_cpg = num_cpg)
