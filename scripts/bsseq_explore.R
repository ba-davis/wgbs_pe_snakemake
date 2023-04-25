
# conda activate bsseq

library(bsseq)
library(tidyverse)

# source my custom functions for pca
source("/home/groups/hoolock2/u0/bd/github_repos/wgbs_pe_snakemake/scripts/general_functions.R")
#source("general_functions.R")

cores <- 8
indir <- "/home/groups/hoolock2/u0/bd/Projects/ECP29/ECP80/wgbs_pe_snakemake/data/meth_extract"
files <- list.files(indir,
    pattern = "*.CpG_report.txt.gz",
    full = TRUE)

# read in cytosine reports (deduped)
# collapse strands (this halves the number of loci)
bs <- bsseq::read.bismark(files = files,
    rmZeroCov = FALSE,
    strandCollapse = TRUE,
    verbose = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(workers = cores,
        progressbar = FALSE),
    nThread = 1)

# clean sample names
sampleNames(bs) <- gsub(paste0(indir, "/"), "", sampleNames(bs))
sampleNames(bs) <- gsub("_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.gz",
    "", sampleNames(bs))

# add metadata to bsseq obj
meta <- data.frame(sample_name = sampleNames(bs),
    Group = c(rep("E12.5", 3), "E15.5-distal",
    "E15.5-proximal", "E15.5-distal", "E15.5-proximal",
    "E15.5-distal", "E15.5-proximal", "E17.5-distal",
    "E17.5-proximal", "E17.5-distal", "E17.5-proximal",
    "E17.5-distal", "E17.5-proximal"),
    location = c("whole", "whole", "whole", "distal",
    "proximal", "distal", "proximal",
    "distal", "proximal", "distal",
    "proximal", "distal", "proximal",
    "distal", "proximal"),
    replicate = c("1", "2", "3", "1", "1", "2", "2",
    "3", "3", "1", "1", "2", "2", "3", "3"))
pData(bs) <- cbind(pData(bs), meta[2:length(meta)])

#--------------------------------------------------#
# Useful functions for manipulating BSseq object

# look at first few genomic locations
head(granges(bs))
# look at methylation value matrix (number of reads showing methylation)
head(getCoverage(bs, type = "M"))
# look at the coverage matrix (total number of reads per loci)
head(getCoverage(bs))
# number of cpgs (regardless of coverage)
length(bs)
# number of cpgs which are covered by at least 1 read in all samples
sum(rowSums(getCoverage(bs) >= 1) == 15)
# number of CpGs with 0 coverage in all samples
sum(rowSums(getCoverage(bs)) == 0)
# the avg coverage of cpgs per sample
round(colMeans(getCoverage(bs)), 1)

#---------------------------------------------------#
# Useful IDE things:

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
plot_cov_dist_bs <- function(bs) {
    # diable scientific notation on the plot
    options(scipen = 999)

    # get coverage df
    loci_cov <- as.data.frame(bsseq::getCoverage(bs, type = "Cov"))

    for (i in seq_len(ncol(loci_cov))) {
        p <- ggplot(as.data.frame(table(loci_cov[, i])),
            aes(x = as.numeric(Var1), y = Freq)) +
            geom_bar(stat = "identity") +
            scale_x_log10() +
            scale_y_log10() +
            theme_minimal() +
            labs(title = colnames(loci_cov)[i],
            x = "Coverage", y = "Number of CpGs") +
            theme(plot.title = element_text(hjust = 0.5))
        ggsave(paste0("cov_plot_", colnames(loci_cov)[i], ".pdf"), p)
    }
}

plot_cov_dist_bs(bs)

#-----------------------#

# Meth Stats Plots

# Plot a methylation percentage distribution barplot for every
# sample in a bsseq obj
plot_meth_stats_bs <- function(bs) {
    # diable scientific notation on the plot
    options(scipen = 999)

    # get meth percent df
    meth_df <- as.data.frame(bsseq::getMeth(bs, type = "raw"))
    # round to 1 decimal place
    meth_df[] <- lapply(meth_df,
        function(x) round(x, 1))

    for (i in seq_len(ncol(meth_df))) {
        p <- ggplot(as.data.frame(table(meth_df[, i])),
            aes(x = Var1, y = Freq)) +
            geom_bar(stat = "identity") +
            theme_minimal() +
            labs(title = colnames(meth_df)[i],
            x = "Methylation Rate", y = "Number of CpGs") +
            theme(plot.title = element_text(hjust = 0.5))
        ggsave(paste0("meth_stats_plot_", colnames(meth_df)[i], ".pdf"), p)
    }
}

plot_meth_stats_bs(bs)

#-----------------------#

# PCA

# subset to cpgs with at least 1X cov in all samples

# Get read coverage for each CpG
cov <- getCoverage(bs)
# Compute the number of samples
n_samples <- dim(cov)[2]
# Subset the BSseq object to only include CpGs with at least 1 read coverage in all samples
keep <- apply(cov, 1, function(x) sum(x >= 1) == n_samples)
bs_sub <- bs[keep, ]

num_cpg <- length(bs_sub)

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
my_pca_df2 <- cbind(my_pca_df, meta)

# get eigen values (standard deviation squared)
eigs <- mypca$sdev^2

plot_pca(my_pca_df2, "PC1", "PC2",
    color_var="Group", shape_var="NULL", label_var="replicate",
    eigs=eigs, num_cpg = num_cpg)

plot_pca(my_pca_df2, "PC1", "PC3",
    color_var="Group", shape_var="NULL", label_var="replicate",
    eigs=eigs, num_cpg = num_cpg)

plot_pca(my_pca_df2, "PC2", "PC3",
    color_var="Group", shape_var="NULL", label_var="replicate",
    eigs=eigs, num_cpg = num_cpg)
    