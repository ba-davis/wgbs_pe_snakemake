#!/usr/bin/Rscript


#---------------#
# Load Libraries

library(bsseq)
library(tidyverse)
library(dmrseq)

#---------------#
# Get Command Line Arguments

args <- commandArgs(trailingOnly = TRUE)

# path to RDS objects from dmrseq comparisons
indir <- args[1]

# desired output directory
outdir <- args[2]

# export qval cutoff
fdr <- as.numeric(args[3])

# test covariate
test_cov <- args[4]

# genome for plotting dmrs
dmr_genome <- args[5]

# max number of sig DMRs to plot
num_plot <- as.numeric(args[6])

#---------------#
# Set up Directory Structure

workdir <- getwd()

#---------------#
# Process Data

#setwd(file.path(paste0(outdir, "/RData")))
print(getwd())

# filtered bsseq object files
bs_files <- list.files(indir, pattern = "filtered_bsseq.RDS$")
print(bs_files)
# regions RDS files
regions_files <- list.files(indir, pattern = "regions.RDS$")
print(regions_files)

# Check that matching pairs of files exist in proper order
identical(gsub(".filtered_bsseq.RDS", "", bs_files),
  gsub(".regions.RDS", "", regions_files))

# Obtain comparison names from RDS file names
comps <- gsub(".filtered_bsseq.RDS", "", bs_files)
print(comps)

# Read in data
setwd(file.path(paste0(outdir, "/RData")))
print(getwd())
bs_list <- lapply(bs_files, readRDS)
regions_list <- lapply(regions_files, readRDS)
print(length(bs_list))
print(length(regions_list))
# Check that the length of both lists is the same
identical(length(bs_list), length(regions_list))

# Create anno track for plotting
print(paste0("Attempting to getAnnot annoTrack for ", dmr_genome))
annoTrack <- getAnnot(dmr_genome)

#-----#

#print(outdir)
#print(getwd())
# Set workdir to outdir
setwd("..")
print(getwd())

# Distribution of methylation changes in all DMRs (not sig DMRs)
# Note that assymmetric changes in the distribution indicate global shifts in methylation between the groups
minMethDiff <- 0.1
print(length(regions_list))
head(as.data.frame(regions_list[[1]]))
head(as.data.frame(regions_list[[2]]))
head(as.data.frame(regions_list[[3]]))
for (i in 1:length(regions_list)) {
  g <- ggplot(as.data.frame(regions_list[[i]], stringsAsFactors=FALSE), aes(x=beta)) +
    geom_histogram() +
    labs(x="Methylation Difference (per-DMR)") +
    geom_vline(xintercept=minMethDiff) + geom_vline(xintercept=-1*minMethDiff) +
    theme_minimal()
  ggsave(paste0(comps[i], ".dist_beta.pdf"), g)
}

#-----#

# Similarly, the test statistic, which is used to compute the p-value, is shown below
# Its interpretation is the same as that above
for (i in 1:length(regions_list)) {
  g <- ggplot(as.data.frame(regions_list[[i]], stringsAsFactors=FALSE), aes(x=stat)) +
    geom_histogram() +
    labs(x="Methylation Difference Statistic") +
    theme_minimal()
  ggsave(paste0(comps[i], ".dist_stat.pdf"), g)
}

#-----#

# Pvalue distribution
for (i in 1:length(regions_list)) {
  g <- ggplot(as.data.frame(regions_list[[i]], stringsAsFactors=FALSE), aes(x=pval)) +
    geom_histogram() +
    labs(x="Unadjusted p-value") +
    theme_minimal()
  ggsave(paste0(comps[i], ".dist_pval.pdf"), g)
}

#-----#

# Subset to sig DMRs for export

for (i in 1:length(regions_list)) {
  print(paste0("Attempting to subset sig DMRs from: ", comps[i]))
  
  export_dmrs <- regions_list[[i]][regions_list[[i]]$qval < fdr,]
  num_sig <- length(export_dmrs)
  print(paste0("Number of regions below fdr ", fdr, ": ", num_sig))

  if (num_sig > 0) {
    # get percent of hyper to hypo methylation of sig DMRs
    print(paste0("Percent of sig DMRs that are hypermethylated: ",
      round(100 * sum(export_dmrs$stat > 0) / length(export_dmrs), 1), "%"))
    print(paste0("Percent of sig DMRs that are hypomethylated: ",
      round(100 * sum(export_dmrs$stat < 0) / length(export_dmrs), 1), "%"))

    # get raw mean methylation differences
    rawDiff <- meanDiff(bs_list[[i]], dmrs = export_dmrs, testCovariate = test_cov)

    sig_dmrs <- as.data.frame(export_dmrs)
    sig_dmrs$mean_meth_diff <- rawDiff

    write.table(sig_dmrs, paste0(comps[i], ".sig_dmrs.txt"),
      sep = "\t", col.names=T, row.names=F, quote=F)

    if (num_sig > num_plot) {
      for (j in 1:num_plot) {
        pdf(paste0(comps[i], "_DMR_", j, "_plot.pdf"))
        plotDMRs(bs_list[[i]], regions=export_dmrs[j,], testCovariate=test_cov,
        annoTrack=annoTrack)
        dev.off()
      }
    } else if (num_sig <= num_plot) {
      for (j in 1:num_sig) {
        pdf(paste0(comps[i], "_DMR_", j, "_plot.pdf"))
        plotDMRs(bs_list[[i]], regions=export_dmrs[j,], testCovariate=test_cov,
        annoTrack=annoTrack)
        dev.off()
      }
    }
  } else if (num_sig < 1) {
    print(paste0("No significant DMRs according to FDR cutoff of ", fdr))
    print("Consider relaxing fdr cutoff.")
  }
}

sink("dmrseq_explore_complete.txt")
print("exploration and export of dmrseq dmrs is complete.")
sink()
