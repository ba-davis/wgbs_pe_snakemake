#!/usr/bin/Rscript

source("scripts/methylkit_functions.R")

#print(sessionInfo())

#---------------#
# Load Libraries

#library(methylKit)
#library(ggplot2)
#library(scales)
#library(ggrepel)
#library(dplyr)

#---------------#
# Get Command Line Arguments

args <- commandArgs(trailingOnly = TRUE)

# path to input directory files (bismark cov files)
inpath <- args[1]
#print(inpath)

# path to metadata.txt file
metadata <- args[2]
#print(metadata)

# name of metadata file column to use as group for methylkit
# For this, I need to pass a column from metadata that represents group
# membership, it's helpful if it has meaningful names because it can
# show up on plots
group_var <- args[3]
#print(group_var)

# minimum cov threshold for object read-in
min_cov <- as.numeric(args[4])
#print(min_cov)

# Set output folder
outdir <- args[5]
#print(outdir)

# minimum cov threshold for merging samples
min_cov_merge <- as.numeric(args[6])
# max cov percentage for merging samples
perc_merge <- args[7]
#print(perc_merge)
# logical on whether to merge on CpG-level or tile-level
merge_regional <- as.logical(args[8])
# if merging regions, minimum number of CpGs required to be covered in the region
min_cpg_region <- as.numeric(args[9])

# min per group for merging
mpg <- args[10]

# merge_dirname (name of directory to be created within data/ide to hold output on merged meth object)
merge_dirname <- args[11]

# whether to repeat initial IDE or not
repeat_initial_ide <- as.logical(args[12])

# Check if the values are "NULL" and change to NULL or numeric
if (perc_merge == "NULL") {
  perc_merge <- NULL
} else {
  # Convert the value to numeric
  perc_merge <- as.numeric(perc_merge)
}
# check mpg param
if (mpg == "NULL") {
  mpg <- NULL
} else {
  # Convert the value to integer
  mpg <- as.integer(mpg)
}

print(mpg)
#---------------#

# If it's the first time running IDE
# Set up directory structure
# Produce plots and descriptions before merging

if (repeat_initial_ide) {
  print("First time running IDE.")
  print("Will read in files to create and save methylRawList object.")
  print("Will explore coverage distribution of samples in the object.")
  
  # Set up ide directories for output

  # set working directory to output folder
  dir.create(file.path(outdir), showWarnings = FALSE)
  setwd(file.path(outdir))
  # create output directory for RDS objects
  dir.create(file.path("RData"), showWarnings = FALSE)

  #---------------#
  # Read in data

  #print(paste0(getwd()))
  setwd("../../") # set work dir to main wgbs_pe_snakemake dir
  #print(getwd())

  # redefine metadata filepath since we set work dir to data/ide
  #newmeta <- paste0("../../", metadata)

  meta <- read.delim(metadata, header = TRUE)
  #head(meta)

  # Set a different temporary directory
  #options(tmpdir = "/home/groups/hoolock2/u0/bd/Projects/mariam_r21_mouse_ad_2022/wgbs/hippocampus/wgbs_pe_snakemake/tmpdir")
  #print(tempdir())

  my_obj <- create_obj(inpath = inpath,
    metadata = metadata,
    group_var = group_var,
    min_cov = min_cov
  )
  #print("finished create_obj")
  #print(paste0("nrow my_obj: ", nrow(my_obj)))

  # export methylRawList object as RDS
  saveRDS(my_obj, paste0(outdir, "/RData/my_obj.RDS"))

  # read in my_obj
  #my_obj <- readRDS("data/ide/RData/my_obj.RDS")

  #---------------#
  # Explore Coverage Distributions

  setwd(file.path(outdir)) # set work dir to "data/ide"

  # create cov data object for plotting
  cov_dist_df <- get_cov_dist_df(my_obj)
  #print(head(cov_dist_df))
  # Plot a boxplot of all samples by number of cpgs covered
  plot_cov_dist_boxplot(cov_dist_df = cov_dist_df,
    filename = "num_cpg_boxplot.pdf",
    fill_col = "lightskyblue1"
  )
  # Plot stacked barplot showing CpG cov per sample
  plot_cov_bin_barplot(cov_dist_df = cov_dist_df,
    filename = "cpg_cov_stacked_barplot.pdf"
  )
  # Plot stacked barplot, faceted by group
  plot_cov_bin_barplot(cov_dist_df = cov_dist_df,
    filename = "cpg_cov_stacked_barplot_facet.pdf",
    group_facet = TRUE,
    metadata = "../../metadata.txt",
    groupbyvar = group_var
  )

  # Plot Methylkit style cov histograms, one plot per sample
  makeCovPlots(my_obj)

  # Plot Methylkit style meth stats histogram, one plot per sample
  makeMethStats(my_obj)

  # filter low coverage cpgs and produce cpg cov table at various cov values
  # the final filtered object will use the cutoffs set by cov5 and perc5
  myObj.filtered <- cpg_cov(my_obj,
    cov1 = 1, perc1 = NULL,
    cov2 = 5, perc2 = NULL,
    cov3 = 10, perc3 = NULL,
    cov4 = 40, perc4 = NULL,
    cov5 = 10, perc5 = 99.9,
    table = TRUE
  )

  # Plot number of CpGs >10X cov barplot
  number_cpg_at_cov_barplot(cpg_cov_table = "CpG.coverage.table.txt")
} else {

  print("Reading in existing methylRawList object.")
  print("Will proceed with merging samples in the object.")
  
  my_obj <- readRDS("data/ide/RData/my_obj.RDS")

  # read in metadata
  meta <- read.delim(metadata, header = TRUE)

  setwd(file.path(outdir))
}

#---------------#
# Merge Data

# create output folder to hold output conducted on merged data
dir.create(file.path(merge_dirname), showWarnings = FALSE)
setwd(file.path(merge_dirname))

print("attempting get_merged_regions now.")
# merge the data and keep cpgs found in all samples per group
meth <- get_merged_regions(myObj = my_obj,
  min = min_cov_merge,
  max = perc_merge,
  normalize = TRUE,
  regional = merge_regional,
  min_cpg = min_cpg_region,
  mpg = mpg,
  destrand = FALSE
)
# export methylRawBase object as RDS
saveRDS(meth, "meth.RDS")

# get the percent methylation matrix
perc_meth = percMethylation(meth)

# export methylRawBase object as RDS
saveRDS(perc_meth, "perc_meth.RDS")

# If mpg is not null, perform imputation on CpGs without coverage in all samples
# PCA requires no missing data
if (is.null(mpg)) {
  print("Min.per.group parameter is NULL. No imputation required.")
} else {
  print("Min.per.group parameter is not NULL. Performing imputation.")
  imputed.mat <- impute.knn(perc_meth)
  # Get percent methylation matrix (CpGs x Samples)
  imp.perc.meth <- imputed.mat$data
  perc_meth <- imp.perc.meth
}

#---------------#
# PCA

# Remove rows with zero variance
perc_meth_filtered <- perc_meth[apply(perc_meth, 1, var) != 0, , drop = FALSE]
sink("ide_pca_log.txt")
print(paste0("Number of zero variance rows removed for PCA: ",
  nrow(perc_meth) - nrow(perc_meth_filtered)))
# obtain new number of cpgs
num_cpg <- nrow(perc_meth_filtered)
print(paste0("Number of CpGs used for PCA: ", num_cpg))
sink()

# perform PCA using prcomp on the transposed perc_meth matrix
mypca <- prcomp(t(perc_meth_filtered),
  center=TRUE,
  scale=TRUE
)

# describe the pca results and capture number of PC's required to explain 80% variance
num_pc <- describe_pca(mypca)
# get eigen values (standard deviation squared)
eigs <- mypca$sdev^2

# grab the principle components as a new df
my_pca_df <- as.data.frame(mypca$x)
# set sample column
my_pca_df$sample_name <- rownames(my_pca_df)
# merge metadata into pca df
my_pca_df2 <- merge(my_pca_df, meta, by='sample_name')

# export the pca df as RDS
saveRDS(my_pca_df2, "pca_df.RDS")

# plot a few pca plots with custom function
plot_pca(my_pca_df2,
  pc_a = "PC1",
  pc_b = "PC2",
  color_var = group_var,
  shape_var = "NULL",
  label_var = "sample_name",
  eigs = eigs,
  num_cpg = num_cpg
)

plot_pca(my_pca_df2,
  pc_a = "PC1",
  pc_b = "PC3",
  color_var = group_var,
  shape_var = "NULL",
  label_var = "sample_name",
  eigs = eigs,
  num_cpg = num_cpg
)

plot_pca(my_pca_df2,
  pc_a = "PC2",
  pc_b = "PC3",
  color_var = group_var,
  shape_var = "NULL",
  label_var = "sample_name",
  eigs = eigs,
  num_cpg = num_cpg
)

setwd("..")

sink("ide_complete.txt")
print("IDE is complete.")
sink()
