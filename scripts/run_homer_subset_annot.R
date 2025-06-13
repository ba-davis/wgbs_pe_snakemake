#!/usr/bin/Rscript

source("scripts/subset_annot_table.R")

#---------------#
# Get Command Line Argumentss

args <- commandArgs(trailingOnly = TRUE)

# path to input directory files (bismark cov files)
homer_out_path <- args[1]

#---------------#

# Assumes the homer output dirs are the same as the annotPeaks output files with
# the suffix ".knownMotifs.all.annot.output.txt" removed"

clean_annotatePeaks_outputs <- function(homer_out_path) {
  # gather file names
  my.annot.files <- list.files(homer_out_path, "*knownMotifs.all.annot.output.txt$")

  # vector of the output directory names, in same order as my.annot.files
  homer.dirs <- gsub(".knownMotifs.all.annot.output.txt", "", my.annot.files)

  # Run the subset function on each file
  for (i in 1:length(my.annot.files)) {
    subset_annot(annot.out=paste0(homer_out_path, my.annot.files[i]),
                 known.res=paste0(homer_out_path, homer.dirs[i], "/", "knownResults.txt"),
                 qval=0.05
    )
  }
}

# Execute the function
clean_annotatePeaks_outputs(homer_out_path)

sink(paste0(homer_out_path, "homer_subset_annot_complete.txt"))
print("HOMER subset annot is complete.")
sink()
