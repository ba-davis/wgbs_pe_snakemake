
#!/usr/bin/Rscript

source("scripts/chipseeker_functions.R")

#print(sessionInfo())

#---------------#
# Load Libraries

library(GenomicFeatures)
library(ChIPseeker)
library(writexl)
library(dplyr)

#---------------#
# Get Command Line Arguments

args <- commandArgs(trailingOnly = TRUE)

# path to DMR txt files to be annotated
inpath <- args[1]

# gtf file
gtf_file <- args[2]

# gene info file
gene_info_file <- args[3]
print(gene_info_file)

data_source <- args[4]
organism <- args[5]
organism <- gsub("_", " ", organism)

# check if gene info is NULL
if (gene_info_file == "NULL") {
  gene_info <- NULL
}
print(gene_info)

#---------------#
# Read in Data

# Read in the sig DMR txt files to be annotated
myfiles <- list.files(inpath, pattern = "sigDMRs.txt$", full.names = TRUE)
# set names of myfiles based on filenames
comps <- gsub(paste0(inpath, "/"), "", myfiles)
comps <- gsub(".sigDMRs.txt", "", comps)

# Read in all files and store as a list of data frames
df_list <- lapply(myfiles, read.table, header = TRUE)
names(df_list) <- comps
# change column2 name to seqnames
df_list <- lapply(df_list, function(df) {
  colnames(df)[2] <- "seqnames"
  return(df)
})

# also read in as GRanges objects for annotation
gr_list <- lapply(myfiles, read_and_convert_to_GRanges)
names(gr_list) <- comps

# make a TxDB object from the annotation gtf file
txdb <- makeTxDbFromGFF(gtf_file,
  format="gtf", data_source, organism)

#---------------#
# Process Data

# Annotate with Chipseeker and use custom function
# to re-annotate input regions that overlap a gene
#res_list <- lapply(gr_list, new_annotate)
#names(res_list) <- comps
res_list <- lapply(names(gr_list), function(gr_name) {
  new_annotate(gr_list[[gr_name]], gr_name, inpath)
})

# combine with bioinfo data to get gene info
#if (gene_info) {
  #combine with gene info
  #ref <- read.delim(gene_info_file, header=T)
  #colnames(ref)[1] <- "geneId"
  #foo2 <- merge(foo, ref, by="geneId")
  #  clean cols
  #foo3 <- foo2[ ,c(7,2,3,4,5,6,8,9,10,11,12,13,14,15,16,1,17,18,22,23,24)]
#} else {
#  print("no gene info provided")
#}

# check for duplicate enrtries for DMR
#foo4 <- subset(foo3, !duplicated(foo3[c("unq_ID", "annotation", "Gene.name")]))

# Check for un-annotated input regions
# Add un-annotated input regions to annotated results
for (i in 1:length(df_list)) {
  df <- df_list[[i]]
  res <- res_list[[i]]
  missing_rows <- anti_join(df, res, by = "unq_ID")
  # Add missing rows to res and fill missing columns with NA values
  if (nrow(missing_rows) > 0) {
    res <- bind_rows(res, missing_rows)
  }
  res_list[[i]] <- res
}

# Order based on unq_ID
# And set unq_ID as first column
for (i in 1:length(res_list)) {
  res <- res_list[[i]]
  # Extract the numeric part of the "unq_ID" column
  numeric_part <- as.numeric(gsub("DMR_", "", res$unq_ID))
  # Order the rows based on the numeric part
  ordered_indices <- order(numeric_part)
  # Reorder the data frame
  res <- res[ordered_indices, ]
  # set unq_ID as first column
  res <- res[, c("unq_ID", setdiff(names(res), "unq_ID"))]
  res_list[[i]] <- res
}

# export annotated results txt files
for (i in 1:length(res_list)) {
  res <- res_list[[i]]
  write.table(res,
  paste0(inpath, "/", names(res_list)[i], ".sigDMRs.annot.txt"),
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)
}

# export excel file
write_xlsx(res_list,
  paste0(inpath, "/annotated_res.xlsx"))

sink("data/diff/methylkit_dmr/annotation_complete.txt")
print("Annotation is complete.")
sink()
