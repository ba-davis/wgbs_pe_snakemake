
#library(data.table)

# annot output table
# output of annotatePeaks.pl, want to keep only sig TF's

# knownResults.txt file (table holding TF's and their qvals)

# qval

# VARIABLES
# annot.out <- annot output table
# known.res <- knownResults.txt
# qval <- 0.05


subset_annot <- function(annot.out, known.res, qval=0.05, remove.homer.annot.cols=TRUE) {

  # read in the annot.out
  anno <- read.delim(annot.out, header=T)
  # read in the known results txt file
  res <- read.delim(known.res, header=T)

  # subset results to only TF's with qval < specified (default 0.05)
  res <- res[res$q.value..Benjamini. < qval, ]

  # Check if any TF's pass qval cutoff
  if (nrow(res) > 0) {
    # convert special characters of motif.name to "."
    res$Motif.Name <- gsub("\\(", ".", res$Motif.Name)
    res$Motif.Name <- gsub("\\)", ".", res$Motif.Name)
    res$Motif.Name <- gsub("\\/", ".", res$Motif.Name)
    res$Motif.Name <- gsub("-", ".", res$Motif.Name)
    res$Motif.Name <- gsub(":", ".", res$Motif.Name)
    res$Motif.Name <- gsub("\\?", ".", res$Motif.Name)
    res$Motif.Name <- gsub(" ", ".", res$Motif.Name)
    res$Motif.Name <- gsub(",", ".", res$Motif.Name)
    res$Motif.Name <- gsub("\\+", ".", res$Motif.Name)

    # remove the suffix from the TF column names of anno result table
    colnames(anno) <- gsub(".Distance.From.Peak.sequence.strand.conservation.", "", colnames(anno))

#----------------#
    # this example has 315 sig TF motifs
    # the html file (and corresponding knownResults/motif files) contains a max of 291 motifs, ordered by qval
    # the anno out table contains 21 columns + 291 motif columns = 312 columns
    # the res table contains 315 rows (motifs with qval < 0.05)
    # there are 24 motifs which appear in res$Motif.Name but NOT in anno columns
    
    # Now, remove columns from anno result table which are not found in the Motif.Name column of res
    # sanity check: if an error results from selecting columns from anno which appear in res$Motif.Name,
    #               then print the input file and columns causing the mismatch
    #if (length(res$Motif.Name[!(res$Motif.Name %in% colnames(anno))]) > 0) {
    #  print(paste0("Error in TF name matching for ", annot.out))
    #  print(res$Motif.Name[!(res$Motif.Name %in% colnames(anno))])
    #}
    
    # the above error seems to result if there are > 291 sig motifs
    # sanity check: if there are > 291 sig motifs, subset res to only include motifs appearing in anno
    # this example has 315 sig TF motifs (taken from subsetting knownResults.txt file)
    if (nrow(res) > 291) {
      # remove rows (motifs) from res which do not appear in anno
      res2 <- res[res$Motif.Name %in% colnames(anno), ]
      res <- res2
    }

#-------------#
    # Now, remove columns from anno result table which are not found in the Motif.Name column of res
    anno.sub <- anno[ ,res$Motif.Name]
    anno.sub2 <- cbind(anno[ ,c(1:21)], anno.sub)
    # convert blank values to NA
    anno.sub3 <- apply(anno.sub2, 2, function(x) gsub("^$|^ $", NA, x))
    # remove rows where there are NA values in every TF column
    num <- (ncol(anno.sub3) - 21) - 1
    anno.sub4 <- anno.sub3[!(rowSums(is.na(anno.sub3[ ,c(22:ncol(anno.sub3))])) > num), ]

    if (remove.homer.annot.cols) {
      anno.sub4 <- anno.sub4[ ,c(1,15,16,18,19,22:ncol(anno.sub4))]
    }
    colnames(anno.sub4)[1] <- "DMR_ID"
    
    # anno result table now contains:
    # rows genes from the motif search analysis which contain at least one TF motif (with qval < qval (0.05)) in the promoter
    # columns: a bunch of default output annot columns and only sig TF motif columns

    outfile <- gsub(".txt", "", annot.out)
    outfile <- paste0(outfile, ".sigTFs.txt")
    write.table(anno.sub4, outfile, sep="\t", col.names=T, row.names=F, quote=F)
  }
  else if (nrow(res) == 0) {
    print("There are zero TFs passing specified qval cutoff.")
    print(paste0("Output for ", annot.out, " remains unchanged, containing all TFs with no qval cutoff."))
  }
}
