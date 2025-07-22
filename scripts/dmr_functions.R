
# design_df is a data frame with the following required 5 columns:
#   Column order and column names matter
#   <sample_name> <group> <comparison> <comparison_designation> <comparison_name>
# Any additional columns will be considered as covariates
# An example design_df is below:

#   design_df <- data.frame(sample_name = c("gr1_r1", "gr1_r2", "gr1_r3", "gr1_r4",
#     "gr2_r1", "gr2_r2", "gr2_r3", "gr2_r4",
#     "gr1_r1", "gr1_r2", "gr1_r3", "gr1_r4",
#     "gr3_r1", "gr3_r2", "gr3_r3", "gr3_r4"),
#     group = c(rep("group1", 4), rep("group2", 4), rep("group1", 4), rep("group3", 4)),
#     comparison = c(rep("comp1", 8), rep("comp2", 8)),
#     comparison_designation = c(rep("control", 4), rep("treatment", 4),
#       rep("control", 4), rep("treatment", 4)),
#     comparison_name = c(rep("GR2_v_GR1", 8), rep("GR3_v_GR1", 8)),
#     sex = c(rep(c("M","M","F","F"), 4))
#   )

# Logic of code and DMR workflow:
#   - design_df is required as stated above.
#   - covariates will be used if they are provided as additional columns in design_df
#   - comparison names are taken from design_df
#   - overdispersion="MN" and test="Chisq" are hardcoded for now
#   - meth.diff and qval cutoff have defaults but can be changed by user

#---------------#
# Load libraries

library(methylKit)
library(ggplot2)

# function to get control samples, treatment samples, and treatment vector
# for use with methylkit's "reorganize" function
# input is:
#   a design_df with all comparisons (assumes specific column names and format)
#   mycomp variable which is a specific name of a comparison
get_reorg_vars <- function(design_df, mycomp) {

  # subset df to just comparison of interest
  df <- design_df[design_df$comparison ==  mycomp, ]

  # get control and treatment samples
  control_samps <- df$sample_name[df$comparison_designation == "control"]
  treat_samps <- df$sample_name[df$comparison_designation == "treatment"]

  # define treatment vector
  treatment_vec <- c(rep(0, length(control_samps)), rep(1, length(treat_samps)))

  # create a list to hold output vectors
  mylist <- list("control_samples" = control_samps,
    "treatment_samples" = treat_samps,
    "treatment_vec" = treatment_vec)
    
  return(mylist)
}

#-----#

# check for covariates
covariate_check <- function(design_df) {
  if (ncol(design_df) == 5) {
    print("No covariates present in design df.")
    return(FALSE)
  } else {
    print(paste0("Number of covariate columns found: ", ncol(design_df) - 5))
    return(TRUE)
  }
}

#-----#

# create covariate df
create_covar <- function(design_df, mycomp) {
  # subset df to just comparison of interest
  df <- design_df[design_df$comparison ==  mycomp, ]

  # subset to covar columns only
  df_covar <- df[, 6:ncol(df), drop = FALSE]

  # convert character columns to factor
  df_covar[] <- lapply(df_covar, function(col) {
    if (is.character(col)) as.factor(col) else col
  })

  return(df_covar)
}

#-----#

# obtain comparison name
get_comp_name <- function(design_df, mycomp) {
  # subset df to just comparison of interest
  df <- design_df[design_df$comparison ==  mycomp, ]

  return(unique(df$comparison_name))
}

#-----#
# Previously used methylKit DMR functions

############--------------------------------------#
# FUNCTION # to calculate DMRs and export results #
############--------------------------------------#
# will utilize methylKit's "MN" overdispersion and "Chisq" test
calc.DMRs <- function(my.meth, covariate=NULL, overdispersion="MN", test="Chisq", comparison, meth.diff=10, qval=0.1, type="DMR", mc=8) {
  options(scipen=999)
  
  print(paste0("Calculating DMRs for ", comparison))
  # Calculate DMRs: Overdispersion:YES, Test:Chisq
  myDiff <- calculateDiffMeth(my.meth,
                              covariates=covariate,
                              overdispersion=overdispersion,
                              test=test,
                              mc.cores=mc
  )
  # convert results to a data frame
  myDiff.df <- getData(myDiff)

  print(paste0("nrow of unfiltered diff results is: ", nrow(myDiff.df)))

  # export full results df (all regions tested for diff analysis)
  write.table(myDiff.df,
    paste0(comparison, ".fullres", type, "s.txt"),
    sep="\t",
    col.names=TRUE,
    row.names=FALSE,
    quote=FALSE
  )
  
  # Plot p-value distribution
  png(paste0(comparison, ".pval_dist.png"))
  hist(myDiff$pvalue)
  dev.off()

  # Subset to significant DMRs
  print(paste0("Number of sig ", type, "s: ", nrow(getMethylDiff(myDiff, difference=meth.diff, qvalue=qval))))
  myDiff.sig <- getMethylDiff(myDiff, difference=meth.diff, qvalue=qval)
  # convert to data frame
  myDiff.sig.df <- getData(myDiff.sig)
  print(paste0("nrow of myDiff.sig.df is: ", nrow(myDiff.sig.df)))

  if (nrow(myDiff.sig.df) > 0) {
    #----------------------------#
    # HYPO and HYPER SEPARATIONS #
    #----------------------------#

    # separate hyper and hypo sig DMRs
    myDiff.hyper <- myDiff.sig.df[myDiff.sig.df$meth.diff > 0, ]
    myDiff.hypo <- myDiff.sig.df[myDiff.sig.df$meth.diff < 0, ]

    print(paste0("number of hyper ", type, "s: ", nrow(myDiff.hyper)))
    print(paste0("number of hypo ", type, "s: ", nrow(myDiff.hypo)))

    # df for pie chart of hyper and hypo sig DMRs
    pie.df <- data.frame(class=c("hyper", "hypo"),
                         percent=c(nrow(myDiff.hyper) / nrow(myDiff.sig.df),
                                   nrow(myDiff.hypo) / nrow(myDiff.sig.df))
    )
    pie.df$label=paste0(round(as.numeric(pie.df$percent)*100, digits=0), '%')

    # Calculate midpoints for each segment
    pie.df$midpoint <- cumsum(pie.df$percent) - pie.df$percent/2

    # define blank theme for prettier pie chart
    blank_theme <- theme_minimal()+
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
    )

    # plot ratio of hyper vs hypo sig DMRs
    myplot <- ggplot(pie.df, aes(x="", y=percent, fill=class)) +
               geom_bar(width = 1, stat = "identity") +
               coord_polar("y", start=0) +
               blank_theme +
               theme(axis.text.x=element_blank()) +
	       geom_text(aes(x = 1, y = midpoint, label = label), size = 6) +
               #geom_text(aes(x=1, y=cumsum(percent/sum(percent)) - percent/sum(percent)/2,
               #  label=label[c(2,1)]), size=6) +  #label = paste(round(percent/sum(percent)*100),"%")), size = 6)
               ggtitle(paste0("Ratio of Hyper and Hypo Methylated sig ", type, "s\ntotal: ", nrow(myDiff.sig.df)))
    ggsave(filename=paste0(comparison, ".hyper.hypo.sig.", type, ".pieChart.png"))

    # -------------------------------#
    # continue with both hyper and hypo export

    # add DMR_ID column
    myDiff.sig.df$unq_ID <- paste(type, seq(1:nrow(myDiff.sig.df)), sep='_')

    # reorder columns
    myDiff.sig.df <- myDiff.sig.df[ ,c(8,1,2,3,5,6,7)]

    # export sig results
    write.table(myDiff.sig.df,
                paste0(comparison, ".sig", type, "s.txt"),
                sep="\t",
                col.names=TRUE,
                row.names=FALSE,
                quote=FALSE
    )

    #--------------#
    # sig dmr beds #
    #--------------#

    # BED file containing all sig DMRs and ensembl contig format
    # subset sig DMR results to a bed file for annotation
    res <- myDiff.sig.df[ ,c(2,3,4,1)]
    res$start <- res$start - 1

    # export sig DMR bed for annotation
    write.table(res,
                paste0(comparison, ".sig", type, "s.bed"),
                sep="\t",
                col.names=FALSE,
                row.names=FALSE,
                quote=FALSE
    )

    # 3 BED files for GREAT
    # UCSC contig format, 1 bed file of hyper only, 1 of hypo only, 1 of both
    # check if any seqnames begin with chr. If none do, assume not UCSC format and add "chr"
    if (!any(grepl("^chr", res$chr))) {
      print("Did not find UCSC format. Appending chr to seqnames.")
      res$chr <- paste0("chr", res$chr)
      ucsc_format <- FALSE
    } else {
      print("Data seqnames already appear to be UCSC format.")
      ucsc_format <- TRUE
    }
    write.table(res,
                paste0(comparison, ".sig", type, "s.both.GREAT.bed"),
                sep="\t",
                col.names=FALSE,
                row.names=FALSE,
                quote=FALSE
    )

    my.hyper <- myDiff.sig.df[myDiff.sig.df$meth.diff > 0, c(2,3,4,1)]
    # check if there are 0 sig hyper DMRs
    if (nrow(my.hyper) > 0) {
      my.hyper$start <- my.hyper$start - 1
      if (!ucsc_format) {
        my.hyper$chr <- paste0("chr", my.hyper$chr)
      }
      write.table(my.hyper,
                  paste0(comparison, ".sig", type, "s.hyper.GREAT.bed"),
                  sep="\t",
                  col.names=FALSE,
                  row.names=FALSE,
                  quote=FALSE
      )
    }
    else if (nrow(my.hyper) < 1) {
      print("Zero sig hyper DMRs for GREAT.")
    }

    my.hypo <- myDiff.sig.df[myDiff.sig.df$meth.diff < 0, c(2,3,4,1)]
    # check if there are 0 sig hypo DMRs
    if (nrow(my.hypo) > 0) {
      my.hypo$start <- my.hypo$start - 1
      if (!ucsc_format) {
        my.hypo$chr <- paste0("chr", my.hypo$chr)
      }
      my.hypo$chr <- paste0("chr", my.hypo$chr)
      write.table(my.hypo,
                  paste0(comparison, ".sig", type, "s.hypo.GREAT.bed"),
                  sep="\t",
                  col.names=FALSE,
                  row.names=FALSE,
                  quote=FALSE
      )
    }
    else if (nrow(my.hypo) < 1) {
      print("Zero sig hypo DMRs for GREAT.")
    }
  }
  else {
    print(paste0("There are zero ", type, "s passing signficance thresholds."))
  }
  return(myDiff.sig.df)
}

#-----#

############---------------------------------------------------------#
# FUNCTION # to make a sig DMR bed file track from results dataframe #
############---------------------------------------------------------#

makeBED <- function(res, comparison, type="DMR") {
  # check if any seqnames begin with chr. If none do, assume not UCSC format and add "chr"
  if (!any(grepl("^chr", res$chr))) {
    res$chr <- paste0("chr", res$chr)
  }

  res2 <- res[ ,c(2,3,4,1,7)]
  res2$start <- res2$start - 1
  res2$chr <- paste0("chr", res2$chr)

  res2[ ,6] <- "."
  res2[ ,7] <- res2$start
  res2[ ,8] <- res2$end
  res2[ ,9] <- ifelse(res2[ ,5] > 0, '255,0,0', ifelse(res2[ ,5] < 0, '0,0,255', '0,0,0'))

  # add track line to exported file
  cat(paste0("track type=bed name=", comparison, " ", "itemRgb=On"), file=paste0(comparison, ".sig", type, "track.bed"))
  cat("\n", file=paste0(comparison, ".sig", type, "track.bed"), append=TRUE)

  # export bed info
  write.table(res2,
              paste0(comparison, ".sig", type, "track.bed"),
              sep="\t",
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE,
              append=TRUE
  )
}

#-----#

run_design_comps <- function(meth, design_df, meth_diff = 10, qval_cutoff = 0.1, dm_type = "DMR") {
  for (i in 1:length(unique(design_df$comparison))) {
    mycomp <- unique(design_df$comparison)[i]
    reorg_var_list <- get_reorg_vars(design_df, mycomp)

    # reorganize
    meth_reorg <- reorganize(meth,
      sample.ids = c(reorg_var_list$control_samples, reorg_var_list$treatment_samples),
      treatment = reorg_var_list$treatment_vec)

    # check for covariates
    covars <- covariate_check(design_df)
    # if covariates are found, create covariate df
    if (covars) {
      print("Creating covariates df.")
      covar_df <- create_covar(design_df, mycomp)
    } else {
      covar_df <- NULL
    }

    # Obtain comparison name
    comp_name <- get_comp_name(design_df, mycomp)  

    # Execute calc.DMRs function
    # exports:
    #   1. Sig DMR results txt file
    #   2. Sig DMR results bed file (ensembl contigs)
    #   3. Hyper vs Hypo DMR pie chart
    #   4. Both GREAT bed file
    #   5. Hyper GREAT bed file
    #   6. Hypo GREAT bed file
    res <- calc.DMRs(meth_reorg,
      covariate = covar_df,
      overdispersion = "MN",
      test = "Chisq",
      comparison = comp_name,
      meth.diff = meth_diff,
      qval = qval_cutoff,
      type = dm_type,
      mc=1)

    # Make bed file track
    makeBED(res, comp_name, dm_type)
  }
}



#result_df <- design_df %>%
#  group_by(comparison, comparison_designation) %>%
#  filter(comparison_designation %in% c("control", "treatment")) %>%
#  summarize(sample_names = paste(sample_name, collapse = ", "))
#  #summarize(sample_names = c(paste(sample_name, sep=",")))
