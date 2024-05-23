
# Filter a bsseq object to contain CpGs with at least min_cov
# coverage in at least perc_samp percent of Samples
filter_bsseq <- function(bs, min_cov=1, perc_samp=0.75) {
  print(paste0("Filtering to include CpGs with at least ", min_cov,
    "X coverage in at least ", perc_samp*100, "% of samples."))
  # get a binary matrix of coverage greater than threshold
  covSample <- (getCoverage(bs) >= min_cov) %>% DelayedMatrixStats::rowSums2()
  # subset
  bs.filtered <- bs[covSample >= (perc_samp * ncol(bs)), ]
  print(paste0("Returning ", nrow(bs.filtered), " CpGs."))
  return(bs.filtered)
}

run_design_comps <- function(bs, design_df, min_cov=1, perc_samp=0.75, cutoff=0.05, test_cov, match_cov, num_cores) {
  for (i in 1:length(unique(design_df$comparison))) {
    mycomp <- unique(design_df$comparison)[i]
    # subset df to just comparison of interest
    df <- design_df[design_df$comparison ==  mycomp, ]

    print(paste0("Beginning comparison: ", unique(df$comparison_name)))
    # Print match covariate
    if (is.null(match_cov)) {
      print(paste0("Match covariate is set to: ", deparse(match_cov)))
    } else {
      print(paste0("Match covariate is set to: ", match_cov))
    }
    
    # get control and treatment samples
    control_samps <- df$sample_name[df$comparison_designation == "control"]
    treat_samps <- df$sample_name[df$comparison_designation == "treatment"]

    # subset bsseq object to only samples from 2 groups used in comparison
    sample.idx <- which(pData(bs)$sample_name %in% c(control_samps, treat_samps))
    bs.sub <- bs[, sample.idx]

    # Filter CpGs
    bs.filtered <- filter_bsseq(bs.sub, min_cov=min_cov, perc_samp=perc_samp)
    print(pData(bs.filtered))
    # save filtered bs object
    saveRDS(bs.filtered, paste0("RData/", unique(df$comparison_name), ".filtered_bsseq.RDS"))

    # Perform Diff Comparison via dmrseq
    regions <- dmrseq(bs = bs.filtered,
                  cutoff = cutoff,
                  testCovariate = test_cov,
		  matchCovariate = match_cov,
                  BPPARAM = BiocParallel::MulticoreParam(workers = num_cores))

    saveRDS(regions, paste0("RData/", unique(df$comparison_name), ".regions.RDS"))
    # or continue with plotting?
    # plotting as separate function?
  }
}
