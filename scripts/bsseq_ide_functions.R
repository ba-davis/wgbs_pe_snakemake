

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
