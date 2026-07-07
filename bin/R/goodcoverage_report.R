# goodcoverage_report.R
# Computes Good's coverage statistics for all samples, saves a PDF of
# diagnostic plots, a TSV of the coverage compilation, and a summary text
# file with the mean ± SD coverage for the minimum singleton threshold.
#
# Author: Paul Dou (adapted by Jeremy Tournayre)
#
# Arguments:
#   args[1] - output directory
#   args[2] - feature table (.qza)

rm(list = ls())

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(gridExtra)
    library(vegan)
    library(cluster)
    library(fossil)
    # Load R library path from environment (set in config.sh).
    user_lib <- Sys.getenv("PIPELINE_R_LIBS")
    if (nzchar(user_lib)) .libPaths(user_lib)
    library(qiime2R)
})

args <- commandArgs(trailingOnly = TRUE)

output_dir <- args[1]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

raw_counts <- read_qza(args[2])

# Compute Good's coverage for a given singleton threshold.
goodcov <- function(x, thr) {
    if (thr == 0) stop("thr argument must be > 0")
    total_count    <- colSums(x)
    Ngleton_bool   <- x * (x > 0 & x <= thr)
    Ngleton_total  <- colSums(Ngleton_bool)
    goodscoverage  <- 100 * (1 - (Ngleton_total / total_count))
    return(data.frame(total = total_count, totalsing = Ngleton_total,
                      good = goodscoverage, sample = names(total_count)))
}

sample_selection <- function(raw_counts,
                             thr                     = 99.5,
                             max_Ngleton             = 100,
                             thr_outliers_goodcov    = 2,
                             thr_outliers_total      = 1.5,
                             upper_bound_low_count   = 50000,
                             res.subset              = as.logical()) {

    raw_counts_subset <- raw_counts
    if (length(res.subset) != 0) {
        if (any(res.subset)) raw_counts_subset <- raw_counts[, which(res.subset)]
    }

    cat("\n\nNumber of selected samples", dim(raw_counts_subset)[2], "\n\n")

    df_sample_info <- data.frame(total = colSums(raw_counts_subset))

    mincounts     <- raw_counts_subset
    mincounts[mincounts == 0] <- NA
    res.goodcov   <- goodcov(x = raw_counts_subset, thr = min(mincounts, na.rm = TRUE))
    min_mincounts <- min(mincounts, na.rm = TRUE)

    p1 <- ggplot(df_sample_info, aes(x = total)) +
        geom_density() + theme_light() +
        labs(caption = "Density of total counts across all samples.")

    p2 <- ggplot(df_sample_info, aes(x = total)) +
        geom_density() + xlim(c(0, upper_bound_low_count)) + theme_light() +
        labs(caption = "Density of total counts (lower range, below 50,000).")

    p3 <- ggplot(res.goodcov, aes(y = good, x = total)) +
        geom_point() + theme_light() + ggtitle(paste(min_mincounts, "tons", sep = "")) +
        labs(caption = "Good's coverage vs. total count.")

    p4 <- ggplot(res.goodcov, aes(y = good, x = total)) +
        geom_point() + xlim(c(0, upper_bound_low_count)) + theme_light() +
        ggtitle(paste(min_mincounts, "tons", sep = "")) +
        labs(caption = "Good's coverage vs. total count (lower range).")

    p5 <- ggplot(res.goodcov, aes(x = good)) +
        geom_density() + theme_light() +
        labs(caption = "Density of Good's coverage across all samples.")

    res.goodcov$richness <- colSums(raw_counts_subset != 0)

    res.goodcov.total <- NULL
    for (i in min_mincounts:max_Ngleton) {
        temp_goodcov <- goodcov(raw_counts_subset, i)
        temp_goodcov$thr      <- i
        temp_goodcov$richness <- colSums(raw_counts_subset != 0)
        temp_goodcov$outliers <- temp_goodcov$good  < quantile(temp_goodcov$good,  0.25) - thr_outliers_goodcov * IQR(temp_goodcov$good) &
                                 temp_goodcov$total < quantile(temp_goodcov$total, 0.25) - thr_outliers_total   * IQR(temp_goodcov$total)
        res.goodcov.total <- rbind(res.goodcov.total, temp_goodcov)
    }

    p6 <- ggplot(res.goodcov.total[res.goodcov.total$thr == min_mincounts, ],
                 aes(y = good, x = total)) +
        geom_point(aes(colour = outliers)) + theme_light() +
        ggtitle(paste(min_mincounts, "tons", sep = "")) +
        labs(caption = "Good's coverage vs. total count.")

    p7 <- ggplot(res.goodcov.total, aes(y = richness, x = total)) +
        geom_point(aes(colour = outliers)) + theme_light() +
        labs(caption = "Richness vs. total count.")

    p8 <- ggplot(res.goodcov.total, aes(x = thr, y = good, group = sample)) +
        geom_line(aes(colour = outliers)) + theme_light() +
        xlab("N-gleton") + ylab("Good's coverage") + theme(legend.position = "top") +
        labs(caption = "Good's coverage vs. rare-count threshold.")

    custom_theme <- theme(
        legend.position = "top",
        legend.text     = element_text(size = 8),
        axis.title      = element_text(size = 8),
        axis.text       = element_text(size = 8),
        plot.title      = element_text(size = 8)
    )

    plots <- list(
        p1 + custom_theme, p2 + custom_theme, p3 + custom_theme, p4 + custom_theme,
        p5 + custom_theme, p6 + custom_theme, p7 + custom_theme, p8 + custom_theme
    )

    pdf(file.path(output_dir, "goodcov_plots.pdf"), width = 10, height = 8)
    grid.arrange(grobs = plots, ncol = 2)
    dev.off()

    return(list(
        sample_out           = NULL,
        info                 = df_sample_info,
        compilation_coverage = res.goodcov.total
    ))
}

res.sel <- sample_selection(
    raw_counts            = raw_counts$data,
    thr                   = 99.98,
    max_Ngleton           = 10,
    thr_outliers_goodcov  = 2,
    thr_outliers_total    = 1.5,
    upper_bound_low_count = 50000,
    res.subset            = NULL
)

write.table(
    res.sel$compilation_coverage,
    file      = file.path(output_dir, "goodcov_compilation_coverage.tsv"),
    sep       = "\t", row.names = FALSE, quote = FALSE
)

good_vals    <- res.sel$compilation_coverage[res.sel$compilation_coverage$thr == min(res.sel$compilation_coverage$thr), ]$good
XX           <- mean(good_vals, na.rm = TRUE)
yy           <- sd(good_vals,   na.rm = TRUE)
summary_text <- sprintf(
    "Mean Good's coverage for %.0f gleton was of %.4f \u00b1 %.4f",
    min(res.sel$compilation_coverage$thr), XX, yy
)
writeLines(summary_text, con = file.path(output_dir, "goodcov_Mean.txt"))