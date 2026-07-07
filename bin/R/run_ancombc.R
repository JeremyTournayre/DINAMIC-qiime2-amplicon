# run_ancombc.R
# Runs differential abundance analysis at each taxonomic level using
# ANCOM-BC2, DESeq2, LEfSe and IndicSpecies.  A patched local copy of
# ANCOMBC is sourced from the directory pointed to by the ANCOMBC_DIR
# environment variable (set in config.sh).
#
# Author: Jeremy Tournayre
#
# Arguments:
#   args[1] - comparison directory (contains table.qza/rooted-tree.qza one level up)
#   args[2] - taxonomy artifact (.qza)
#   args[3] - metadata column for group comparison
#   args[4] - "yes" to produce boxplots, anything else to skip
#   args[5] - prevalence filter "N,frac" (e.g. "100,0.5")

rm(list = ls())
options(warn = -1)

# Fixed seed for reproducibility of multipatt's (IndicSpecies) random
# permutations below.
set.seed(42)

# Load R library path from environment (set in config.sh).
user_lib <- Sys.getenv("PIPELINE_R_LIBS")
if (nzchar(user_lib)) .libPaths(user_lib)

library(cowplot)

# Source the local patched copy of ANCOMBC.
ancombc_dir <- Sys.getenv("ANCOMBC_DIR")
if (!nzchar(ancombc_dir)) stop("ANCOMBC_DIR is not set. Check config.sh.")
tmp_getwd   <- getwd()
setwd(ancombc_dir)
script_files <- list.files(path = "R", pattern = "\\.R$", full.names = TRUE)
for (file in script_files) source(file)
setwd(tmp_getwd)

library(lefser)
library(Cairo)
library(indicspecies)
library(DESeq2)
library(qiime2R)
library(data.table)
library(doSNOW)
library(phyloseq)
library(stringr)
library(doRNG)
suppressMessages(library(microbiome))

args <- commandArgs(trailingOnly = TRUE)

cl <- parallel::makeCluster(1)
registerDoSNOW(cl)

physeq <- qza_to_phyloseq(
    features  = paste(args[1], "/../table.qza",       sep = ""),
    tree      = paste(args[1], "/../rooted-tree.qza",  sep = ""),
    args[2],
    metadata  = "../mapping.tsv"
)

var <- args[3]

meta <- sample_data(physeq)
meta[[var]] <- as.character(meta[[var]])
meta[[var]][is.na(meta[[var]])] <- "NA"
meta[[var]] <- factor(meta[[var]])
sample_data(physeq) <- meta


taxa <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Boxplot function: coloured by group variable.
toplot2 <- function() {
    p <<- ggplot(data, aes_string(x = args[3], y = "name", fill = args[3])) +
        geom_boxplot() +
        geom_jitter(color = "black", height = 0, size = 0.4, alpha = 0.9) +
        theme(
            legend.position = "none",
            axis.text.x     = element_text(angle = 45, hjust = 1)
        ) +
        ylab(title_ylab)
}

groups    <- unique(sample_data(physeq)[, colnames(sample_data(physeq)) == args[3]])[[1]]
groups    <- factor(ifelse(is.na(as.character(groups)), "NA", as.character(groups)))
dir_out   <- paste(args[1], '/', args[5], sep = '')
if (!file.exists(dir_out)) dir.create(dir_out)

i_group  <- 1
dir_out2 <- paste(dir_out, "/ANCOM_BC_ref_is_", levels(groups)[i_group], sep = "")
if (file.exists(dir_out2)) unlink(dir_out2, recursive = TRUE)
dir.create(dir_out2)

if (exists("res"))         remove(res)
if (exists("respairwise")) remove(respairwise)

if (nrow(tax_table(physeq)) < 50) {
    warning_size_taxo <- paste(
        "The number of taxa is too small for reliable ANCOMBC results.",
        "Typically > 50 taxa are required; current count:",
        nrow(tax_table(physeq)))
    outfile <- paste(dir_out2, "/warning_size_taxo", sep = "")
    write.table("warning_size_taxo", outfile, sep = "\t", row.names = FALSE)
}

meta_df   <- as.data.frame(sample_data(physeq))
group_col <- args[3]

# qza_to_phyloseq() reads mapping.tsv with its own defaults, which treat the
# literal text "NA" as a genuine missing value (R's standard convention).
# In this pipeline's metadata, "NA" can also be a deliberate, meaningful
# category value (e.g. "no particular treatment"), not necessarily a
# missing measurement. Re-read the comparison column directly as raw text
# and restore any value that qza_to_phyloseq() turned into a real NA back
# to the literal string "NA", so it is treated as its own group rather than
# silently dropped or breaking group-membership tests further down.
raw_mapping <- read.table("../mapping.tsv", header = TRUE, sep = "\t",
                          row.names = 1, comment.char = "", colClasses = "character",
                          check.names = FALSE, na.strings = character(0))
# Drop the QIIME2 "#q2:types" row, if present (its row name starts with
# "#q2:types" once read as the first data row here since header=TRUE only
# consumes the very first line).
if (nrow(raw_mapping) > 0 && grepl("^#q2:types", rownames(raw_mapping)[1])) {
    raw_mapping <- raw_mapping[-1, , drop = FALSE]
}
if (group_col %in% colnames(raw_mapping)) {
    raw_values <- raw_mapping[rownames(meta_df), group_col]
    meta_df[[group_col]] <- ifelse(is.na(meta_df[[group_col]]) & raw_values == "NA",
                                   "NA", as.character(meta_df[[group_col]]))
    meta_df[[group_col]] <- factor(meta_df[[group_col]])
}

for (i in 2:length(taxa)) {
    aggregate_taxa_physeq <- aggregate_taxa(physeq, taxa[i])

    taxonomy_df <- as.data.frame(tax_table(aggregate_taxa_physeq))
    if ("unique" %in% colnames(taxonomy_df))
        taxonomy_df <- taxonomy_df[, !colnames(taxonomy_df) %in% "unique"]
    taxonomy_df$taxon_full <- apply(taxonomy_df, 1, function(x) paste(na.omit(x), collapse = "|"))

    otu_mat       <- as.data.frame(otu_table(aggregate_taxa_physeq))
    otu_mat$taxon_full <- taxonomy_df$taxon_full

    groups_local <- unique(meta_df[[group_col]])

    # Prevalence filtering: keep taxa present in at least frac of samples in each group.
    taxa_to_keep <- c()
    strsplit_args <- as.numeric(strsplit(args[5], ",")[[1]])

    for (group in groups_local) {
        samples_in_group  <- rownames(meta_df[meta_df[[group_col]] == group, ])
        otu_subset        <- otu_mat[, samples_in_group, drop = FALSE]
        taxa_presence     <- rowSums(otu_subset >= strsplit_args[1])
        threshold         <- ceiling(strsplit_args[2] * length(samples_in_group))
        taxa_to_keep      <- unique(c(taxa_to_keep, otu_mat$taxon_full[taxa_presence >= threshold]))
    }

    tax_table_df <- as.data.frame(tax_table(physeq))
    if ("unique" %in% colnames(tax_table_df))
        tax_table_df <- tax_table_df[, !colnames(tax_table_df) %in% "unique"]
    level_index          <- which(colnames(tax_table_df) == taxa[i])
    tax_table_df$taxon_full <- apply(tax_table_df[, 1:level_index, drop = FALSE], 1,
                                     function(x) paste(na.omit(x), collapse = "|"))

    asv_to_keep        <- rownames(tax_table_df[tax_table_df$taxon_full %in% taxa_to_keep, , drop = FALSE])
    physeq_asv_filtered <- prune_taxa(asv_to_keep, physeq)
    physeq_filtered_tmp <- aggregate_taxa(physeq_asv_filtered, taxa[i])

    tax_asv_filtered <- as.data.frame(tax_table(physeq_filtered_tmp))
    if ("unique" %in% colnames(tax_asv_filtered))
        tax_asv_filtered <- tax_asv_filtered[, !colnames(tax_asv_filtered) %in% "unique"]
    level_index <- which(colnames(tax_asv_filtered) == taxa[i])
    tax_asv_filtered$taxon_full <- apply(tax_asv_filtered[, 1:level_index, drop = FALSE], 1,
                                         function(x) paste(na.omit(x), collapse = "|"))

    message_DESEQ2 <- c()
    physeq_filtered <- aggregate_taxa(physeq_asv_filtered, taxa[i])

    if (nrow(tax_table(physeq_filtered)) > 1) {
        otu_physeq_filtered_aggr_matr <- as(otu_table(physeq_filtered), "matrix")
        # mapping <- data.frame(sample_data(physeq_filtered))
        # Use meta_df's already-corrected group column (see the "NA" as a
        # literal category fix above) rather than re-reading it from
        # physeq_filtered's sample_data, which would reintroduce the same
        # literal-"NA"-turned-into-real-NA problem for every DESeq2/
        # IndicSpecies/LEfSe comparison below.

        mapping <- as(sample_data(physeq_filtered), "data.frame")

        design_formula      <- as.formula(paste("~", args[3]))

        cat("[", taxa[i], "] group counts for DESeq2 design: ",
            paste(names(table(mapping[, args[3]])), table(mapping[, args[3]]), sep = "=", collapse = ", "),
            "\n", sep = "")

        # DESeq2 analysis. Both DESeqDataSetFromMatrix() and DESeq() are
        # wrapped: DESeqDataSetFromMatrix() itself fails (not just DESeq())
        # when, after filtering, only one group value remains for this
        # taxonomic level (e.g. "TemperatureConservation" reduced to a
        # single value among the samples retained) — a real, expected
        # per-level data limitation, not a bug, and should be skipped like
        # any other per-level DESeq2 failure rather than halting the script.
        dds <- tryCatch({

            dds_tmp <- DESeqDataSetFromMatrix(
                countData = otu_physeq_filtered_aggr_matr,
                colData   = mapping,
                design    = design_formula
            )

            DESeq(dds_tmp)

        }, error = function(e) {

            message("DESeq failed, switching to poscounts: ", e$message)

            dds_tmp <- DESeqDataSetFromMatrix(
                countData = otu_physeq_filtered_aggr_matr,
                colData   = mapping,
                design    = design_formula
            )

            dds_tmp <- estimateSizeFactors(dds_tmp, type = "poscounts")

            DESeq(dds_tmp)
        })

        if (!is.null(dds)) {
            results         <- results(dds)
            rownames(results) <- tax_asv_filtered$taxon_full
            res_cbind       <- cbind(rownames(results), taxa[i], results)
            if (exists("deseq2_results_global")) {
                deseq2_results_global <- rbind(deseq2_results_global, res_cbind)
            } else {
                deseq2_results_global <- res_cbind
            }

            comparisons        <- combn(levels(mapping[, args[3]]), 2, simplify = FALSE)
            pairwise_results_d <- lapply(comparisons, function(x) results(dds, contrast = c(args[3], x[1], x[2])))
            names(pairwise_results_d) <- sapply(comparisons, function(x) paste(x, collapse = "_vs_"))

            if (length(comparisons) > 1) {
                formatted_results <- list()
                for (comparison in names(pairwise_results_d)) {
                    grps   <- unlist(strsplit(comparison, "_vs_"))
                    res_tmp <- pairwise_results_d[[comparison]]
                    res_tmp$features <- tax_asv_filtered$taxon_full
                    res_tmp$group1   <- grps[1]
                    res_tmp$group2   <- grps[2]
                    res_tmp$level    <- taxa[i]
                    res_tmp <- res_tmp[, c("features", "level", "group1", "group2",
                                           colnames(res_tmp)[1:(ncol(res_tmp) - 4)])]
                    formatted_results[[comparison]] <- res_tmp
                }
                new_pw_deseq2 <- do.call(rbind, formatted_results)
                if (exists("deseq2_results_pairwise")) {
                    deseq2_results_pairwise <- rbind(deseq2_results_pairwise, new_pw_deseq2)
                } else {
                    deseq2_results_pairwise <- new_pw_deseq2
                }
            }
        }

        # IndicSpecies analysis.
        cat("[", taxa[i], "] group counts for multipatt (IndicSpecies): ",
            paste(names(table(mapping[, args[3]])), table(mapping[, args[3]]), sep = "=", collapse = ", "),
            " | n samples in otu matrix: ", ncol(otu_physeq_filtered_aggr_matr),
            " | n samples in mapping: ", nrow(mapping),
            "\n", sep = "")
        results_ind    <- multipatt(t(otu_physeq_filtered_aggr_matr), mapping[, args[3]], func = "IndVal.g")
        results_ind    <- results_ind$sign
        rownames(results_ind) <- tax_asv_filtered$taxon_full
        res_cbind_ind  <- cbind(features = rownames(results_ind), level = taxa[i], results_ind)
        if (exists("IndicSpecies_results_global")) {
            IndicSpecies_results_global <- rbind(IndicSpecies_results_global, res_cbind_ind)
        } else {
            IndicSpecies_results_global <- res_cbind_ind
        }

        groups_is  <- unique(mapping[, args[3]])
        formatted_is <- list()
        if (length(groups_is) > 2) {
            for (ii in 1:(length(groups_is) - 1)) {
                for (j in (ii + 1):length(groups_is)) {
                    g1           <- groups_is[ii]
                    g2           <- groups_is[j]
                    subset_idx   <- mapping[, args[3]] %in% c(g1, g2)
                    subset_data  <- t(otu_physeq_filtered_aggr_matr)[subset_idx, ]
                    subset_grps  <- as.factor(as.character(mapping[subset_idx, args[3]]))
                    res_pw_is    <- multipatt(subset_data, subset_grps, func = "IndVal.g")
                    pvals        <- res_pw_is$sign
                    rownames(pvals) <- tax_asv_filtered$taxon_full
                    colnames(pvals)[colnames(pvals) == paste0("s.", g1)] <- "s.group1"
                    colnames(pvals)[colnames(pvals) == paste0("s.", g2)] <- "s.group2"
                    pvals$group1 <- g1; pvals$group2 <- g2; pvals$level <- taxa[i]
                    pvals        <- pvals[, c("level", "group1", "group2", "s.group1", "s.group2",
                                              setdiff(colnames(pvals), c("level", "group1", "group2", "s.group1", "s.group2")))]
                    pvals        <- cbind(features = rownames(pvals), pvals)
                    formatted_is[[paste(g1, g2, sep = "_vs_")]] <- pvals
                }
            }
            new_is <- do.call(rbind, formatted_is)
            rownames(new_is) <- NULL
            if (exists("IndicSpecies_results_pairwise")) {
                IndicSpecies_results_pairwise <- rbind(IndicSpecies_results_pairwise, new_is)
            } else {
                IndicSpecies_results_pairwise <- new_is
            }
        }

        # LEfSe pairwise analysis.
        pairwise_lefse_results <- list()
        for (ii in 1:(length(groups_is) - 1)) {
            for (j in (ii + 1):length(groups_is)) {
                g1           <- groups_is[ii]
                g2           <- groups_is[j]
                subset_idx   <- mapping[, args[3]] %in% c(g1, g2)
                subset_data  <- otu_physeq_filtered_aggr_matr[, subset_idx]
                subset_grps  <- as.factor(as.character(mapping[subset_idx, args[3]]))
                rownames(subset_data) <- tax_asv_filtered$taxon_full
                se <- SummarizedExperiment(
                    assays  = list(counts = subset_data),
                    colData = data.frame(class = subset_grps)
                )
                res_lefse <- lefser(se, class = "class",
                                    kruskal.threshold = 0.05, wilcox.threshold = 0.05)
                if (nrow(res_lefse) > 0) {
                    res_lefse$group1 <- g1
                    res_lefse$group2 <- g2
                    res_lefse$level  <- taxa[i]
                    pairwise_lefse_results[[paste(g1, g2, sep = "_vs_")]] <- res_lefse
                }
            }
        }
        new_lefse_f <- do.call(rbind, pairwise_lefse_results)
        new_lefse_f <- new_lefse_f[, c("features", "level", "group1", "group2", "scores")]
        if (exists("lefse_results_pairwise")) {
            lefse_results_pairwise <- rbind(lefse_results_pairwise, new_lefse_f)
        } else {
            lefse_results_pairwise <- new_lefse_f
        }
    }

    # ANCOM-BC2: add pseudocount to completely absent taxa to avoid structural zeros.
    pseudovalue <- 1
    otu_table_v <- otu_table(physeq_asv_filtered)

    for (group in groups_local) {
        samples_in_group <- rownames(meta_df[meta_df[[group_col]] == group, ])
        zero_taxa        <- rowSums(otu_table(physeq_filtered_tmp)[, samples_in_group]) == 0
        names(zero_taxa) <- tax_asv_filtered$taxon_full
        zero_asv_list    <- rownames(tax_table_df)[tax_table_df$taxon_full %in% names(zero_taxa)[zero_taxa]]
        zero_taxon_list  <- tax_table_df$taxon_full[tax_table_df$taxon_full %in% names(zero_taxa)[zero_taxa]]
        already_done     <- c()
        zz <- 0
        for (taxon in zero_taxon_list) {
            zz <- zz + 1
            if (taxon %in% already_done) next
            selected_asv                          <- zero_asv_list[zz]
            otu_table_v[selected_asv, samples_in_group] <- pseudovalue
            already_done <- c(already_done, taxon)
        }
    }
    physeq_asv_filtered@otu_table <- otu_table_v

    physeq_filtered <- aggregate_taxa(physeq_asv_filtered, taxa[i])

    if (nrow(tax_table(physeq_filtered)) > 1) {
        out <- suppressMessages(
            ancombc2(data = physeq_asv_filtered, fix_formula = args[3],
                     tax_level = taxa[i], p_adj_method = "BH",
                     global = TRUE, group = args[3], pairwise = TRUE)
        )

        abundance_mat    <- abundances(physeq_filtered)
        mean_abundances  <- sapply(groups_local, function(g) {
            samples_in_group <- rownames(meta_df[meta_df[[group_col]] == g, ])
            rowMeans(abundance_mat[, samples_in_group, drop = FALSE], na.rm = TRUE)
        }, simplify = TRUE)
        mean_abundances_df           <- as.data.frame(mean_abundances)
        colnames(mean_abundances_df) <- paste('Mean_', groups_local, sep = '')
        mean_abundances_df$Taxon     <- tax_asv_filtered$taxon_full

        if (length(out$res_global)) {
            res_cbind        <- cbind(out$res_global[, 1], taxa[i], out$res_global[, -1])
            respairwise_cbind <- cbind(out$res_pair[, 1],  taxa[i], out$res_pair[, -1])
            res_cbind        <- merge(res_cbind,        mean_abundances_df,
                                      by.x = "out$res_global[, 1]", by.y = "Taxon", all.x = TRUE)
            respairwise_cbind <- merge(respairwise_cbind, mean_abundances_df,
                                       by.x = "out$res_pair[, 1]",   by.y = "Taxon", all.x = TRUE)
        } else {
            res_cbind <- cbind(out$res[, 1], taxa[i], out$res[, -1])
            res_cbind <- merge(res_cbind, mean_abundances_df,
                               by.x = "out$res[, 1]", by.y = "Taxon", all.x = TRUE)
        }

        if (exists("res")) {
            res <- rbind(res, res_cbind)
            if (length(out$res_global)) respairwise <- rbind(respairwise, respairwise_cbind)
        } else {
            res <- res_cbind
            if (length(out$res_global)) respairwise <- respairwise_cbind
        }

        # Boxplots of bias-corrected abundances (only when requested).
        if (args[4] == "yes") {
            if (is.null(out$res_global)) {
                q_cols    <- grep("^q_", names(out$res), value = TRUE)
                q_col     <- setdiff(q_cols, "q_(Intercept)")
                out$res_global$q_val  <- out$res[[q_col]]
                out$res_global$taxon  <- out$res$taxon
            }

            samp_frac         <- out$samp_frac
            samp_frac[is.na(samp_frac)] <- 0
            log_obs_abn       <- log(abundances(physeq_asv_filtered) + 1)
            log_obs_abn_adj   <- t(t(log_obs_abn) - samp_frac)
            unlog_obs_abn_adj <- exp(log_obs_abn_adj) - 1

            physeq_log_obs_abn_adj   <- physeq_asv_filtered
            physeq_log_obs_abn       <- physeq_asv_filtered
            physeq_unlog_obs_abn_adj <- physeq_asv_filtered

            otu_table(physeq_log_obs_abn)       <- otu_table(log_obs_abn,       taxa_are_rows = TRUE)
            otu_table(physeq_log_obs_abn_adj)   <- otu_table(log_obs_abn_adj,   taxa_are_rows = TRUE)
            otu_table(physeq_unlog_obs_abn_adj) <- otu_table(unlog_obs_abn_adj, taxa_are_rows = TRUE)

            physeq_log_obs_abn_adj   <- aggregate_taxa(physeq_log_obs_abn_adj,   taxa[i])
            physeq_log_obs_abn       <- aggregate_taxa(physeq_log_obs_abn,       taxa[i])
            physeq_unlog_obs_abn_adj <- aggregate_taxa(physeq_unlog_obs_abn_adj, taxa[i])

            log_obs_abn_adj   <- otu_table(physeq_log_obs_abn_adj)
            log_obs_abn       <- otu_table(physeq_log_obs_abn)
            unlog_obs_abn_adj <- otu_table(physeq_unlog_obs_abn_adj)

            rownames(log_obs_abn_adj)   <- tax_asv_filtered$taxon_full
            rownames(log_obs_abn)       <- tax_asv_filtered$taxon_full
            rownames(unlog_obs_abn_adj) <- tax_asv_filtered$taxon_full

            tmp_aggregate <- abundances(physeq_filtered)
            rownames(tmp_aggregate) <- tax_asv_filtered$taxon_full

            k <- 0
            for (j in order(out$res_global$q_val)) {
                k    <- k + 1
                name <- out$res_global$taxon[j]
                qvalue <- if (length(out$res_global)) out$res_global$q_val[j] else out$res_global$q_val[[1]][j]

                dir_box <- if (qvalue > 0.05) "Abundance_NS_Globaltest" else "Abundance_Significant_Globaltest"
                dir_out3 <- paste(dir_out2, "/", dir_box, sep = "")
                if (!file.exists(dir_out3))            dir.create(dir_out3)
                if (!file.exists(paste0(dir_out3, "/")))       dir.create(paste0(dir_out3, "/"))
                if (!file.exists(paste0(dir_out3, "/", taxa[i]))) dir.create(paste0(dir_out3, "/", taxa[i]))

                if (all(rownames(log_obs_abn_adj) != name)) next

                tmp_log_adj   <- t(as.matrix(c(log_obs_abn_adj[rownames(log_obs_abn_adj)     == name, ])))
                tmp_unlog_adj <- t(as.matrix(c(unlog_obs_abn_adj[rownames(unlog_obs_abn_adj) == name, ])))
                tmp_log_raw   <- t(as.matrix(c(log_obs_abn[rownames(log_obs_abn)             == name, ])))
                tmp_raw       <- t(as.matrix(tmp_aggregate[rownames(tmp_aggregate)            == name, ]))

                features <- 1
                sample_data_corrected <- meta_df[rownames(sample_data(physeq)), , drop = FALSE]
                data <- cbind(features, tmp_log_adj[features, ], sample_data_corrected)
                data$features <- factor(data$features); colnames(data)[2] <- "name"
                title_ylab <- "Bias-adjusted log observed abundances"; toplot2(); p1_2 <- p

                data <- cbind(features, tmp_unlog_adj[features, ], sample_data_corrected)
                data$features <- factor(data$features); colnames(data)[2] <- "name"
                title_ylab <- "Bias-adjusted unlog observed abundances"; toplot2(); p2_2 <- p

                data <- cbind(features, tmp_log_raw[features, ], sample_data_corrected)
                data$features <- factor(data$features); colnames(data)[2] <- "name"
                title_ylab <- "Log observed abundances"; toplot2(); p3_2 <- p

                data <- cbind(features, tmp_raw[features, ], sample_data_corrected)
                data$features <- factor(data$features); colnames(data)[2] <- "name"
                title_ylab <- "Observed abundances"; toplot2(); p4_2 <- p

                a            <- plot_grid(p1_2, p2_2, p3_2, p4_2, ncol = 2, nrow = 2)
                title_string <- paste0(name, " (", taxa[i], ") qval: ", round(qvalue, digits = 4))
                title        <- ggdraw() +
                    draw_label(title_string, size = 10, x = 0, hjust = 0) +
                    theme(plot.margin = margin(0, 0, 0, 7))
                b <- plot_grid(title, a, ncol = 1, rel_heights = c(0.1, 1))

                taxon_short <- sub(".*\\|(.*)", "\\1", name)
                taxon_short <- substr(taxon_short, 1, 50)
                outfile_png <- paste0(dir_out3, "/", taxa[i], "/", k, "_", taxon_short, ".png")
                CairoPNG(outfile_png, width = 10, height = 10, units = "in", res = 300)
                print(b)
                dev.off()
            }
        }
    }
}

# Write global ANCOM-BC results.
colnames(res)[1] <- "name"
colnames(res)[2] <- "Stat on taxa"
if (exists("respairwise")) res <- res[order(as.numeric(res[, 5])), ]
write.table(res, paste0(dir_out2, "/ANCOM-BC.tsv"), sep = "\t", row.names = FALSE)

if (is.null(lefse_results_pairwise)) lefse_results_pairwise <- c("No significant features found.")

if (exists("respairwise")) {
    colnames(respairwise)[1] <- "name"
    colnames(respairwise)[2] <- "Stat on taxa"
    respairwise[, 3] <- as.character(respairwise[, 3])
    write.table(respairwise, paste0(dir_out2, "/ANCOM-BC_pairwise.tsv"), sep = "\t", row.names = FALSE)

    if (exists("deseq2_results_pairwise")) {
        colnames(deseq2_results_pairwise)[1] <- "name"
        colnames(deseq2_results_pairwise)[2] <- "Stat on taxa"
        write.table(deseq2_results_pairwise, paste0(dir_out, "/DESEQ2_pairwise.tsv"), sep = "\t", row.names = FALSE)
    } else {
        write.table(message_DESEQ2, paste0(dir_out, "/DESEQ2_pairwise.tsv"), sep = "\t", row.names = FALSE)
    }

    if (exists("IndicSpecies_results_pairwise")) {
        colnames(IndicSpecies_results_pairwise)[1] <- "name"
        colnames(IndicSpecies_results_pairwise)[2] <- "Stat on taxa"
        write.table(IndicSpecies_results_pairwise, paste0(dir_out, "/IndicSpecies_pairwise.tsv"), sep = "\t", row.names = FALSE)
    } else {
        write.table(c("Error."), paste0(dir_out, "/IndicSpecies_pairwise.tsv"), sep = "\t", row.names = FALSE)
    }

    if (!is.null(dim(lefse_results_pairwise))) {
        colnames(lefse_results_pairwise)[1] <- "name"
        colnames(lefse_results_pairwise)[2] <- "Stat on taxa"
    }
    write.table(lefse_results_pairwise, paste0(dir_out, "/lefse_pairwise.tsv"), sep = "\t", row.names = FALSE)
} else {
    if (!is.null(dim(lefse_results_pairwise))) {
        colnames(lefse_results_pairwise)[1] <- "name"
        colnames(lefse_results_pairwise)[2] <- "Stat on taxa"
    }
    write.table(lefse_results_pairwise, paste0(dir_out, "/lefse.tsv"), sep = "\t", row.names = FALSE)
}

if (exists("deseq2_results_global")) {
    colnames(deseq2_results_global)[1] <- "name"
    colnames(deseq2_results_global)[2] <- "Stat on taxa"
    write.table(deseq2_results_global, paste0(dir_out, "/DESEQ2.tsv"), sep = "\t", row.names = FALSE)
} else {
    write.table(message_DESEQ2, paste0(dir_out, "/DESEQ2.tsv"), sep = "\t", row.names = FALSE)
}

if (exists("IndicSpecies_results_global")) {
    colnames(IndicSpecies_results_global)[1] <- "name"
    colnames(IndicSpecies_results_global)[2] <- "Stat on taxa"
    write.table(IndicSpecies_results_global, paste0(dir_out, "/IndicSpecies.tsv"), sep = "\t", row.names = FALSE)
} else {
    write.table(c("Error."), paste0(dir_out, "/IndicSpecies.tsv"), sep = "\t", row.names = FALSE)
}

if (args[4] == "yes")
    write.table("done", paste0(dir_out, "/done_ANCOM-BC_boxplot_2"), sep = "\t", row.names = FALSE)

write.table("done", paste0(dir_out, "/done_DifferentialAbundance"), sep = "\t", row.names = FALSE)