# diversity_vegan_tss.R
# Computes alpha and beta diversity statistics using vegan and phyloseq
# under three normalisations: non-normalised, rarefied, and TSS (total-sum
# scaling). Results are written to TSV files in a diversity_vegan_TSS/
# subdirectory.
#
# Author: Jeremy Tournayre
#
# Arguments:
#   args[1] - comparison directory (contains table.qza and rooted-tree.qza)
#   args[2] - metadata column for group comparison
#   args[3] - rarefaction depth
#   args[4] - QIIME2 artifact UUID (biom subfolder name)
#   args[5] - optional blocking factor: mapping.tsv column (e.g. Animal) or
#             "none" (default). Used only for the *_block_* / alpha_lm
#             outputs below; the historical one-factor outputs
#             (permanova_results_*, pairwise_results_*, kruskal_results_*,
#             ...) are unchanged because summary_results.pl parses them by
#             column position.
#
# Outputs added 2026-11 (reviewer requests: donor effect, compositionality):
#   - Aitchison distance (CLR + Euclidean, pseudocount 1) added to the
#     one-factor PERMANOVA lists, non_normalized table only (Aitchison is
#     scale invariant, so TSS would give the same result, and rarefying
#     before CLR only adds zeros).
#   - betadisper_<norm>.tsv: dispersion test (betadisper + permutest).
#   - permanova_block_<norm>.tsv: adonis2(dist ~ block + group,
#     by = "margin"), once with free permutations and once with
#     permutations restricted within blocks (or free permutations again,
#     labelled "fallback", when blocks allow < 100 distinct permutations,
#     see perm_scheme). by = "margin" is explicit:
#     vegan >= 2.7 defaults to by = NULL (omnibus model test only).
#   - pairwise_block_<norm>.tsv: same blocked model on each pair of groups
#     (only when the comparison has > 2 groups), BH-adjusted.
#   - alpha_lm_<norm>.tsv: alpha ~ block + group linear models, marginal
#     F tests (drop1), BH across indices for each term.

rm(list = ls())
options(warn = -1)

# Fixed seed for reproducibility of permutation-based tests below
# (adonis2, pairwise.adonis use random permutations; rarefy_even_depth
# already had its own explicit rngseed).
set.seed(42)

library(phyloseq)
library(ggplot2)
library(vegan)
library(pairwiseAdonis)
library(ape)

# Load R library path from environment (set in config.sh).
user_lib <- Sys.getenv("PIPELINE_R_LIBS")
if (nzchar(user_lib)) .libPaths(user_lib)

library(qiime2R)

args <- commandArgs(trailingOnly = TRUE)

biom_file   <- paste(args[1], '/table/', args[4], '/data/feature-table.biom', sep = "")
table       <- read_q2biom(biom_file)
list_sample <- colnames(table)

feature_table <- read_qza(paste(args[1], "/table.qza",       sep = ""))
otu_table     <- otu_table(feature_table$data, taxa_are_rows = TRUE)
phy_tree      <- read_qza(paste(args[1], "/rooted-tree.qza", sep = ""))$data

# Resolve multichotomies to avoid warnings from UniFrac.
new_tre <- ape::multi2di(phy_tree)

# mapping.tsv's "NA" values can be a deliberate, meaningful category (e.g.
# "no particular treatment"), not necessarily a missing measurement.
# na.strings = character(0) prevents read.table from silently turning the
# literal text "NA" into a real NA, which would later crash vegan/adonis2
# (na.fail on any NA in the grouping factor) or drop it as a group.
metadata <- read.table("../mapping.tsv", header = TRUE, sep = "\t", row.names = 1,
                       comment.char = "", na.strings = character(0))
# Drop the QIIME2 "#q2:types" row, if present.
if (nrow(metadata) > 0 && grepl("^#q2:types", rownames(metadata)[1])) {
    metadata <- metadata[-1, , drop = FALSE]
}
sample_data <- sample_data(metadata)

group_var        <- args[2]
block_var        <- if (length(args) >= 5 && nzchar(args[5])) args[5] else "none"
output_directory <- paste(args[1], "/diversity_vegan_TSS/", sep = "")
if (!dir.exists(output_directory)) dir.create(output_directory)

n_perm <- 999

# Returns the blocking factor (aligned on `metadata` rows) when it can be
# used for this comparison, NULL otherwise (with the reason printed).
# Usable = a block was given, differs from the comparison column, exists in
# the metadata, has > 1 level here, and block + group is not confounded
# (full-rank design, e.g. each group not restricted to a single animal).
get_block <- function(metadata, groups) {
    if (block_var == "none") return(NULL)
    if (block_var == group_var) {
        cat("Block '", block_var, "' is the comparison column: block analyses skipped.\n", sep = "")
        return(NULL)
    }
    if (!block_var %in% colnames(metadata)) {
        cat("Block column '", block_var, "' not found in mapping.tsv: block analyses skipped.\n", sep = "")
        return(NULL)
    }
    blk <- factor(as.character(metadata[[block_var]]))
    if (nlevels(blk) < 2) {
        cat("Block '", block_var, "' has a single level here: block analyses skipped.\n", sep = "")
        return(NULL)
    }
    mm <- model.matrix(~ blk + groups)
    if (qr(mm)$rank < ncol(mm)) {
        cat("Block '", block_var, "' is confounded with '", group_var,
            "' (rank-deficient design): block analyses skipped.\n", sep = "")
        return(NULL)
    }
    return(blk)
}

# Permutation scheme for a blocked test. Permutations restricted within
# blocks are used when they allow at least min_distinct_perm distinct
# permutations; otherwise (e.g. 3 animals x 2 samples = 8 permutations, so
# p >= 0.125) the test falls back to free permutations, the block staying
# in the model formula. Returns the `permutations` argument and a label.
min_distinct_perm <- 100
perm_scheme <- function(blk) {
    ctrl <- how(nperm = n_perm, blocks = blk)
    np   <- numPerms(length(blk), control = ctrl)
    if (np < min_distinct_perm) {
        return(list(perm = n_perm,
                    label = paste0("free (fallback: only ", np, " within-block permutations)")))
    }
    return(list(perm = ctrl, label = "within_block"))
}

# Tidy one adonis2 table (rows = terms) for writing.
tidy_adonis <- function(a, method, permutations, label) {
    res <- as.data.frame(a)
    res <- cbind(Term = rownames(res), res)
    res$Permutations  <- permutations
    res$Method        <- method
    res$Normalization <- label
    return(res)
}

analyze_diversity <- function(ps_object, label) {
    metadata <- data.frame(sample_data(ps_object))
    groups   <- factor(as.character(metadata[[group_var]]))

    # Beta diversity distances.
    unifrac   <- distance(ps_object, method = "unifrac",  weighted = FALSE)
    wunifrac  <- distance(ps_object, method = "unifrac",  weighted = TRUE)
    jaccard   <- distance(ps_object, method = "jaccard")
    bray_curtis <- distance(ps_object, method = "bray")

    # Recalculate UniFrac with the resolved tree.
    ps_new_tree        <- phyloseq(otu_table(ps_object), sample_data(ps_object), new_tre)
    unifrac_new_tree   <- distance(ps_new_tree, method = "unifrac", weighted = FALSE)
    wunifrac_new_tree  <- distance(ps_new_tree, method = "unifrac", weighted = TRUE)

    dists <- list(
        Unifrac                    = unifrac,
        Weighted_Unifrac           = wunifrac,
        Jaccard                    = jaccard,
        Bray_Curtis                = bray_curtis,
        Unifrac_New_Tree           = unifrac_new_tree,
        Weighted_Unifrac_New_Tree  = wunifrac_new_tree
    )

    # Align the metadata (and thus `groups`) on the sample order of the
    # distance matrices, instead of assuming both orders are the same.
    sample_order <- labels(unifrac)
    # The alpha-diversity part below binds estimate_richness() (phyloseq
    # sample order) to `metadata` by position: both orders must be equal.
    stopifnot(identical(sample_order, sample_names(ps_object)))
    metadata     <- metadata[sample_order, , drop = FALSE]
    groups       <- factor(as.character(metadata[[group_var]]))

    # Aitchison distance (CLR + Euclidean), raw counts only (see header).
    if (label == "non_normalized") {
        otu_m <- as(otu_table(ps_object), "matrix")
        if (!taxa_are_rows(ps_object)) otu_m <- t(otu_m)
        otu_m <- otu_m[, sample_order, drop = FALSE]
        dists$Aitchison <- vegdist(t(otu_m), method = "aitchison", pseudocount = 1)
    }

    # Safety check: every distance must share that sample order.
    for (dn in names(dists)) {
        if (!identical(labels(dists[[dn]]), sample_order)) {
            stop("Sample order of distance ", dn, " differs from the other distances.")
        }
    }

    # One-factor PERMANOVA (historical output, by = NULL -> "Model" row,
    # parsed by summary_results.pl).
    permanova_results <- lapply(dists, function(d) adonis2(d ~ groups, permutations = n_perm))

    pairwise_results <- lapply(dists, function(d) pairwise.adonis(d, factors = groups))

    pairwise_results_corrected <- lapply(pairwise_results, function(res) {
        res_df           <- as.data.frame(res)
        res_df$q_values_BH <- p.adjust(res_df$p.value, method = "BH")
        return(res_df)
    })

    if (label != "TSS_normalized") {
        alpha_diversity <- estimate_richness(ps_object,
            measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
        alpha_diversity <- cbind(metadata, alpha_diversity)

        kruskal_results <- list()
        for (measure in colnames(alpha_diversity)[-(1:ncol(metadata))]) {
            if (grepl("^se\\.", measure)) next
            kruskal_test <- kruskal.test(alpha_diversity[[measure]] ~ metadata[[group_var]], data = alpha_diversity)
            kruskal_results[[measure]] <- data.frame(
                Test      = "Kruskal-Wallis",
                Measure   = measure,
                Statistic = kruskal_test$statistic,
                p.value   = kruskal_test$p.value,
                Parameter = kruskal_test$parameter,
                Method    = kruskal_test$method,
                Normalization = label
            )
        }

        means_results <- list()
        sem_results   <- list()
        pairwise_wilcox_results <- list()

        for (measure in colnames(alpha_diversity)[-(1:ncol(metadata))]) {
            if (grepl("^se\\.", measure)) next
            clean_data <- na.omit(data.frame(
                value = alpha_diversity[[measure]],
                group = metadata[[group_var]]
            ))
            if (length(unique(metadata[[group_var]])) == length(unique(clean_data$group))) {
                mean_values <- aggregate(value ~ group, data = clean_data, mean)
                sem_values  <- aggregate(value ~ group, data = clean_data,
                                         function(x) sd(x) / sqrt(length(x)))
                mean_values$Measure    <- measure
                sem_values$Measure     <- measure
                colnames(sem_values)[2] <- "sem"
                means_results[[measure]] <- mean_values
                sem_results[[measure]]   <- sem_values

                pairwise_wilcox <- pairwise.wilcox.test(
                    clean_data$value, clean_data$group, p.adjust.method = "fdr")
                pairwise_wilcox_results[[measure]] <- data.frame(
                    Test    = "Pairwise Wilcoxon",
                    Measure = measure,
                    p.values = pairwise_wilcox$p.value,
                    Normalization = label
                )
            }
        }

        all_kruskal         <- do.call(rbind, kruskal_results)
        all_pairwise_wilcox <- do.call(rbind, pairwise_wilcox_results)
        all_means           <- do.call(rbind, means_results)
        all_sem             <- do.call(rbind, sem_results)

        write.table(all_kruskal,
            file.path(output_directory, paste0("kruskal_results_",         label, ".tsv")),
            sep = "\t", col.names = NA, quote = FALSE)
        write.table(all_pairwise_wilcox,
            file.path(output_directory, paste0("pairwise_wilcox_results_", label, ".tsv")),
            sep = "\t", col.names = NA, quote = FALSE)
        write.table(all_means,
            file.path(output_directory, paste0("mean_values_",             label, ".tsv")),
            sep = "\t", col.names = NA, quote = FALSE)
        write.table(all_sem,
            file.path(output_directory, paste0("sem_values_",              label, ".tsv")),
            sep = "\t", col.names = NA, quote = FALSE)
    }

    write.table(do.call(rbind, lapply(names(permanova_results), function(x) {
        res           <- as.data.frame(permanova_results[[x]])
        res$Test      <- "PERMANOVA"
        res$Method    <- x
        res$Normalization <- label
        return(res)
    })), file.path(output_directory, paste0("permanova_results_", label, ".tsv")),
    sep = "\t", col.names = NA, quote = FALSE)

    write.table(do.call(rbind, lapply(names(pairwise_results_corrected), function(x) {
        res               <- pairwise_results_corrected[[x]]
        res$Test          <- "Pairwise PERMANOVA"
        res$Method        <- x
        res$Normalization <- label
        return(res)
    })), file.path(output_directory, paste0("pairwise_results_", label, ".tsv")),
    sep = "\t", col.names = NA, quote = FALSE)

    # ---- Dispersion test (independent of the block) ----------------------
    # A significant PERMANOVA can come from a difference in dispersion
    # rather than in location, especially with unbalanced groups.
    if (nlevels(groups) > 1) {
        betadisper_res <- do.call(rbind, lapply(names(dists), function(dn) {
            bd <- betadisper(dists[[dn]], groups)
            pt <- permutest(bd, permutations = n_perm)
            tab <- as.data.frame(pt$tab)
            data.frame(Method = dn, Normalization = label,
                       Df = tab[1, "Df"], F = tab[1, "F"],
                       p.value = tab[1, "Pr(>F)"],
                       # Actual number of permutations (fewer than n_perm
                       # when all of them can be enumerated).
                       Permutations = if (NROW(pt$perm) > 0) NROW(pt$perm) else n_perm,
                       mean_dist_to_centroid = paste(names(tapply(bd$distances, groups, mean)),
                                                     signif(tapply(bd$distances, groups, mean), 4),
                                                     sep = "=", collapse = ";"))
        }))
        write.table(betadisper_res,
            file.path(output_directory, paste0("betadisper_", label, ".tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)
    }

    # ---- Analyses with the blocking factor -------------------------------
    blk <- get_block(metadata, groups)
    if (!is.null(blk)) {
        df <- data.frame(blk = blk, grp = groups)

        # PERMANOVA group + block, marginal tests. Free permutations, then
        # permutations restricted within blocks (the p-value of the block
        # term itself is not meaningful in the blocked version: block labels
        # are never exchanged, so it is not reported).
        ps_all <- perm_scheme(df$blk)
        cat("Blocked PERMANOVA (", label, "): ", ps_all$label, "\n", sep = "")
        perm_block <- do.call(rbind, lapply(names(dists), function(dn) {
            d <- dists[[dn]]
            a_free  <- adonis2(d ~ blk + grp, data = df, by = "margin",
                               permutations = n_perm)
            a_block <- adonis2(d ~ blk + grp, data = df, by = "margin",
                               permutations = ps_all$perm)
            t_free  <- tidy_adonis(a_free,  dn, "free",       label)
            t_block <- tidy_adonis(a_block, dn, ps_all$label, label)
            if (ps_all$label == "within_block") t_block[t_block$Term == "blk", "Pr(>F)"] <- NA
            rbind(t_free, t_block)
        }))
        perm_block$Term[perm_block$Term == "blk"] <- block_var
        perm_block$Term[perm_block$Term == "grp"] <- group_var
        write.table(perm_block,
            file.path(output_directory, paste0("permanova_block_", label, ".tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

        # Pairwise blocked PERMANOVA (only meaningful with > 2 groups).
        if (nlevels(groups) > 2) {
            pairs <- combn(levels(groups), 2, simplify = FALSE)
            pw <- do.call(rbind, lapply(names(dists), function(dn) {
                do.call(rbind, lapply(pairs, function(pr) {
                    keep  <- groups %in% pr
                    d_sub <- as.dist(as.matrix(dists[[dn]])[keep, keep])
                    df_s  <- data.frame(blk = factor(as.character(blk[keep])),
                                        grp = factor(as.character(groups[keep])))
                    mm <- model.matrix(~ blk + grp, data = df_s)
                    if (nlevels(df_s$blk) < 2 || qr(mm)$rank < ncol(mm)) {
                        return(data.frame(Method = dn, Normalization = label,
                                          group1 = pr[1], group2 = pr[2],
                                          R2 = NA, F = NA, p.value = NA,
                                          note = "block not estimable in this pair"))
                    }
                    ps_pair <- perm_scheme(df_s$blk)
                    a <- adonis2(d_sub ~ blk + grp, data = df_s, by = "margin",
                                 permutations = ps_pair$perm)
                    data.frame(Method = dn, Normalization = label,
                               group1 = pr[1], group2 = pr[2],
                               R2 = a["grp", "R2"], F = a["grp", "F"],
                               p.value = a["grp", "Pr(>F)"], note = ps_pair$label)
                }))
            }))
            # BH within each distance.
            pw$q.value.BH <- ave(pw$p.value, pw$Method,
                                 FUN = function(p) p.adjust(p, method = "BH"))
            write.table(pw,
                file.path(output_directory, paste0("pairwise_block_", label, ".tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
        }

        # Alpha diversity ~ block + group (not on TSS: richness estimators
        # need integer counts).
        if (label != "TSS_normalized") {
            alpha_div <- estimate_richness(ps_object,
                measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson"))
            # estimate_richness() may mangle sample names (make.names), so
            # align by the phyloseq sample order, not by row names.
            alpha_div <- alpha_div[match(sample_order, sample_names(ps_object)), , drop = FALSE]
            alpha_lm <- do.call(rbind, lapply(colnames(alpha_div), function(ms) {
                if (grepl("^se\\.", ms)) return(NULL)
                dfa <- data.frame(value = alpha_div[[ms]], blk = blk, grp = groups)
                dfa <- droplevels(dfa[is.finite(dfa$value), , drop = FALSE])
                d1  <- tryCatch(drop1(lm(value ~ blk + grp, data = dfa), test = "F"),
                                error = function(e) NULL)
                if (is.null(d1)) return(NULL)
                data.frame(Measure = ms, Normalization = label,
                           Term = c(block_var, group_var),
                           Df = d1[c("blk", "grp"), "Df"],
                           F = d1[c("blk", "grp"), "F value"],
                           p.value = d1[c("blk", "grp"), "Pr(>F)"])
            }))
            if (!is.null(alpha_lm)) {
                alpha_lm$q.value.BH <- ave(alpha_lm$p.value, alpha_lm$Term,
                                           FUN = function(p) p.adjust(p, method = "BH"))
                write.table(alpha_lm,
                    file.path(output_directory, paste0("alpha_lm_", label, ".tsv")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
            }
        }
    }
}

# Run analyses under three normalisations.
ps_original <- phyloseq(otu_table, sample_data, phy_tree)
analyze_diversity(ps_original, "non_normalized")

ps_rarefied <- rarefy_even_depth(ps_original, sample.size = args[3], rngseed = 10, replace = FALSE)
analyze_diversity(ps_rarefied, "rarefied")

# TSS normalisation (total-sum scaling to 1e6).
otu_matrix <- as(otu_table(ps_original), "matrix")
tss_norm   <- sweep(otu_matrix, 2, colSums(otu_matrix), "/") * 1e6
otu_table(ps_original) <- otu_table(tss_norm, taxa_are_rows = TRUE)
analyze_diversity(ps_original, "TSS_normalized")

# Write completion flag.
outfile <- paste(args[1], "/done_diversity_vegan_TSS", sep = "")
write.table("done", outfile, sep = "\t", row.names = FALSE)