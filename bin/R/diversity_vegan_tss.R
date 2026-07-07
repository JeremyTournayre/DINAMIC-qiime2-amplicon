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
output_directory <- paste(args[1], "/diversity_vegan_TSS/", sep = "")
if (!dir.exists(output_directory)) dir.create(output_directory)

analyze_diversity <- function(ps_object, label) {
    metadata <- data.frame(sample_data(ps_object))
    groups   <- metadata[[group_var]]

    # Beta diversity distances.
    unifrac   <- distance(ps_object, method = "unifrac",  weighted = FALSE)
    wunifrac  <- distance(ps_object, method = "unifrac",  weighted = TRUE)
    jaccard   <- distance(ps_object, method = "jaccard")
    bray_curtis <- distance(ps_object, method = "bray")

    # Recalculate UniFrac with the resolved tree.
    ps_new_tree        <- phyloseq(otu_table(ps_object), sample_data(ps_object), new_tre)
    unifrac_new_tree   <- distance(ps_new_tree, method = "unifrac", weighted = FALSE)
    wunifrac_new_tree  <- distance(ps_new_tree, method = "unifrac", weighted = TRUE)

    permanova_results <- list(
        Unifrac                    = adonis2(unifrac          ~ groups),
        Weighted_Unifrac           = adonis2(wunifrac         ~ groups),
        Jaccard                    = adonis2(jaccard          ~ groups),
        Bray_Curtis                = adonis2(bray_curtis      ~ groups),
        Unifrac_New_Tree           = adonis2(unifrac_new_tree  ~ groups),
        Weighted_Unifrac_New_Tree  = adonis2(wunifrac_new_tree ~ groups)
    )

    pairwise_results <- list(
        Unifrac                    = pairwise.adonis(unifrac,          factors = groups),
        Weighted_Unifrac           = pairwise.adonis(wunifrac,         factors = groups),
        Jaccard                    = pairwise.adonis(jaccard,          factors = groups),
        Bray_Curtis                = pairwise.adonis(bray_curtis,      factors = groups),
        Unifrac_New_Tree           = pairwise.adonis(unifrac_new_tree,  factors = groups),
        Weighted_Unifrac_New_Tree  = pairwise.adonis(wunifrac_new_tree, factors = groups)
    )

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