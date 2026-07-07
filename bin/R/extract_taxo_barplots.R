# extract_taxo_barplots.R
# Exports taxonomy barplots (absolute and relative abundances at all levels)
# to an Excel workbook, with one sheet per level and per type.
#
# Author: Jeremy Tournayre
#
# Arguments:
#   args[1] - feature table (.qza)
#   args[2] - taxonomy (.qza)
#   args[3] - output Excel file (.xlsx)

rm(list = ls())

# Load R library path from environment (set in config.sh).
user_lib <- Sys.getenv("PIPELINE_R_LIBS")
if (nzchar(user_lib)) .libPaths(user_lib)

library(purrr)
library(qiime2R)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)

args <- commandArgs(trailingOnly = TRUE)

# Remove existing output file if present.
unlink(args[3])

output_dir    <- dirname(args[3])
plots_tmp_dir <- file.path(output_dir, "plots_tmp")

if (dir.exists(plots_tmp_dir)) unlink(plots_tmp_dir, recursive = TRUE)
dir.create(plots_tmp_dir, showWarnings = FALSE, recursive = TRUE)

# Load feature table and taxonomy.
feature_table <- read_qza(args[1])$data
taxonomy      <- read_qza(args[2])$data

feature_df <- as.data.frame(as.matrix(feature_table)) %>%
    tibble::rownames_to_column("FeatureID")

# Split taxonomy string into the 7 standard levels.
taxonomy_df <- taxonomy %>%
    select(FeatureID = Feature.ID, Taxon) %>%
    separate(Taxon,
             into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
             sep = ";", fill = "right") %>%
    mutate(across(Kingdom:Species, ~ trimws(.)))

taxonomy_df <- taxonomy_df %>%
    mutate(Full_Taxon = do.call(paste, c(across(Kingdom:Species), sep = "; ")))

merged_df <- left_join(feature_df, taxonomy_df, by = "FeatureID")

feature_long <- merged_df %>%
    pivot_longer(cols = -c(FeatureID, Full_Taxon, Kingdom:Species),
                 names_to = "SampleID", values_to = "Abundance")

total_reads <- feature_long %>%
    group_by(SampleID) %>%
    summarise(Total_Reads = sum(Abundance), .groups = "drop")

levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
relative_abundances <- list()
absolute_abundances <- list()

for (level in levels) {
    levels_to_concat <- levels[1:which(levels == level)]

    feature_long <- feature_long %>%
        mutate(
            Full_Taxon_Level = pmap_chr(
                select(., all_of(levels_to_concat)),
                ~ {
                    vals <- c(...)
                    vals[is.na(vals) | vals == "" | vals == "Unassigned"] <- "Unassigned"
                    paste(vals, collapse = "; ")
                }
            )
        ) %>%
        mutate(Taxon = sub(".*;\\s*", "", Full_Taxon_Level))

    rel_df <- feature_long %>%
        group_by(SampleID) %>%
        mutate(RelAbundance = Abundance / sum(Abundance)) %>%
        group_by(SampleID, Full_Taxon_Level, Taxon) %>%
        summarise(RelAbundance = sum(RelAbundance), .groups = "drop") %>%
        mutate(Level = level)

    relative_abundances[[level]] <- rel_df

    abs_df <- feature_long %>%
        group_by(SampleID, Full_Taxon_Level, Taxon) %>%
        summarise(Abundance = sum(Abundance), .groups = "drop") %>%
        mutate(Level = level)

    absolute_abundances[[level]] <- abs_df
}

wb <- createWorkbook()

addWorksheet(wb, sheetName = "Total_Reads")
writeData(wb, sheet = "Total_Reads", x = total_reads)

for (level in levels) {
    df_abs <- absolute_abundances[[level]] %>%
        pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0) %>%
        relocate(Taxon, .before = everything()) %>%
        select(-Level)

    addWorksheet(wb, sheetName = paste0("AbsAb_", level))
    writeData(wb, sheet = paste0("AbsAb_", level), x = df_abs)
}

for (level in levels) {
    df_rel <- relative_abundances[[level]] %>%
        pivot_wider(names_from = SampleID, values_from = RelAbundance, values_fill = 0) %>%
        relocate(Taxon, .before = everything()) %>%
        select(-Level)

    addWorksheet(wb, sheetName = paste0("RelAb_", level))
    writeData(wb, sheet = paste0("RelAb_", level), x = df_rel)
}

# Plot sheet: top-10 taxa at each level (or all if <= 10).
for (level in levels) {
    sheet_name <- paste0("plot_RelAb_", level)
    addWorksheet(wb, sheetName = sheet_name)

    n_samples <- length(unique(relative_abundances[[level]]$SampleID))
    n_taxa    <- n_distinct(relative_abundances[[level]]$Full_Taxon_Level)

    if (n_taxa > 10) {
        top_taxa <- relative_abundances[[level]] %>%
            group_by(Taxon) %>%
            summarise(TotalAbundance = sum(RelAbundance)) %>%
            arrange(desc(TotalAbundance)) %>%
            slice_head(n = 10) %>%
            pull(Taxon)

        df_plot    <- relative_abundances[[level]] %>%
            mutate(Taxon = ifelse(Taxon %in% top_taxa, Taxon, "Other"))
        plot_title <- paste("Relative Abundance -", level, "(Top 10 + Other)")
    } else {
        df_plot    <- relative_abundances[[level]]
        plot_title <- paste("Relative Abundance -", level)
    }

    p <- ggplot(df_plot, aes(x = SampleID, y = RelAbundance, fill = Taxon)) +
        geom_bar(stat = "identity") +
        theme_minimal(base_size = 14) +
        ylab("Relative Abundance") +
        ggtitle(plot_title) +
        theme(
            plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
            axis.text.x   = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),
            axis.text.y   = element_text(size = 20),
            axis.title    = element_text(size = 20),
            legend.position = "bottom",
            legend.text   = element_text(size = 20),
            legend.title  = element_text(size = 20),
            legend.key.size = unit(1.25, "cm")
        ) +
        guides(fill = guide_legend(nrow = 2))

    plot_width  <- max(5, n_samples * 0.25)
    plot_height <- 10
    plot_file   <- file.path(plots_tmp_dir, paste0("plot_", level, ".png"))
    ggsave(plot_file, plot = p, width = plot_width, height = plot_height, dpi = 300)

    insertImage(wb, sheet = sheet_name, file = plot_file,
                startRow = 1, startCol = 1,
                width = plot_width * 1, height = plot_height * 1, units = "cm")
}

saveWorkbook(wb, args[3], overwrite = TRUE)
if (dir.exists(plots_tmp_dir)) unlink(plots_tmp_dir, recursive = TRUE)