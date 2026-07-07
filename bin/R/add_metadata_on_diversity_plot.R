# add_metadata_on_diversity_plot.R
# Reads exported PCoA ordination coordinates and overlays numeric metadata
# variables using envfit (vegan). Saves a PNG coloured by a categorical group.
#
# Metadata comes from two files: mapping.tsv (used by QIIME2 itself) and
# mapping_extra.tsv (additional numeric variables, e.g. physico-chemical
# measurements, kept separate because sparse/NA-heavy columns can make
# QIIME2's own metadata-driven commands fail). The two are merged here by
# sample ID before projecting variables onto the ordination.
#
# Arguments:
#   args[1] - ordination.txt exported by qiime tools export
#   args[2] - mapping.tsv metadata file (used by QIIME2)
#   args[3] - metadata column to use for colouring (categorical group)
#   args[4] - output file base name (without extension)
#   args[5] - mapping_extra.tsv (additional numeric metadata; optional)

rm(list = ls())

# Author: Jeremy Tournayre

# Fixed seed for reproducibility of envfit's random permutations below.
set.seed(42)

# Load R library path from environment (set in config.sh).
user_lib <- Sys.getenv("PIPELINE_R_LIBS")
if (nzchar(user_lib)) .libPaths(user_lib)

library(vegan)
library(Cairo)

args <- commandArgs(trailingOnly = TRUE)

lines <- readLines(args[1])

# Extract variance explained by each axis.
prop_line_idx <- grep("^Proportion explained", lines)
prop_line     <- lines[prop_line_idx + 1]
prop_expl     <- as.numeric(unlist(strsplit(prop_line, "\t")))

# Extract site scores.
start_idx <- grep("^Site", lines) + 1
end_idx   <- grep("^(Biplot|Site constraints|$)", lines)
end_idx   <- end_idx[end_idx > start_idx][1] - 1
site_lines <- lines[start_idx:end_idx]

ordination          <- read.table(text = site_lines, header = FALSE, sep = "\t", row.names = 1)
colnames(ordination) <- c("PCoA1", "PCoA2")

# Load metadata used by QIIME2.
metadata <- read.table(args[2], header = TRUE, sep = "\t", row.names = 1, comment.char = "",
                       na.strings = c("NA", "N/A", ""))

# Merge in the additional numeric metadata file, if provided, matching rows
# by sample ID. Columns from mapping_extra.tsv are added alongside those
# from mapping.tsv; sample IDs not present in mapping_extra.tsv simply get
# NA for those columns.
if (length(args) >= 5 && nzchar(args[5]) && file.exists(args[5])) {
    metadata_extra <- read.table(args[5], header = TRUE, sep = "\t", row.names = 1, comment.char = "",
                                 na.strings = c("NA", "N/A", ""))
    metadata <- cbind(metadata, metadata_extra[rownames(metadata), , drop = FALSE])
}

common_samples <- intersect(rownames(ordination), rownames(metadata))
ordination     <- ordination[common_samples, ]
metadata       <- metadata[common_samples, ]

# Convert numeric-looking character columns to numeric.
metadata[] <- lapply(metadata, function(x) {
    if (is.character(x) && all(grepl("^\\s*-?\\d*\\.?\\d+\\s*$", x[!is.na(x)]))) {
        as.numeric(x)
    } else {
        x
    }
})

# Project numeric variables onto the ordination with envfit.
num_vars <- metadata[sapply(metadata, is.numeric)]
if (ncol(num_vars) == 0) {
    message("No valid numeric variable to project with envfit. Ending script.")
    quit(save = "no", status = 0)
}
num_vars  <- num_vars[sapply(num_vars, function(x) length(unique(x)) > 1)]
fit_all   <- envfit(ordination, num_vars, permutations = 999)

df_scores      <- ordination
groups         <- factor(metadata[rownames(df_scores), args[3]])
palette_colors <- rainbow(length(levels(groups)))

xlab <- paste0("PCoA1 (", round(prop_expl[1] * 100, 1), "%)")
ylab <- paste0("PCoA2 (", round(prop_expl[2] * 100, 1), "%)")

CairoPNG(paste0(args[4], ".png"), width = 700, height = 700, res = 120)

max_title_len <- 90
titre_plot <- if (nchar(args[4]) > max_title_len) {
    paste0("...", substr(args[4], nchar(args[4]) - max_title_len + 1, nchar(args[4])))
} else {
    args[4]
}

plot(df_scores$PCoA1, df_scores$PCoA2,
     col  = palette_colors[groups], pch = 19, cex = 1.5,
     xlab = xlab, ylab = ylab,
     main = titre_plot, cex.main = 0.7)
legend("topright", legend = levels(groups), col = palette_colors, pch = 19, bty = "n")
plot(fit_all, col = "blue", add = TRUE)

dev.off()

print(fit_all)