# get_lowest_depth.R
# Reads a feature table in biom format and prints the minimum total read
# count across all samples. The result is used by qiime2.sh to set the
# rarefaction depth.
#
# Author: Jeremy Tournayre
#
# Arguments:
#   args[1] - base directory containing the unzipped table folder
#   args[2] - QIIME2 artifact UUID (subfolder name inside the table/ directory)

# Load R library path from environment (set in config.sh).
user_lib <- Sys.getenv("PIPELINE_R_LIBS")
if (nzchar(user_lib)) .libPaths(user_lib)

library(qiime2R)

args <- commandArgs(trailingOnly = TRUE)

biom_file  <- paste(args[1], '/table/', args[2], '/data/feature-table.biom', sep = "")
table      <- read_q2biom(biom_file)

get_lowest <- min(apply(table, 2, sum))
print(get_lowest)