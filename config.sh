#!/usr/bin/env bash
# config.sh - central configuration for the holoruminant QIIME2 pipeline.
#
# Author: Jeremy Tournayre
#
# Every entry point sources this file first:  source "$(dirname "$0")/config.sh"
# This is the ONLY place machine-specific paths should live.

# --- Pipeline root (directory that contains this file) -----------------------
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR

# --- Conda -------------------------------------------------------------------
# Conda profile so that `conda activate` works in a NON-interactive shell.
CONDA_PROFILE="${CONDA_PROFILE:-/usr/local/genome/miniforge3/etc/profile.d/conda.sh}"
export CONDA_PROFILE

# QIIME2 environment used by every script in this repository (main analysis
# pipeline via launch.pl -> qiime2.sh, and the database-building scripts
# under db/). Before rerunning an existing project under a different QIIME2
# version, make sure this value reflects the version you intend to use.
CONDA_ENV_QIIME2="${CONDA_ENV_QIIME2:-qiime2-amplicon-2024.5}"
export CONDA_ENV_QIIME2

# --- R -------------------------------------------------------------------
# R library path used by helper scripts after `conda deactivate` (system R).
PIPELINE_R_LIBS="/home/jtournayre/R/x86_64-pc-linux-gnu-library/4.5"
export PIPELINE_R_LIBS
export R_LIBS_USER="$PIPELINE_R_LIBS"

# --- Vendored / patched dependencies ----------------------------------------
# Local patched copy of ANCOMBC sourced by run_ancombc.R.
ANCOMBC_DIR="$PIPELINE_DIR/bin/tools/ANCOMBC-2025_02_20"
export ANCOMBC_DIR

# libstdc++ shim injected via LD_LIBRARY_PATH. Optional: set to "" to disable.
LIBSTDCXX_LIB="/home/jtournayre/work/.conda/envs/myenv_libstdcxx/lib"
export LIBSTDCXX_LIB

# --- Databases ---------------------------------------------------------------
# Classifiers live under db/Data/<Database>/<QIIME2_env>/<Primer>/classifier.qza
# — segregated by QIIME2 environment because .qza artifact formats can differ
# between QIIME2 releases (see db/formatters/README.md).
DB_DIR="$PIPELINE_DIR/db"
export DB_DIR

# --- Helper scripts -----------------------------------------------------
RSCRIPTS_DIR="$PIPELINE_DIR/bin/R"
PLSCRIPTS_DIR="$PIPELINE_DIR/bin/perl"
export RSCRIPTS_DIR PLSCRIPTS_DIR

# --- Locale --------------------------------------------------------------
export LANGUAGE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# --- Activation helper ----------------------------------------------------
activate_qiime2() {
    # shellcheck disable=SC1090
    source "$CONDA_PROFILE"
    conda activate "$CONDA_ENV_QIIME2"
    [[ -n "$LIBSTDCXX_LIB" ]] && export LD_LIBRARY_PATH="$LIBSTDCXX_LIB:$LD_LIBRARY_PATH"
}