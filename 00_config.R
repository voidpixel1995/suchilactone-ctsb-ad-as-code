# =============================================================================
# 00_config.R - shared settings for every script in this repository
# -----------------------------------------------------------------------------
# This file is loaded automatically by the bootstrap block at the top of each
# analysis script. You normally do not need to edit anything here.
#
# Directory convention
#   data/    small tables that are versioned in git (metadata, reference lists)
#   work/    working directories of the analysis steps; NOT versioned, created
#            automatically, and safe to delete at any time
#   results/ curated figures and tables that you keep for the record; NOT
#            versioned unless you remove it from .gitignore
# =============================================================================

if (!exists("PROJECT_ROOT")) {
  PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

DATA_DIR    <- file.path(PROJECT_ROOT, "data")
WORK_DIR    <- file.path(PROJECT_ROOT, "work")
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")

invisible(lapply(c(WORK_DIR, RESULTS_DIR),
                 dir.create, showWarnings = FALSE, recursive = TRUE))

# --- reproducibility ---------------------------------------------------------
SEED <- 12345          # single-cell analyses
SEED_ML <- 123          # machine-learning models

# --- analysis parameters reported in the manuscript --------------------------
LOGFC_DEG     <- 0.585   # |log2FC| cut-off for bulk DEGs (limma), equivalent to FC 1.5
PADJ_DEG      <- 0.05    # Benjamini-Hochberg adjusted p-value cut-off
LOGFC_PERTURB <- 1.0     # |log2FC| cut-off for the in silico perturbation DEGs
PCA_DIMS      <- 50      # number of principal components computed in scRNA-seq
DIMS_USE      <- 20      # components retained for UMAP / SNN graph
RES_CLUSTER   <- 0.5     # FindClusters resolution

options(stringsAsFactors = FALSE)
set.seed(SEED)

message("PROJECT_ROOT: ", PROJECT_ROOT)
