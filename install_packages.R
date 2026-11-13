# =============================================================================
# install_packages.R - install the R packages used in this study
# -----------------------------------------------------------------------------
# Tested with R 4.4.3. Run once from the repository root:
#     source("install_packages.R")
# The script skips packages that are already installed.
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran <- c(
  # general
  "ggplot2", "dplyr", "tidyr", "stringr", "cowplot", "patchwork",
  "RColorBrewer", "viridis", "scales", "openxlsx",
  # bulk transcriptomics
  "pheatmap", "ggpubr", "ggrepel", "VennDiagram", "ggvenn", "circlize",
  "msigdbr", "igraph",
  # machine learning
  "caret", "randomForest", "kernlab", "xgboost", "pROC",
  "DALEX", "shapviz", "ggbeeswarm", "e1071", "nnet", "gbm", "ranger",
  # single cell
  "Seurat", "harmony", "scTenifoldKnk", "Matrix", "future",
  "scRNAtoolVis"   # GitHub only, handled below
)
cran <- setdiff(cran, "scRNAtoolVis")

bioc <- c(
  "limma", "sva", "GEOquery",
  "org.Hs.eg.db", "clusterProfiler", "enrichplot", "DOSE",
  "ReactomePA", "GOSemSim", "pathview", "ComplexHeatmap",
  "SingleR", "AUCell", "GSVA", "monocle3", "celldex"
)

github <- c(
  "satijalab/seurat-wrappers" = "SeuratWrappers",
  "junjunlab/scRNAtoolVis"    = "scRNAtoolVis",
  "jinworks/CellChat"         = "CellChat"
)

miss_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss_cran)) install.packages(miss_cran)

miss_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss_bioc)) BiocManager::install(miss_bioc, update = FALSE, ask = FALSE)

if (length(github)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  for (repo in names(github)) {
    pkg <- github[[repo]]
    if (!requireNamespace(pkg, quietly = TRUE)) remotes::install_github(repo)
  }
}

# --- record the exact versions that were used --------------------------------
installed <- c(cran, bioc, unname(github))
vers <- vapply(installed, function(p) {
  tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
}, character(1))
vers <- vers[!is.na(vers)]
cat("R version:", R.version.string, "\n")
write.table(data.frame(package = names(vers), version = unname(vers)),
            file = "package_versions.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote package_versions.txt\n")
