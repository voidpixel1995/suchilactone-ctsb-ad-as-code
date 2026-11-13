# =============================================================================
# 00_download_data.R - fetch the public data used in this study
# -----------------------------------------------------------------------------
# The processed expression matrices are not stored in this repository because of
# their size. This script downloads the raw series from GEO/ArrayExpress and
# writes one folder per cohort under data/download/, which can then be copied to
# the working directories described in README.md.
#
# Requires: BiocManager::install("GEOquery")
# =============================================================================

local({
  d <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(d, "00_config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "00_config.R")))
    stop("00_config.R not found - run this script from inside the repository.")
  assign("PROJECT_ROOT", d, envir = .GlobalEnv)
  source(file.path(d, "00_config.R"))
  invisible(NULL)
})

suppressPackageStartupMessages({
  if (!requireNamespace("GEOquery", quietly = TRUE))
    stop("Please install GEOquery: BiocManager::install('GEOquery')")
  library(GEOquery)
})

outdir <- file.path(DATA_DIR, "download")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- bulk transcriptomes (five AD cohorts, three AS cohorts, 553 samples) ----
bulk <- c("GSE138260", "GSE29378", "GSE36980", "GSE37263", "GSE5281",      # AD
          "GSE100927", "GSE28829", "GSE43292")                            # AS

for (g in bulk) {
  dest <- file.path(outdir, g)
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  message("== bulk: ", g)
  try({
    eset <- getGEO(g, GSEMatrix = TRUE, destdir = dest, getGPL = TRUE)
    exprs <- Biobase::exprs(eset[[1]])
    fdata <- Biobase::fData(eset[[1]])
    pdata <- Biobase::pData(eset[[1]])

    # collapse probes to gene symbols, keeping the probe with the largest mean
    sym_col <- intersect(c("Gene symbol", "GENE_SYMBOL", "Symbol",
                           "Gene Symbol", "gene_assignment"), colnames(fdata))
    if (length(sym_col) > 0) {
      sym <- toupper(sub("///.*$", "", fdata[[sym_col[1]]]))
      keep <- !is.na(sym) & nzchar(sym)
      exprs <- exprs[keep, , drop = FALSE]
      sym <- sym[keep]
      ord <- order(rowMeans(exprs), decreasing = TRUE)
      exprs <- exprs[ord, , drop = FALSE]
      sym <- sym[ord]
      exprs <- exprs[!duplicated(sym), , drop = FALSE]
      rownames(exprs) <- sym[!duplicated(sym)]
    }
    out <- data.frame(Gene = rownames(exprs), exprs, check.names = FALSE)
    write.table(out, file.path(dest, "geneMatrix.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.csv(pdata, file.path(dest, "phenotype.csv"), row.names = FALSE)

    # the group definition used in this study is kept in data/metadata/
    meta <- file.path(DATA_DIR, "metadata",
                      ifelse(g %in% c("GSE100927", "GSE28829", "GSE43292"), "AS", "AD"), g)
    if (dir.exists(meta)) {
      file.copy(file.path(meta, "control_samples.txt"), file.path(dest, "s1.txt"), overwrite = TRUE)
      file.copy(file.path(meta, "case_samples.txt"),    file.path(dest, "s2.txt"), overwrite = TRUE)
      message("   group files copied from data/metadata/")
    }
  }, silent = FALSE)
}

# --- single-cell transcriptomes ---------------------------------------------
# GSE157827 (AD, 21 samples) and GSE159677 (AS, 6 samples) are distributed as
# tar archives of 10x Genomics matrices; download them from the GEO ftp page:
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677
# and untar into data/download/GSE157827/ and data/download/GSE159677/.
message("\nSingle-cell series must be downloaded manually - see the comment above.")

message("\nDone. Data are in: ", outdir)
