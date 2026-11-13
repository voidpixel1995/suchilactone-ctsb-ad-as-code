# =============================================================================
# Step 3.8 - Publication-quality expression plots for the selected genes
# -----------------------------------------------------------------------------
# Part of the reproducibility package for:
#   Suchilactone engages cathepsin B to attenuate NF-kB/NLRP3-driven
#   cerebrovascular endothelial injury
# Run this script from anywhere inside the repository; the block below locates
# the repository root and loads the shared settings in 00_config.R.
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

#install.packages("clustree")
#install.packages("Seurat")
#install.packages("harmony")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("monocle")

#install.packages("devtools")
#devtools::install_github('immunogenomics/presto')

#install.packages("assertthat")
#install.packages("SCpubr")


#引用包
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(clustree)
library(harmony)
library(assertthat)
library(SCpubr)

#设置工作目录
workDir="E:\\1SX\\202604网络药理+心脑血管\\AD\\28.geneView"
setwd(file.path(WORK_DIR, "singlecell"))   # set by 00_config.R
load("Seurat.Rdata")

#读取交集基因的文件
geneRT=read.table("interGenes.List.txt", header=F, sep="\t", check.names=F)
interGenes=as.vector(geneRT[,1])
interGenes=interGenes[interGenes %in% rownames(pbmc)]

#绘制散点图
pdf(file="Scatter.pdf", width=8.5, height=5)
FeaturePlot(object = pbmc, features = interGenes, cols = c("green", "red"), ncol=3)
dev.off()

#绘制小提琴图
pdf(file="Violin.pdf", width=9.5, height=7)
VlnPlot(object = pbmc, features = interGenes, ncol=3)
dev.off()

#绘制气泡图
pdf(file="Bubble.pdf", width=8, height=5)
DotPlot(object = pbmc, features = interGenes) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
dev.off()
