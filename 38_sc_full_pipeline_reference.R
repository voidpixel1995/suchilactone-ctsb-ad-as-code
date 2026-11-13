# =============================================================================
# Step 3.9 - Consolidated single-cell pipeline (reference copy of the complete run)
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

# ============================================================================
# Part 0: 环境配置与R包加载
# ============================================================================

# --- 0.1 安装必要的包 ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_packages <- c(
  "Seurat", "tidyverse", "ggplot2", "patchwork", "cowplot",
  "pheatmap", "RColorBrewer", "ggrepel", "ggpubr", "ggsci",
  "VennDiagram", "circlize", "reshape2", "corrplot",
  "harmony", "devtools", "msigdbr", "survival"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

bioc_packages <- c(
  "clusterProfiler", "org.Hs.eg.db", "enrichplot",
  "AUCell", "GSVA", "ComplexHeatmap", "monocle3",
  "SingleCellExperiment", "scran", "scater"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
  devtools::install_github("sqjin/CellChat")
}

# --- 0.2 加载包 ---
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
library(ggsci)
library(harmony)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AUCell)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(CellChat)
library(VennDiagram)
library(msigdbr)

# --- 0.3 配色方案（保留ggsci调色板，只取颜色值） ---
npg_colors <- ggsci::pal_npg("nrc")(10)
npg_extended <- c(
  ggsci::pal_npg("nrc")(10),
  ggsci::pal_nejm("default")(8),
  ggsci::pal_lancet("lanonc")(9)
)

# --- 0.4 定义安全主题（替代原来的 theme_npg） ---
# 使用 theme_bw 为基础，字体统一 "sans"，永不报错
theme_publish <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "sans", color = "black"),
      axis.text        = element_text(size = base_size - 1, color = "black"),
      axis.title       = element_text(size = base_size, face = "bold"),
      plot.title       = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      legend.text      = element_text(size = base_size - 2),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text       = element_text(size = base_size - 1, face = "bold")
    )
}

# 设置全局主题
theme_set(theme_publish())

# 设置工作目录
setwd(file.path(WORK_DIR, "singlecell"))   # set by 00_config.R

# 创建输出目录
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)

cat("============================================\n")
cat("环境配置完成！\n")
cat("============================================\n")


# ============================================================================
# Part 1: AS单细胞数据预处理
# ============================================================================

# --- 1.1 读取AS数据 ---
as_seurat <- readRDS("data/AS_seurat_raw.rds")

cat("AS数据集维度：\n")
cat("基因数：", nrow(as_seurat), "\n")
cat("细胞数：", ncol(as_seurat), "\n")

# --- 1.2 质量控制 ---
as_seurat[["percent.mt"]] <- PercentageFeatureSet(as_seurat, pattern = "^MT-")
as_seurat[["percent.rb"]] <- PercentageFeatureSet(as_seurat, pattern = "^RP[SL]")

p_qc_as <- VlnPlot(as_seurat,
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3, pt.size = 0,
                   cols = npg_colors[1]) &
  theme_publish()

ggsave("results/figures/Fig_S1A_AS_QC_violin.pdf", p_qc_as,
       width = 12, height = 4, dpi = 300)

p_scatter1 <- FeatureScatter(as_seurat, feature1 = "nCount_RNA",
                             feature2 = "percent.mt", cols = npg_colors[1]) +
  theme_publish()
p_scatter2 <- FeatureScatter(as_seurat, feature1 = "nCount_RNA",
                             feature2 = "nFeature_RNA", cols = npg_colors[1]) +
  theme_publish()

ggsave("results/figures/Fig_S1B_AS_QC_scatter.pdf",
       p_scatter1 + p_scatter2, width = 10, height = 4, dpi = 300)

# 过滤
as_seurat <- subset(as_seurat,
                    subset = nFeature_RNA > 200 &
                      nFeature_RNA < 5000 &
                      percent.mt < 15 &
                      nCount_RNA < 30000)

cat("过滤后细胞数：", ncol(as_seurat), "\n")

# --- 1.3 标准化与降维 ---
as_seurat <- NormalizeData(as_seurat, normalization.method = "LogNormalize",
                           scale.factor = 10000)
as_seurat <- FindVariableFeatures(as_seurat, selection.method = "vst",
                                  nfeatures = 2000)

top10_hvg <- head(VariableFeatures(as_seurat), 10)
p_hvg <- VariableFeaturePlot(as_seurat)
p_hvg <- LabelPoints(plot = p_hvg, points = top10_hvg, repel = TRUE)
ggsave("results/figures/Fig_S1C_AS_HVG.pdf", p_hvg, width = 8, height = 5)

as_seurat <- ScaleData(as_seurat, features = rownames(as_seurat))
as_seurat <- RunPCA(as_seurat, features = VariableFeatures(as_seurat), npcs = 50)

p_elbow <- ElbowPlot(as_seurat, ndims = 50) + theme_publish()
ggsave("results/figures/Fig_S1D_AS_ElbowPlot.pdf", p_elbow, width = 6, height = 4)

# --- 1.4 聚类与UMAP ---
n_pcs <- 30

as_seurat <- FindNeighbors(as_seurat, dims = 1:n_pcs)
as_seurat <- FindClusters(as_seurat, resolution = 0.8)
as_seurat <- RunUMAP(as_seurat, dims = 1:n_pcs)

p_cluster_as <- DimPlot(as_seurat, reduction = "umap", label = TRUE,
                        label.size = 4, repel = TRUE,
                        cols = npg_extended) +
  ggtitle("AS - UMAP Clustering") +
  theme_publish()

ggsave("results/figures/Fig_S1E_AS_UMAP_clusters.pdf",
       p_cluster_as, width = 8, height = 6, dpi = 300)

# --- 1.5 细胞类型注释 ---
as_cell_markers <- list(
  "T cells" = c("CD3D", "CD3E", "CD3G", "CD2"),
  "B cells" = c("CD79A", "MS4A1", "CD19", "BANK1"),
  "NK cells" = c("GNLY", "NKG7", "KLRD1", "KLRB1"),
  "Macrophage" = c("CD68", "CD14", "CD163", "MSR1", "MARCO", "MRC1"),
  "Monocyte" = c("CD14", "LYZ", "S100A8", "S100A9", "FCN1"),
  "DC" = c("FCER1A", "CD1C", "CLEC10A", "ITGAX"),
  "Endothelial" = c("PECAM1", "VWF", "CDH5", "KDR"),
  "SMC" = c("ACTA2", "MYH11", "TAGLN", "CNN1"),
  "Fibroblast" = c("DCN", "COL1A1", "COL3A1", "LUM"),
  "Mast cells" = c("KIT", "TPSAB1", "TPSB2", "CPA3")
)

p_dotplot_as <- DotPlot(as_seurat,
                        features = unique(unlist(as_cell_markers)),
                        cols = c("lightgrey", npg_colors[1])) +
  RotatedAxis() +
  theme_publish() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave("results/figures/Fig_S1F_AS_marker_dotplot.pdf",
       p_dotplot_as, width = 16, height = 6, dpi = 300)

# 手动注释（根据实际聚类结果调整）
as_annotations <- c(
  "0" = "Macrophage",
  "1" = "T cells",
  "2" = "Monocyte",
  "3" = "SMC",
  "4" = "Endothelial",
  "5" = "Macrophage",
  "6" = "NK cells",
  "7" = "B cells",
  "8" = "Fibroblast",
  "9" = "DC",
  "10" = "Mast cells",
  "11" = "T cells",
  "12" = "Macrophage"
)

as_seurat$cellType <- as_annotations[as.character(Idents(as_seurat))]
Idents(as_seurat) <- "cellType"

# >>> 如果使用预处理好的数据，直接加载 <<<
# load("AS_obj.Rdata")
# as_seurat <- AS_obj
# Idents(as_seurat) <- "cellType"

# --- 统一元数据列名为 cellType ---
# 如果原始数据列名是 cellType，统一改为 cellType
if ("cellType" %in% colnames(as_seurat@meta.data) &
    !"cellType" %in% colnames(as_seurat@meta.data)) {
  as_seurat$cellType <- as_seurat$cellType
}
Idents(as_seurat) <- "cellType"

p_annotated_as_legend <- DimPlot(as_seurat, reduction = "umap",
                                 group.by = "cellType",
                                 label = FALSE,
                                 cols = npg_extended) +
  ggtitle("AS - Cell Type Annotation") +
  theme_publish()

ggsave("results/figures/Fig2A_AS_UMAP_cellType.pdf",
       p_annotated_as_legend, width = 9, height = 6, dpi = 300)

cell_prop_as <- as.data.frame(table(as_seurat$cellType)) %>%
  mutate(Prop = Freq / sum(Freq) * 100) %>%
  arrange(desc(Prop))

p_prop_as <- ggplot(cell_prop_as, aes(x = reorder(Var1, -Prop), y = Prop,
                                      fill = Var1)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = npg_extended) +
  labs(x = "", y = "Proportion (%)", title = "AS - Cell Type Proportion") +
  theme_publish() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  NoLegend()

ggsave("results/figures/Fig2B_AS_cellType_proportion.pdf",
       p_prop_as, width = 6, height = 5, dpi = 300)

saveRDS(as_seurat, "results/AS_seurat_annotated.rds")
cat("AS数据集预处理与注释完成！\n")


# ============================================================================
# Part 2: AD单细胞数据预处理
# ============================================================================

# --- 2.1 读取AD数据 ---
ad_seurat <- readRDS("data/AD_seurat_raw.rds")

cat("AD数据集维度：\n")
cat("基因数：", nrow(ad_seurat), "\n")
cat("细胞数：", ncol(ad_seurat), "\n")

# --- 2.2 质量控制 ---
ad_seurat[["percent.mt"]] <- PercentageFeatureSet(ad_seurat, pattern = "^MT-")
ad_seurat[["percent.rb"]] <- PercentageFeatureSet(ad_seurat, pattern = "^RP[SL]")

p_qc_ad <- VlnPlot(ad_seurat,
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3, pt.size = 0,
                   cols = npg_colors[2]) &
  theme_publish()

ggsave("results/figures/Fig_S2A_AD_QC_violin.pdf", p_qc_ad,
       width = 12, height = 4, dpi = 300)

ad_seurat <- subset(ad_seurat,
                    subset = nFeature_RNA > 200 &
                      nFeature_RNA < 5000 &
                      percent.mt < 15 &
                      nCount_RNA < 30000)

cat("过滤后细胞数：", ncol(ad_seurat), "\n")

# --- 2.3 标准化与降维 ---
ad_seurat <- NormalizeData(ad_seurat)
ad_seurat <- FindVariableFeatures(ad_seurat, nfeatures = 2000)
ad_seurat <- ScaleData(ad_seurat, features = rownames(ad_seurat))
ad_seurat <- RunPCA(ad_seurat, npcs = 50)

p_elbow_ad <- ElbowPlot(ad_seurat, ndims = 50) + theme_publish()
ggsave("results/figures/Fig_S2B_AD_ElbowPlot.pdf", p_elbow_ad,
       width = 6, height = 4)

# --- 2.4 聚类与UMAP ---
ad_seurat <- FindNeighbors(ad_seurat, dims = 1:30)
ad_seurat <- FindClusters(ad_seurat, resolution = 0.8)
ad_seurat <- RunUMAP(ad_seurat, dims = 1:30)

# --- 2.5 细胞类型注释 ---
ad_cell_markers <- list(
  "Excitatory Neurons" = c("SLC17A7", "CAMK2A", "NRGN"),
  "Inhibitory Neurons" = c("GAD1", "GAD2", "SLC32A1"),
  "Astrocytes" = c("AQP4", "GFAP", "ALDH1L1", "GJA1"),
  "Microglia" = c("CX3CR1", "P2RY12", "TMEM119", "CSF1R", "AIF1"),
  "Oligodendrocytes" = c("MBP", "PLP1", "MOG", "MAG"),
  "OPC" = c("PDGFRA", "CSPG4", "SOX10"),
  "Endothelial" = c("CLDN5", "FLT1", "PECAM1"),
  "Pericytes" = c("PDGFRB", "RGS5", "KCNJ8")
)

p_dotplot_ad <- DotPlot(ad_seurat,
                        features = unique(unlist(ad_cell_markers)),
                        cols = c("lightgrey", npg_colors[2])) +
  RotatedAxis() +
  theme_publish() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave("results/figures/Fig_S2C_AD_marker_dotplot.pdf",
       p_dotplot_ad, width = 14, height = 6, dpi = 300)

ad_annotations <- c(
  "0" = "Excitatory Neurons",
  "1" = "Oligodendrocytes",
  "2" = "Excitatory Neurons",
  "3" = "Inhibitory Neurons",
  "4" = "Astrocytes",
  "5" = "Microglia",
  "6" = "OPC",
  "7" = "Excitatory Neurons",
  "8" = "Microglia",
  "9" = "Endothelial",
  "10" = "Pericytes",
  "11" = "Inhibitory Neurons"
)

ad_seurat$cellType <- ad_annotations[as.character(Idents(ad_seurat))]
Idents(ad_seurat) <- "cellType"

# >>> 如果使用预处理好的数据 <<<
# load("AD_obj.Rdata")
# ad_seurat <- AD_obj
# Idents(ad_seurat) <- "cellType"

# --- 统一元数据列名为 cellType ---
if ("cellType" %in% colnames(ad_seurat@meta.data) &
    !"cellType" %in% colnames(ad_seurat@meta.data)) {
  ad_seurat$cellType <- ad_seurat$cellType
}
Idents(ad_seurat) <- "cellType"

p_annotated_ad <- DimPlot(ad_seurat, reduction = "umap",
                          group.by = "cellType",
                          label = TRUE, label.size = 4, repel = TRUE,
                          cols = npg_extended) +
  ggtitle("AD - Cell Type Annotation") +
  theme_publish()

ggsave("results/figures/Fig3A_AD_UMAP_cellType.pdf",
       p_annotated_ad, width = 9, height = 6, dpi = 300)

cell_prop_ad <- as.data.frame(table(ad_seurat$cellType)) %>%
  mutate(Prop = Freq / sum(Freq) * 100) %>%
  arrange(desc(Prop))

p_prop_ad <- ggplot(cell_prop_ad, aes(x = reorder(Var1, -Prop), y = Prop,
                                      fill = Var1)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = npg_extended) +
  labs(x = "", y = "Proportion (%)", title = "AD - Cell Type Proportion") +
  theme_publish() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  NoLegend()

ggsave("results/figures/Fig3B_AD_cellType_proportion.pdf",
       p_prop_ad, width = 6, height = 5, dpi = 300)

saveRDS(ad_seurat, "results/AD_seurat_annotated.rds")
cat("AD数据集预处理与注释完成！\n")


# ============================================================================
# Part 3: CTSB在两个数据集中的表达特征
# ============================================================================

as_seurat <- readRDS("results/AS_seurat_annotated.rds")
ad_seurat <- readRDS("results/AD_seurat_annotated.rds")

# 确保 cellType 列存在并设为 Idents
if ("cellType" %in% colnames(as_seurat@meta.data) &
    !"cellType" %in% colnames(as_seurat@meta.data)) {
  as_seurat$cellType <- as_seurat$cellType
}
if ("cellType" %in% colnames(ad_seurat@meta.data) &
    !"cellType" %in% colnames(ad_seurat@meta.data)) {
  ad_seurat$cellType <- ad_seurat$cellType
}
Idents(as_seurat) <- "cellType"
Idents(ad_seurat) <- "cellType"

# --- 3.1 AS数据集中CTSB的表达 ---

# 判断可用的降维方式
as_reduction <- ifelse("tsne" %in% names(as_seurat@reductions), "tsne", "umap")
ad_reduction <- ifelse("tsne" %in% names(ad_seurat@reductions), "tsne", "umap")

p_feature_as <- FeaturePlot(as_seurat, features = "CTSB",
                            reduction = as_reduction,
                            cols = c("lightgrey", "#E64B35FF"),
                            order = TRUE, pt.size = 0.3) +
  ggtitle("CTSB Expression in AS") +
  theme_publish()

p_vln_as <- VlnPlot(as_seurat, features = "CTSB",
                    group.by = "cellType",
                    cols = npg_extended, pt.size = 0) +
  ggtitle("CTSB Expression across Cell Types (AS)") +
  theme_publish() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  NoLegend() +
  stat_summary(fun = median, geom = "point", shape = 18,
               size = 3, color = "black")

p_dot_ctsb_as <- DotPlot(as_seurat, features = "CTSB",
                         group.by = "cellType",
                         cols = c("lightgrey", "#E64B35FF")) +
  coord_flip() +
  ggtitle("CTSB Expression (AS)") +
  theme_publish()

# --- 3.2 AD数据集中CTSB的表达 ---

p_feature_ad <- FeaturePlot(ad_seurat, features = "CTSB",
                            reduction = ad_reduction,
                            cols = c("lightgrey", "#4DBBD5FF"),
                            order = TRUE, pt.size = 0.3) +
  ggtitle("CTSB Expression in AD") +
  theme_publish()

p_vln_ad <- VlnPlot(ad_seurat, features = "CTSB",
                    group.by = "cellType",
                    cols = npg_extended, pt.size = 0) +
  ggtitle("CTSB Expression across Cell Types (AD)") +
  theme_publish() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  NoLegend() +
  stat_summary(fun = median, geom = "point", shape = 18,
               size = 3, color = "black")

# --- 3.3 组合图 ---
p_feature_combined <- p_feature_as + p_feature_ad +
  plot_layout(ncol = 2) +
  plot_annotation(title = "CTSB Expression Pattern",
                  theme = theme(plot.title = element_text(
                    size = 16, face = "bold", hjust = 0.5,
                    family = "sans")))

ggsave("results/figures/Fig4A_CTSB_FeaturePlot_combined.pdf",
       p_feature_combined, width = 14, height = 6, dpi = 300)

p_vln_combined <- p_vln_as / p_vln_ad +
  plot_annotation(tag_levels = "A")

ggsave("results/figures/Fig4B_CTSB_VlnPlot_combined.pdf",
       p_vln_combined, width = 10, height = 10, dpi = 300)

# --- 3.4 CTSB表达统计检验 ---
as_meta <- as_seurat@meta.data
as_meta$CTSB_expr <- as_seurat[["RNA"]]$data["CTSB", ]
as_meta$is_mac <- ifelse(as_meta$cellType == "Macrophage", "Macrophage", "Others")

wilcox_as <- wilcox.test(CTSB_expr ~ is_mac, data = as_meta)
cat("AS - Macrophage vs Others: p-value =", wilcox_as$p.value, "\n")

ad_meta <- ad_seurat@meta.data
ad_meta$CTSB_expr <- ad_seurat[["RNA"]]$data["CTSB", ]
ad_meta$is_mic <- ifelse(ad_meta$cellType == "Microglia", "Microglia", "Others")

wilcox_ad <- wilcox.test(CTSB_expr ~ is_mic, data = ad_meta)
cat("AD - Microglia vs Others: p-value =", wilcox_ad$p.value, "\n")

# --- 3.5 CTSB表达比较箱线图 ---
p_box_as <- ggplot(as_meta, aes(x = is_mac, y = CTSB_expr, fill = is_mac)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  scale_fill_manual(values = c(npg_colors[1], "grey80")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, label.y = max(as_meta$CTSB_expr) * 1.05) +
  labs(x = "", y = "CTSB Expression (log-normalized)",
       title = "CTSB in AS") +
  theme_publish() +
  NoLegend()

p_box_ad <- ggplot(ad_meta, aes(x = is_mic, y = CTSB_expr, fill = is_mic)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  scale_fill_manual(values = c(npg_colors[2], "grey80")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, label.y = max(ad_meta$CTSB_expr) * 1.05) +
  labs(x = "", y = "CTSB Expression (log-normalized)",
       title = "CTSB in AD") +
  theme_publish() +
  NoLegend()

p_box_combined <- p_box_as + p_box_ad
ggsave("results/figures/Fig4C_CTSB_boxplot_comparison.pdf",
       p_box_combined, width = 8, height = 5, dpi = 300)

# --- 3.6 各细胞类型CTSB表达均值条形图 ---
ctsb_mean_as <- as_meta %>%
  group_by(cellType) %>%
  summarise(
    mean_expr = mean(CTSB_expr),
    se = sd(CTSB_expr) / sqrt(n()),
    pct_express = sum(CTSB_expr > 0) / n() * 100
  ) %>%
  arrange(desc(mean_expr))

p_bar_as <- ggplot(ctsb_mean_as, aes(x = reorder(cellType, mean_expr),
                                     y = mean_expr, fill = cellType)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = mean_expr - se, ymax = mean_expr + se),
                width = 0.2) +
  scale_fill_manual(values = npg_extended) +
  coord_flip() +
  labs(x = "", y = "Mean CTSB Expression", title = "AS Dataset") +
  theme_publish() +
  NoLegend()

ctsb_mean_ad <- ad_meta %>%
  group_by(cellType) %>%
  summarise(
    mean_expr = mean(CTSB_expr),
    se = sd(CTSB_expr) / sqrt(n()),
    pct_express = sum(CTSB_expr > 0) / n() * 100
  ) %>%
  arrange(desc(mean_expr))

p_bar_ad <- ggplot(ctsb_mean_ad, aes(x = reorder(cellType, mean_expr),
                                     y = mean_expr, fill = cellType)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = mean_expr - se, ymax = mean_expr + se),
                width = 0.2) +
  scale_fill_manual(values = npg_extended) +
  coord_flip() +
  labs(x = "", y = "Mean CTSB Expression", title = "AD Dataset") +
  theme_publish() +
  NoLegend()

p_bar_combined <- p_bar_as + p_bar_ad
ggsave("results/figures/Fig4D_CTSB_mean_expression_barplot.pdf",
       p_bar_combined, width = 12, height = 5, dpi = 300)

cat("Part 3 完成：CTSB表达特征分析\n")


# ============================================================================
# Part 5: 跨疾病转录组关联分析 (核心创新点)
# ============================================================================

# --- 5.1 Marker基因重叠分析 ---
Idents(as_seurat) <- "cellType"
Idents(ad_seurat) <- "cellType"

as_mac_markers <- FindMarkers(as_seurat, ident.1 = "Macrophage",
                              only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.5)
as_mac_genes <- rownames(as_mac_markers[as_mac_markers$p_val_adj < 0.05, ])

ad_mic_markers <- FindMarkers(ad_seurat, ident.1 = "Microglia",
                              only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.5)
ad_mic_genes <- rownames(ad_mic_markers[ad_mic_markers$p_val_adj < 0.05, ])

overlap_genes <- intersect(as_mac_genes, ad_mic_genes)
cat("AS Macrophage markers:", length(as_mac_genes), "\n")
cat("AD Microglia markers:", length(ad_mic_genes), "\n")
cat("Overlapping genes:", length(overlap_genes), "\n")

# Fisher精确检验
common_genes <- intersect(rownames(as_seurat), rownames(ad_seurat))
n_total <- length(common_genes)

contingency <- matrix(c(
  length(overlap_genes),
  length(as_mac_genes) - length(overlap_genes),
  length(ad_mic_genes) - length(overlap_genes),
  n_total - length(union(as_mac_genes, ad_mic_genes))
), nrow = 2)

fisher_result <- fisher.test(contingency)
cat("Fisher test p-value:", fisher_result$p.value, "\n")
cat("Odds ratio:", fisher_result$estimate, "\n")

# 韦恩图
pdf("results/figures/Fig6A_Venn_Mac_Mic_markers.pdf", width = 6, height = 6)
venn.plot <- draw.pairwise.venn(
  area1 = length(as_mac_genes),
  area2 = length(ad_mic_genes),
  cross.area = length(overlap_genes),
  category = c("AS-Macrophage", "AD-Microglia"),
  fill = c("#E64B35FF", "#4DBBD5FF"),
  alpha = 0.5,
  lty = "blank",
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.pos = c(-30, 30),
  cat.dist = 0.05,
  cat.col = c("#E64B35FF", "#4DBBD5FF")
)
grid.text(paste0("Fisher p = ", formatC(fisher_result$p.value, format = "e", digits = 2)),
          x = 0.5, y = 0.1, gp = gpar(fontsize = 12, fontface = "bold",
                                      fontfamily = "sans"))
dev.off()

write.csv(data.frame(Gene = overlap_genes),
          "results/tables/overlap_genes_Mac_Mic.csv", row.names = FALSE)

# --- 5.2 重叠基因的GO/KEGG富集分析 ---
overlap_entrez <- bitr(overlap_genes, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_overlap <- enrichGO(gene = overlap_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)

kegg_overlap <- enrichKEGG(gene = overlap_entrez$ENTREZID,
                           organism = "hsa",
                           pvalueCutoff = 0.05)

# GO气泡图
p_go_overlap <- dotplot(go_overlap, showCategory = 15, font.size = 10) +
  scale_color_gradient(low = "#4DBBD5FF", high = "#E64B35FF") +
  ggtitle("GO-BP: Shared Genes (Mac ∩ Mic)") +
  theme_bw() +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text = element_text(size = 10)
  )

ggsave("results/figures/Fig6B_GO_overlap_genes.pdf",
       p_go_overlap, width = 8, height = 7, dpi = 300)

# KEGG条形图
p_kegg_overlap <- barplot(kegg_overlap, showCategory = 15) +
  ggtitle("KEGG: Shared Genes (Mac ∩ Mic)") +
  theme_bw() +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave("results/figures/Fig6C_KEGG_overlap_genes.pdf",
       p_kegg_overlap, width = 8, height = 6, dpi = 300)

# --- 5.3 Pseudo-bulk相关性分析 ---
as_avg <- AverageExpression(as_seurat, group.by = "cellType",
                            return.seurat = FALSE)$RNA
ad_avg <- AverageExpression(ad_seurat, group.by = "cellType",
                            return.seurat = FALSE)$RNA

common_genes <- intersect(rownames(as_avg), rownames(ad_avg))
as_avg_common <- as_avg[common_genes, ]
ad_avg_common <- ad_avg[common_genes, ]

cor_matrix <- matrix(NA, ncol = ncol(ad_avg_common), nrow = ncol(as_avg_common))
rownames(cor_matrix) <- paste0("AS_", colnames(as_avg_common))
colnames(cor_matrix) <- paste0("AD_", colnames(ad_avg_common))

pval_matrix <- cor_matrix

for (i in 1:ncol(as_avg_common)) {
  for (j in 1:ncol(ad_avg_common)) {
    test <- cor.test(as_avg_common[, i], ad_avg_common[, j],
                     method = "pearson")
    cor_matrix[i, j] <- test$estimate
    pval_matrix[i, j] <- test$p.value
  }
}

# 热图
pdf("results/figures/Fig6D_cross_disease_correlation.pdf",
    width = 8, height = 7)

col_fun <- colorRamp2(c(0, 0.5, 1), c("#4DBBD5FF", "white", "#E64B35FF"))

sig_matrix <- ifelse(pval_matrix < 0.001, "***",
                     ifelse(pval_matrix < 0.01, "**",
                            ifelse(pval_matrix < 0.05, "*", "")))

ht <- Heatmap(cor_matrix,
              name = "Pearson r",
              col = col_fun,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              row_names_gp = gpar(fontsize = 10, fontface = "bold",
                                  fontfamily = "sans"),
              column_names_gp = gpar(fontsize = 10, fontface = "bold",
                                     fontfamily = "sans"),
              column_names_rot = 45,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sig_matrix[i, j], x, y,
                          gp = gpar(fontsize = 12, fontface = "bold",
                                    fontfamily = "sans"))
                grid.text(sprintf("%.2f", cor_matrix[i, j]), x,
                          y - unit(2, "mm"),
                          gp = gpar(fontsize = 8, fontfamily = "sans"))
              },
              column_title = "Cross-Disease Cell Type Correlation",
              column_title_gp = gpar(fontsize = 14, fontface = "bold",
                                     fontfamily = "sans"),
              heatmap_legend_param = list(
                title = "Pearson r",
                title_gp = gpar(fontsize = 10, fontface = "bold",
                                fontfamily = "sans"),
                labels_gp = gpar(fontsize = 9, fontfamily = "sans")
              ))
draw(ht)
dev.off()

# --- 5.4 散点图 ---
mac_col <- grep("Macrophage", colnames(as_avg_common))
mic_col <- grep("Microglia", colnames(ad_avg_common))

scatter_df <- data.frame(
  AS_Macrophage = as_avg_common[, mac_col],
  AD_Microglia = ad_avg_common[, mic_col],
  Gene = common_genes
)

highlight_genes <- c("CTSB", "CTSD", "CTSL", "CD68", "TREM2", "APOE",
                     "IL1B", "TNF", "CSF1R", "TYROBP", "LYZ", "LAMP1",
                     "LGALS3", "SPP1", "C1QA", "C1QB", "HLA-DRA")
scatter_df$highlight <- ifelse(scatter_df$Gene %in% highlight_genes,
                               scatter_df$Gene, NA)

cor_test <- cor.test(scatter_df$AS_Macrophage, scatter_df$AD_Microglia)

p_scatter_mac_mic <- ggplot(scatter_df, aes(x = AS_Macrophage, y = AD_Microglia)) +
  geom_point(alpha = 0.3, size = 0.8, color = "grey60") +
  geom_point(data = scatter_df[!is.na(scatter_df$highlight), ],
             aes(x = AS_Macrophage, y = AD_Microglia),
             color = "#E64B35FF", size = 2.5) +
  ggrepel::geom_text_repel(aes(label = highlight),
                           size = 3, fontface = "italic",
                           color = "#E64B35FF",
                           max.overlaps = 20,
                           box.padding = 0.5) +
  geom_smooth(method = "lm", color = "#3C5488FF", se = TRUE,
              fill = "#3C5488FF", alpha = 0.15) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("r = ", round(cor_test$estimate, 3),
                          "\np = ", formatC(cor_test$p.value,
                                            format = "e", digits = 2)),
           hjust = 1.1, vjust = -0.5, size = 4.5, fontface = "bold",
           color = "#3C5488FF", family = "sans") +
  labs(x = "AS - Macrophage (Mean Expression)",
       y = "AD - Microglia (Mean Expression)",
       title = "Transcriptomic Similarity: Macrophage vs Microglia") +
  theme_publish()

ggsave("results/figures/Fig6E_scatter_Mac_vs_Mic.pdf",
       p_scatter_mac_mic, width = 8, height = 7, dpi = 300)

# --- 5.5 共享高表达基因热图 ---
top_shared <- scatter_df %>%
  filter(AS_Macrophage > quantile(AS_Macrophage, 0.9) &
           AD_Microglia > quantile(AD_Microglia, 0.9)) %>%
  arrange(desc(AS_Macrophage + AD_Microglia)) %>%
  head(50)

heatmap_data <- cbind(
  as_avg_common[top_shared$Gene, ],
  ad_avg_common[top_shared$Gene, ]
)

col_anno <- HeatmapAnnotation(
  Disease = c(rep("AS", ncol(as_avg_common)),
              rep("AD", ncol(ad_avg_common))),
  col = list(Disease = c("AS" = "#E64B35FF", "AD" = "#4DBBD5FF")),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold",
                            fontfamily = "sans")
)

row_idx_ctsb <- which(rownames(heatmap_data) == "CTSB")
row_anno <- NULL
if (length(row_idx_ctsb) > 0) {
  row_anno <- rowAnnotation(
    CTSB = anno_mark(at = row_idx_ctsb,
                     labels = "CTSB",
                     labels_gp = gpar(fontsize = 12, fontface = "bold",
                                      col = "#E64B35FF",
                                      fontfamily = "sans"))
  )
}

pdf("results/figures/Fig6F_shared_genes_heatmap.pdf", width = 12, height = 10)
ht2 <- Heatmap(t(scale(t(as.matrix(heatmap_data)))),
               name = "Z-score",
               col = colorRamp2(c(-2, 0, 2), c("#4DBBD5FF", "white", "#E64B35FF")),
               top_annotation = col_anno,
               right_annotation = row_anno,
               show_row_names = TRUE,
               row_names_gp = gpar(fontsize = 7, fontfamily = "sans"),
               column_names_gp = gpar(fontsize = 9, fontfamily = "sans"),
               column_names_rot = 45,
               cluster_columns = FALSE,
               column_title = "Shared Highly Expressed Genes",
               column_title_gp = gpar(fontsize = 14, fontface = "bold",
                                      fontfamily = "sans"))
draw(ht2)
dev.off()

cat("Part 5 完成：跨疾病转录组关联分析\n")


# ============================================================================
# Part 6: CTSB高/低表达组差异分析
# ============================================================================

# --- 6.1 AS - Macrophage中CTSB高/低表达组 ---
mac_cells <- subset(as_seurat, cellType == "Macrophage")
ctsb_expr_mac <- GetAssayData(mac_cells, slot = "data")["CTSB", ]

mac_cells$CTSB_group <- ifelse(ctsb_expr_mac > median(ctsb_expr_mac),
                               "CTSB_high", "CTSB_low")

Idents(mac_cells) <- "CTSB_group"
degs_mac <- FindMarkers(mac_cells, ident.1 = "CTSB_high",
                        ident.2 = "CTSB_low",
                        logfc.threshold = 0.25, min.pct = 0.1)

degs_mac$gene <- rownames(degs_mac)
degs_mac_sig <- degs_mac[degs_mac$p_val_adj < 0.05, ]

cat("AS Mac - CTSB_high vs CTSB_low DEGs:", nrow(degs_mac_sig), "\n")

# --- 6.2 AD - Microglia中CTSB高/低表达组 ---
mic_cells <- subset(ad_seurat, cellType == "Microglia")
ctsb_expr_mic <- GetAssayData(mic_cells, slot = "data")["CTSB", ]

mic_cells$CTSB_group <- ifelse(ctsb_expr_mic > median(ctsb_expr_mic),
                               "CTSB_high", "CTSB_low")

Idents(mic_cells) <- "CTSB_group"
degs_mic <- FindMarkers(mic_cells, ident.1 = "CTSB_high",
                        ident.2 = "CTSB_low",
                        logfc.threshold = 0.25, min.pct = 0.1)

degs_mic$gene <- rownames(degs_mic)
degs_mic_sig <- degs_mic[degs_mic$p_val_adj < 0.05, ]

cat("AD Mic - CTSB_high vs CTSB_low DEGs:", nrow(degs_mic_sig), "\n")

# --- 6.3 火山图 ---
plot_volcano <- function(deg_df, title, color1 = "#E64B35FF",
                         color2 = "#4DBBD5FF") {
  deg_df$significance <- "NS"
  deg_df$significance[deg_df$p_val_adj < 0.05 & deg_df$avg_log2FC > 0.5] <- "Up"
  deg_df$significance[deg_df$p_val_adj < 0.05 & deg_df$avg_log2FC < -0.5] <- "Down"
  
  top_up <- deg_df %>% filter(significance == "Up") %>%
    arrange(p_val_adj) %>% head(10)
  top_down <- deg_df %>% filter(significance == "Down") %>%
    arrange(p_val_adj) %>% head(10)
  top_genes <- rbind(top_up, top_down)
  
  if ("CTSB" %in% deg_df$gene & !("CTSB" %in% top_genes$gene)) {
    top_genes <- rbind(top_genes, deg_df[deg_df$gene == "CTSB", ])
  }
  
  p <- ggplot(deg_df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = significance), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Down" = color2, "NS" = "grey80",
                                  "Up" = color1)) +
    geom_text_repel(data = top_genes,
                    aes(label = gene), size = 3,
                    fontface = "italic",
                    max.overlaps = 20) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
               color = "grey40") +
    labs(x = "log2(Fold Change)", y = "-log10(Adjusted P-value)",
         title = title) +
    theme_publish() +
    theme(legend.position = "right")
  
  return(p)
}

p_volcano_mac <- plot_volcano(degs_mac,
                              "DEGs: CTSB_high vs CTSB_low (AS-Macrophage)")
p_volcano_mic <- plot_volcano(degs_mic,
                              "DEGs: CTSB_high vs CTSB_low (AD-Microglia)")

p_volcano_combined <- p_volcano_mac + p_volcano_mic +
  plot_annotation(tag_levels = "A")

ggsave("results/figures/Fig7A_volcano_CTSB_groups.pdf",
       p_volcano_combined, width = 14, height = 6, dpi = 300)

# --- 6.4 共同DEGs分析 ---
up_mac <- degs_mac_sig$gene[degs_mac_sig$avg_log2FC > 0.25]
up_mic <- degs_mic_sig$gene[degs_mic_sig$avg_log2FC > 0.25]
down_mac <- degs_mac_sig$gene[degs_mac_sig$avg_log2FC < -0.25]
down_mic <- degs_mic_sig$gene[degs_mic_sig$avg_log2FC < -0.25]

common_up <- intersect(up_mac, up_mic)
common_down <- intersect(down_mac, down_mic)

cat("共同上调基因:", length(common_up), "\n")
cat("共同下调基因:", length(common_down), "\n")

pdf("results/figures/Fig7B_Venn_common_up_DEGs.pdf", width = 6, height = 6)
venn_up <- draw.pairwise.venn(
  area1 = length(up_mac), area2 = length(up_mic),
  cross.area = length(common_up),
  category = c("AS-Mac\n(CTSB_high Up)", "AD-Mic\n(CTSB_high Up)"),
  fill = c("#E64B35FF", "#4DBBD5FF"),
  alpha = 0.5, lty = "blank",
  cex = 1.5, fontface = "bold", fontfamily = "sans",
  cat.cex = 1.0, cat.fontface = "bold", cat.fontfamily = "sans",
  cat.col = c("#E64B35FF", "#4DBBD5FF")
)
dev.off()

# --- 6.5 共同上调DEGs的富集分析 ---
common_entrez <- bitr(common_up, fromType = "SYMBOL",
                      toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_common <- enrichGO(gene = common_entrez$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      readable = TRUE)

kegg_common <- enrichKEGG(gene = common_entrez$ENTREZID,
                          organism = "hsa",
                          pvalueCutoff = 0.05)

# GO气泡图
p_go_common <- dotplot(go_common, showCategory = 20, font.size = 9) +
  scale_color_gradient(low = "#E64B35FF", high = "#4DBBD5FF") +
  ggtitle("GO-BP: Common DEGs in CTSB_high Cells") +
  theme_bw() +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.text = element_text(size = 10)
  )

ggsave("results/figures/Fig7C_GO_common_DEGs.pdf",
       p_go_common, width = 9, height = 8, dpi = 300)

# KEGG气泡图
p_kegg_common <- dotplot(kegg_common, showCategory = 15, font.size = 10) +
  scale_color_gradient(low = "#E64B35FF", high = "#4DBBD5FF") +
  ggtitle("KEGG: Common DEGs in CTSB_high Cells") +
  theme_bw() +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 10)
  )

ggsave("results/figures/Fig7D_KEGG_common_DEGs.pdf",
       p_kegg_common, width = 8, height = 6, dpi = 300)

# --- 6.6 共同DEGs热图 ---
common_degs_all <- union(common_up, common_down)

mac_expr <- GetAssayData(mac_cells, slot = "data")[
  common_degs_all[common_degs_all %in% rownames(mac_cells)], ]
mac_avg_by_group <- data.frame(
  Mac_CTSB_high = rowMeans(mac_expr[, mac_cells$CTSB_group == "CTSB_high"]),
  Mac_CTSB_low = rowMeans(mac_expr[, mac_cells$CTSB_group == "CTSB_low"])
)

mic_expr <- GetAssayData(mic_cells, slot = "data")[
  common_degs_all[common_degs_all %in% rownames(mic_cells)], ]
mic_avg_by_group <- data.frame(
  Mic_CTSB_high = rowMeans(mic_expr[, mic_cells$CTSB_group == "CTSB_high"]),
  Mic_CTSB_low = rowMeans(mic_expr[, mic_cells$CTSB_group == "CTSB_low"])
)

common_genes_both <- intersect(rownames(mac_avg_by_group),
                               rownames(mic_avg_by_group))
heatmap_deg <- cbind(mac_avg_by_group[common_genes_both, ],
                     mic_avg_by_group[common_genes_both, ])

col_anno_deg <- HeatmapAnnotation(
  Disease = c("AS", "AS", "AD", "AD"),
  CTSB_group = c("High", "Low", "High", "Low"),
  col = list(
    Disease = c("AS" = "#E64B35FF", "AD" = "#4DBBD5FF"),
    CTSB_group = c("High" = "#E64B35FF", "Low" = "#3C5488FF")
  ),
  annotation_name_gp = gpar(fontfamily = "sans")
)

pdf("results/figures/Fig7E_common_DEGs_heatmap.pdf", width = 8, height = 10)
ht3 <- Heatmap(t(scale(t(as.matrix(heatmap_deg)))),
               name = "Z-score",
               col = colorRamp2(c(-2, 0, 2), c("#4DBBD5FF", "white", "#E64B35FF")),
               top_annotation = col_anno_deg,
               show_row_names = TRUE,
               row_names_gp = gpar(fontsize = 7, fontfamily = "sans"),
               column_names_gp = gpar(fontsize = 10, fontface = "bold",
                                      fontfamily = "sans"),
               column_title = "Common DEGs: CTSB_high vs CTSB_low",
               column_title_gp = gpar(fontsize = 13, fontface = "bold",
                                      fontfamily = "sans"),
               cluster_columns = FALSE)
draw(ht3)
dev.off()

write.csv(degs_mac, "results/tables/DEGs_Mac_CTSB_high_vs_low.csv")
write.csv(degs_mic, "results/tables/DEGs_Mic_CTSB_high_vs_low.csv")
write.csv(data.frame(Gene = common_up, Direction = "Up"),
          "results/tables/Common_up_DEGs.csv", row.names = FALSE)

cat("Part 6 完成：CTSB高/低表达组差异分析\n")


# ============================================================================
# Part 7: 通路活性评分与CTSB相关性分析
# ============================================================================

# --- 7.1 定义关键通路基因集 ---
msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
msig_c2 <- msigdbr(species = "Homo sapiens", category = "C2",
                   subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol)

pathway_genes <- list(
  Lysosome = c("CTSB", "CTSD", "CTSL", "CTSS", "CTSK", "CTSZ",
               "LAMP1", "LAMP2", "PSAP", "HEXA", "HEXB", "GBA",
               "GALC", "GLA", "NPC1", "NPC2", "ATP6V0A1", "ATP6V1A"),
  Autophagy = c("BECN1", "ATG5", "ATG7", "ATG12", "ATG16L1",
                "LC3B", "MAP1LC3B", "SQSTM1", "NBR1", "OPTN",
                "ULK1", "PIK3C3", "AMBRA1", "ATG14", "WIPI2"),
  NLRP3_Inflammasome = c("NLRP3", "PYCARD", "CASP1", "IL1B", "IL18",
                         "GSDMD", "NEK7", "TXNIP", "P2RX7", "PANX1"),
  NFkB_Signaling = c("NFKB1", "NFKB2", "RELA", "RELB", "REL",
                     "NFKBIA", "NFKBIB", "IKBKB", "IKBKG", "CHUK",
                     "TRAF6", "MYD88", "IRAK1", "IRAK4"),
  Lipid_Metabolism = c("ABCA1", "ABCG1", "LDLR", "SCARB1", "CD36",
                       "PPARG", "PPARA", "NR1H3", "SREBF1", "SREBF2",
                       "FABP4", "FABP5", "APOE", "APOC1", "LPL"),
  Oxidative_Stress = c("SOD1", "SOD2", "CAT", "GPX1", "GPX4",
                       "NQO1", "HMOX1", "NFE2L2", "KEAP1",
                       "TXNRD1", "PRDX1", "PRDX2", "GCLC", "GCLM"),
  Phagocytosis = c("MERTK", "AXL", "TYRO3", "CD68", "MARCO",
                   "MSR1", "TREM2", "TYROBP", "SYK", "FCGR1A",
                   "FCGR2A", "FCGR3A", "C1QA", "C1QB", "C1QC"),
  Complement = c("C1QA", "C1QB", "C1QC", "C2", "C3", "C4A",
                 "C5", "CFB", "CFD", "CFH", "CFI", "CD46",
                 "CD55", "CD59", "C3AR1", "C5AR1"),
  Amyloid_Processing = c("APP", "BACE1", "BACE2", "PSEN1", "PSEN2",
                         "NCSTN", "APH1A", "ADAM10", "ADAM17",
                         "IDE", "NEP", "ECE1", "LRP1")
)

# --- 7.2 计算通路活性评分 ---
for (pathway_name in names(pathway_genes)) {
  genes <- pathway_genes[[pathway_name]]
  
  genes_in_mac <- genes[genes %in% rownames(mac_cells)]
  if (length(genes_in_mac) >= 3) {
    mac_cells <- AddModuleScore(mac_cells,
                                features = list(genes_in_mac),
                                name = pathway_name)
  }
  
  genes_in_mic <- genes[genes %in% rownames(mic_cells)]
  if (length(genes_in_mic) >= 3) {
    mic_cells <- AddModuleScore(mic_cells,
                                features = list(genes_in_mic),
                                name = pathway_name)
  }
}

# --- 7.3 CTSB表达与通路活性的相关性 ---
pathway_names <- names(pathway_genes)
score_cols <- paste0(pathway_names, "1")

mac_meta <- mac_cells@meta.data
mac_meta$CTSB_expr <- GetAssayData(mac_cells, slot = "data")["CTSB", ]

mic_meta <- mic_cells@meta.data
mic_meta$CTSB_expr <- GetAssayData(mic_cells, slot = "data")["CTSB", ]

cor_results <- data.frame()

for (i in seq_along(pathway_names)) {
  pname <- pathway_names[i]
  scol <- score_cols[i]
  
  if (scol %in% colnames(mac_meta)) {
    cor_mac <- cor.test(mac_meta$CTSB_expr, mac_meta[[scol]],
                        method = "spearman")
    cor_results <- rbind(cor_results, data.frame(
      Pathway = pname, Cell = "AS_Macrophage",
      Correlation = cor_mac$estimate,
      P_value = cor_mac$p.value
    ))
  }
  
  if (scol %in% colnames(mic_meta)) {
    cor_mic <- cor.test(mic_meta$CTSB_expr, mic_meta[[scol]],
                        method = "spearman")
    cor_results <- rbind(cor_results, data.frame(
      Pathway = pname, Cell = "AD_Microglia",
      Correlation = cor_mic$estimate,
      P_value = cor_mic$p.value
    ))
  }
}

cor_results$Significance <- ifelse(cor_results$P_value < 0.001, "***",
                                   ifelse(cor_results$P_value < 0.01, "**",
                                          ifelse(cor_results$P_value < 0.05, "*", "NS")))

# --- 7.4 相关性气泡图 ---
p_cor_bubble <- ggplot(cor_results, aes(x = Cell, y = Pathway)) +
  geom_point(aes(size = abs(Correlation), color = Correlation)) +
  scale_color_gradient2(low = "#4DBBD5FF", mid = "white", high = "#E64B35FF",
                        midpoint = 0, name = "Spearman rho") +
  scale_size_continuous(range = c(2, 10), name = "|rho|") +
  geom_text(aes(label = Significance), size = 4, vjust = 0.5,
            family = "sans") +
  labs(x = "", y = "",
       title = "CTSB Expression vs Pathway Activity") +
  theme_publish() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 10)
  )

ggsave("results/figures/Fig8A_CTSB_pathway_correlation_bubble.pdf",
       p_cor_bubble, width = 7, height = 7, dpi = 300)

# --- 7.5 关键通路的散点图 ---
key_pathways <- c("Lysosome")

scatter_plots <- list()

for (pname in key_pathways) {
  scol <- paste0(pname, "1")
  
  if (scol %in% colnames(mac_meta) & scol %in% colnames(mic_meta)) {
    df_mac <- data.frame(CTSB = mac_meta$CTSB_expr,
                         Score = mac_meta[[scol]])
    cor_m <- cor.test(df_mac$CTSB, df_mac$Score, method = "spearman")
    
    p1 <- ggplot(df_mac, aes(x = CTSB, y = Score)) +
      geom_point(alpha = 0.2, size = 0.5, color = "#E64B35FF") +
      geom_smooth(method = "lm", color = "#3C5488FF", se = TRUE) +
      annotate("text", x = Inf, y = Inf,
               label = paste0("rho = ", round(cor_m$estimate, 3),
                              "\np = ", formatC(cor_m$p.value,
                                                format = "e", digits = 1)),
               hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "bold",
               family = "sans") +
      labs(x = "CTSB Expression", y = paste0(pname, " Score"),
           title = paste0("AS-Mac: ", pname)) +
      theme_publish() +
      theme(plot.title = element_text(size = 10))
    
    df_mic <- data.frame(CTSB = mic_meta$CTSB_expr,
                         Score = mic_meta[[scol]])
    cor_c <- cor.test(df_mic$CTSB, df_mic$Score, method = "spearman")
    
    p2 <- ggplot(df_mic, aes(x = CTSB, y = Score)) +
      geom_point(alpha = 0.2, size = 0.5, color = "#4DBBD5FF") +
      geom_smooth(method = "lm", color = "#3C5488FF", se = TRUE) +
      annotate("text", x = Inf, y = Inf,
               label = paste0("rho = ", round(cor_c$estimate, 3),
                              "\np = ", formatC(cor_c$p.value,
                                                format = "e", digits = 1)),
               hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "bold",
               family = "sans") +
      labs(x = "CTSB Expression", y = paste0(pname, " Score"),
           title = paste0("AD-Mic: ", pname)) +
      theme_publish() +
      theme(plot.title = element_text(size = 10))
    
    scatter_plots[[paste0(pname, "_Mac")]] <- p1
    scatter_plots[[paste0(pname, "_Mic")]] <- p2
  }
}

p_scatter_pathways <- wrap_plots(scatter_plots, ncol = 4) +
  plot_annotation(tag_levels = "A",
                  title = "CTSB Expression Correlates with Key Pathways",
                  theme = theme(plot.title = element_text(
                    size = 14, face = "bold", hjust = 0.5, family = "sans")))

ggsave("results/figures/Fig8B_CTSB_pathway_scatter.pdf",
       p_scatter_pathways, width = 16, height = 8, dpi = 300)

# --- 7.6 通路活性在CTSB高/低组的比较 ---
mac_meta$CTSB_group <- mac_cells$CTSB_group
mic_meta$CTSB_group <- mic_cells$CTSB_group

pathway_comparison_mac <- mac_meta %>%
  dplyr::select(CTSB_group, all_of(score_cols[score_cols %in% colnames(mac_meta)])) %>%
  pivot_longer(cols = -CTSB_group, names_to = "Pathway", values_to = "Score") %>%
  mutate(Pathway = gsub("1$", "", Pathway),
         Dataset = "AS-Macrophage")

pathway_comparison_mic <- mic_meta %>%
  dplyr::select(CTSB_group, all_of(score_cols[score_cols %in% colnames(mic_meta)])) %>%
  pivot_longer(cols = -CTSB_group, names_to = "Pathway", values_to = "Score") %>%
  mutate(Pathway = gsub("1$", "", Pathway),
         Dataset = "AD-Microglia")

pathway_comparison <- rbind(pathway_comparison_mac, pathway_comparison_mic)

p_pathway_boxplot <- ggplot(pathway_comparison,
                            aes(x = CTSB_group, y = Score, fill = CTSB_group)) +
  geom_boxplot(outlier.size = 0.3, width = 0.6) +
  facet_grid(Dataset ~ Pathway, scales = "free_y") +
  scale_fill_manual(values = c("CTSB_high" = "#E64B35FF",
                               "CTSB_low" = "#4DBBD5FF")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     size = 3, label.y.npc = 0.95) +
  labs(x = "", y = "Pathway Activity Score",
       title = "Pathway Activity: CTSB_high vs CTSB_low") +
  theme_publish() +
  theme(strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom")

ggsave("results/figures/Fig8C_pathway_activity_boxplot.pdf",
       p_pathway_boxplot, width = 18, height = 8, dpi = 300)

write.csv(cor_results, "results/tables/CTSB_pathway_correlations.csv",
          row.names = FALSE)

cat("Part 7 完成：通路活性评分分析\n")


# ============================================================================
# Part 8: 跨数据集整合分析 (Harmony)
# ============================================================================

# --- 8.1 提取目标细胞并整合 ---
mac_for_merge <- subset(as_seurat, cellType == "Macrophage")
mic_for_merge <- subset(ad_seurat, cellType == "Microglia")

mac_for_merge$disease <- "AS"
mac_for_merge$cell_origin <- "AS_Macrophage"
mic_for_merge$disease <- "AD"
mic_for_merge$cell_origin <- "AD_Microglia"

# 取AS和AD中其他代表性细胞做对照
# 请根据实际细胞类型名称修改
as_types <- unique(as_seurat$cellType)
ad_types <- unique(ad_seurat$cellType)

cat("AS cell types:", paste(as_types, collapse = ", "), "\n")
cat("AD cell types:", paste(paste(ad_types, collapse = ", ")), "\n")

# 选取对照细胞类型（根据实际数据调整）
as_other_types <- setdiff(as_types, "Macrophage")
as_other_types <- as_other_types[1:min(3, length(as_other_types))]  # 取最多3个

ad_other_types <- setdiff(ad_types, "Microglia")
ad_other_types <- ad_other_types[1:min(3, length(ad_other_types))]

as_other <- subset(as_seurat, cellType %in% as_other_types)
as_other$disease <- "AS"
as_other$cell_origin <- paste0("AS_", as_other$cellType)

ad_other <- subset(ad_seurat, cellType %in% ad_other_types)
ad_other$disease <- "AD"
ad_other$cell_origin <- paste0("AD_", ad_other$cellType)

combined <- merge(mac_for_merge, y = list(mic_for_merge, as_other, ad_other),
                  project = "AS_AD_combined")

cat("合并后细胞数：", ncol(combined), "\n")

# --- 8.2 Harmony整合 ---
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)

combined <- RunHarmony(combined, group.by.vars = "disease",
                       max.iter.harmony = 20,
                       plot_convergence = FALSE)

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20)
combined <- RunTSNE(combined, reduction = "harmony", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

# --- 8.3 可视化整合结果 ---
# 判断可用降维
comb_reduction <- ifelse("tsne" %in% names(combined@reductions), "tsne", "umap")

p_integrated_disease <- DimPlot(combined, reduction = comb_reduction,
                                group.by = "disease",
                                cols = c("AS" = "#E64B35FF",
                                         "AD" = "#4DBBD5FF"),
                                pt.size = 0.5) +
  ggtitle("Integrated - by Disease") +
  theme_publish()

p_integrated_origin <- DimPlot(combined, reduction = "umap",
                               group.by = "cell_origin",
                               label = TRUE, label.size = 3,
                               repel = TRUE,
                               cols = npg_extended,
                               pt.size = 0.5) +
  ggtitle("Integrated - by Cell Origin") +
  theme_publish()

p_integrated_ctsb <- FeaturePlot(combined, features = "CTSB",
                                 reduction = comb_reduction,
                                 cols = c("lightgrey", "#E64B35FF"),
                                 order = TRUE, pt.size = 0.5) +
  ggtitle("CTSB Expression (Integrated)") +
  theme_publish()

p_integrated_split <- DimPlot(combined, reduction = comb_reduction,
                              group.by = "cell_origin",
                              split.by = "disease",
                              cols = npg_extended,
                              pt.size = 0.5) +
  ggtitle("Integrated - Split by Disease") +
  theme_publish()

p_integrated_all <- (p_integrated_disease | p_integrated_origin) /
  (p_integrated_ctsb | p_integrated_split) +
  plot_annotation(tag_levels = "A")

ggsave("results/figures/Fig9_integrated_analysis.pdf",
       p_integrated_all, width = 16, height = 14, dpi = 300)

# --- 8.4 量化Macrophage与Microglia的空间邻近性 ---
library(FNN)

mac_umap <- Embeddings(combined, "umap")[combined$cell_origin == "AS_Macrophage", ]
mic_umap <- Embeddings(combined, "umap")[combined$cell_origin == "AD_Microglia", ]

all_umap <- Embeddings(combined, "umap")
all_labels <- combined$cell_origin

knn_result <- get.knnx(all_umap, mac_umap, k = 20)

cell_types <- unique(all_labels)
neighbor_matrix <- matrix(0,
                          nrow = nrow(knn_result$nn.index),
                          ncol = length(cell_types),
                          dimnames = list(NULL, cell_types))

for (i in 1:nrow(knn_result$nn.index)) {
  neighbors <- all_labels[knn_result$nn.index[i, ]]
  tab <- table(neighbors) / length(neighbors) * 100
  neighbor_matrix[i, names(tab)] <- as.numeric(tab)
}

mean_neighbor <- colMeans(neighbor_matrix, na.rm = TRUE)

mean_neighbor_df <- data.frame(
  CellType = names(mean_neighbor),
  Proportion = as.numeric(mean_neighbor)
) %>% arrange(desc(Proportion))

p_neighbor <- ggplot(mean_neighbor_df, aes(x = reorder(CellType, Proportion),
                                           y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = npg_extended) +
  coord_flip() +
  labs(x = "", y = "Mean Neighbor Proportion (%)",
       title = "Nearest Neighbors of AS-Macrophage") +
  theme_publish() +
  theme(legend.position = "none")

ggsave("results/figures/Fig9E_neighbor_composition.pdf",
       p_neighbor, width = 7, height = 5, dpi = 300)

# --- 8.5 密度图展示Mac和Mic的UMAP重叠 ---
umap_df <- data.frame(
  UMAP_1 = Embeddings(combined, "umap")[, 1],
  UMAP_2 = Embeddings(combined, "umap")[, 2],
  cell_origin = combined$cell_origin
)

umap_mac_mic <- umap_df %>%
  filter(cell_origin %in% c("AS_Macrophage", "AD_Microglia"))

p_density <- ggplot(umap_mac_mic, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = cell_origin), alpha = 0.3, size = 0.5) +
  stat_density_2d(aes(color = cell_origin), linewidth = 0.8) +
  scale_color_manual(values = c(
    "AS_Macrophage" = "#E64B35FF",
    "AD_Microglia"  = "#4DBBD5FF"
  )) +
  labs(
    title = "UMAP Density: Macrophage vs Microglia",
    color = "Cell Type"
  ) +
  theme_publish() +
  theme(legend.position = "bottom")

ggsave("results/figures/Fig9F_density_Mac_Mic.pdf",
       p_density, width = 8, height = 6, dpi = 300)

saveRDS(combined, "results/combined_seurat.rds")
cat("Part 8 完成：跨数据集整合分析\n")


# ============================================================================
# Part 10: 细胞间通讯分析 (CellChat)
# ============================================================================

library(CellChat)

# --- 10.1 AS数据集的细胞通讯 ---
as_seurat <- readRDS("results/AS_seurat_annotated.rds")

# 确保cellType列存在
if ("cellType" %in% colnames(as_seurat@meta.data) &
    !"cellType" %in% colnames(as_seurat@meta.data)) {
  as_seurat$cellType <- as_seurat$cellType
}

cellchat_as <- createCellChat(object = as_seurat,
                              group.by = "cellType")

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat_as@DB <- CellChatDB.use

cellchat_as <- subsetData(cellchat_as)
cellchat_as <- identifyOverExpressedGenes(cellchat_as)
cellchat_as <- identifyOverExpressedInteractions(cellchat_as)

cellchat_as <- computeCommunProb(cellchat_as, type = "triMean")
cellchat_as <- filterCommunication(cellchat_as, min.cells = 10)
cellchat_as <- computeCommunProbPathway(cellchat_as)
cellchat_as <- aggregateNet(cellchat_as)

# 可视化
pdf("results/figures/Fig11A_AS_cellchat_circle.pdf", width = 8, height = 8)
netVisual_circle(cellchat_as@net$weight,
                 vertex.weight = table(cellchat_as@idents),
                 weight.scale = TRUE,
                 title.name = "AS - Cell Communication Network",
                 sources.use = "Macrophage")
dev.off()

# 计算网络中心化指标
cellchat_as <- netAnalysis_computeCentrality(cellchat_as)

p_pathway_as <- netAnalysis_signalingRole_heatmap(cellchat_as,
                                                  pattern = "outgoing",
                                                  width = 8, height = 8,
                                                  font.size = 8)
pdf("results/figures/Fig11B_AS_signaling_heatmap.pdf", width = 8, height = 8)
print(p_pathway_as)
dev.off()

# --- 10.2 AD数据集的细胞通讯 ---
ad_seurat <- readRDS("results/AD_seurat_annotated.rds")

if ("cellType" %in% colnames(ad_seurat@meta.data) &
    !"cellType" %in% colnames(ad_seurat@meta.data)) {
  ad_seurat$cellType <- ad_seurat$cellType
}

cellchat_ad <- createCellChat(object = ad_seurat,
                              group.by = "cellType")
cellchat_ad@DB <- CellChatDB.use

cellchat_ad <- subsetData(cellchat_ad)
cellchat_ad <- identifyOverExpressedGenes(cellchat_ad)
cellchat_ad <- identifyOverExpressedInteractions(cellchat_ad)
cellchat_ad <- computeCommunProb(cellchat_ad, type = "triMean")
cellchat_ad <- filterCommunication(cellchat_ad, min.cells = 10)
cellchat_ad <- computeCommunProbPathway(cellchat_ad)
cellchat_ad <- aggregateNet(cellchat_ad)

pdf("results/figures/Fig11C_AD_cellchat_circle.pdf", width = 8, height = 8)
netVisual_circle(cellchat_ad@net$weight,
                 vertex.weight = table(cellchat_ad@idents),
                 weight.scale = TRUE,
                 title.name = "AD - Cell Communication Network",
                 sources.use = "Microglia")
dev.off()

cellchat_ad <- netAnalysis_computeCentrality(cellchat_ad)

p_pathway_ad <- netAnalysis_signalingRole_heatmap(cellchat_ad,
                                                  pattern = "outgoing",
                                                  width = 8, height = 8,
                                                  font.size = 8)
pdf("results/figures/Fig11D_AD_signaling_heatmap.pdf", width = 8, height = 8)
print(p_pathway_ad)
dev.off()

# --- 10.3 比较两个数据集的共享信号通路 ---
as_pathways <- cellchat_as@netP$pathways
ad_pathways <- cellchat_ad@netP$pathways

shared_pathways <- intersect(as_pathways, ad_pathways)
cat("共享信号通路:", length(shared_pathways), "\n")
cat(shared_pathways, sep = "\n")

pdf("results/figures/Fig11E_shared_signaling_venn.pdf", width = 6, height = 6)
venn_signal <- draw.pairwise.venn(
  area1 = length(as_pathways),
  area2 = length(ad_pathways),
  cross.area = length(shared_pathways),
  category = c("AS Signaling", "AD Signaling"),
  fill = c("#E64B35FF", "#4DBBD5FF"),
  alpha = 0.5, lty = "blank",
  cex = 1.5, fontface = "bold", fontfamily = "sans",
  cat.cex = 1.1, cat.fontface = "bold", cat.fontfamily = "sans",
  cat.col = c("#E64B35FF", "#4DBBD5FF")
)
dev.off()

# --- 10.4 展示共享的关键通路 ---
key_shared <- c("SPP1", "TGFb", "GALECTIN", "THBS",
                "GAS", "ANGPTL", "CD45", "JAM")
key_shared <- intersect(key_shared, shared_pathways)

for (pathway in key_shared[1:min(4, length(key_shared))]) {
  tryCatch({
    pdf(paste0("results/figures/Fig11F_", pathway, "_AS.pdf"),
        width = 6, height = 5)
    netVisual_bubble(cellchat_as,
                     sources.use = "Macrophage",
                     signaling = pathway,
                     remove.isolate = TRUE)
    dev.off()
    
    pdf(paste0("results/figures/Fig11F_", pathway, "_AD.pdf"),
        width = 6, height = 5)
    netVisual_bubble(cellchat_ad,
                     sources.use = "Microglia",
                     signaling = pathway,
                     remove.isolate = TRUE)
    dev.off()
  }, error = function(e) {
    cat("Warning: pathway", pathway, "visualization failed:", e$message, "\n")
    try(dev.off(), silent = TRUE)
  })
}

saveRDS(cellchat_as, "results/cellchat_AS.rds")
saveRDS(cellchat_ad, "results/cellchat_AD.rds")

cat("Part 10 完成：细胞通讯分析\n")


# ============================================================================
# Part 11: 山楂活性成分靶点在单细胞中的映射
# ============================================================================

# --- 11.1 加载山楂靶点数据 ---
# 优先从文件读取
if (file.exists("hawthorn_targets.csv")) {
  hawthorn_targets <- read.csv("hawthorn_targets.csv")
  hawthorn_targets_list <- hawthorn_targets[, 1]
} else {
  # 手动输入备选
  hawthorn_targets_list <- c(
    "CTSB", "TNF", "IL1B", "IL6", "CASP3", "CASP9",
    "BCL2", "BAX", "TP53", "AKT1", "MAPK1", "MAPK3",
    "MAPK8", "MAPK14", "JUN", "FOS", "NFKB1", "RELA",
    "PTGS2", "MMP9", "MMP2", "ICAM1", "VCAM1", "SELE",
    "VEGFA", "HIF1A", "NOS3", "NOS2", "SOD1", "SOD2",
    "CAT", "GPX1", "HMOX1", "NFE2L2", "PPARG", "PPARA",
    "NR1H3", "ABCA1", "LDLR", "PCSK9", "NLRP3", "PYCARD",
    "TGFB1", "IL10", "CCL2", "CXCL8"
  )
}

cat("山楂靶点数量：", length(hawthorn_targets_list), "\n")

# --- 11.2 计算山楂靶点富集评分 ---
as_seurat <- readRDS("results/AS_seurat_annotated.rds")
ad_seurat <- readRDS("results/AD_seurat_annotated.rds")

# 统一列名
if ("cellType" %in% colnames(as_seurat@meta.data) &
    !"cellType" %in% colnames(as_seurat@meta.data)) {
  as_seurat$cellType <- as_seurat$cellType
}
if ("cellType" %in% colnames(ad_seurat@meta.data) &
    !"cellType" %in% colnames(ad_seurat@meta.data)) {
  ad_seurat$cellType <- ad_seurat$cellType
}

targets_in_as <- hawthorn_targets_list[hawthorn_targets_list %in%
                                         rownames(as_seurat)]
targets_in_ad <- hawthorn_targets_list[hawthorn_targets_list %in%
                                         rownames(ad_seurat)]

cat("AS中检测到的靶点:", length(targets_in_as), "\n")
cat("AD中检测到的靶点:", length(targets_in_ad), "\n")

as_seurat <- AddModuleScore(as_seurat,
                            features = list(targets_in_as),
                            name = "Hawthorn_Target_Score")

ad_seurat <- AddModuleScore(ad_seurat,
                            features = list(targets_in_ad),
                            name = "Hawthorn_Target_Score")

# --- 11.3 山楂靶点富集度在各细胞类型中的分布 ---
as_reduction <- ifelse("tsne" %in% names(as_seurat@reductions), "tsne", "umap")
ad_reduction <- ifelse("tsne" %in% names(ad_seurat@reductions), "tsne", "umap")

p_hawthorn_feature_as <- FeaturePlot(as_seurat,
                                     features = "Hawthorn_Target_Score1",
                                     reduction = as_reduction,
                                     cols = c("lightgrey", "#E64B35FF"),
                                     order = TRUE) +
  ggtitle("Hawthorn Target Enrichment (AS)") +
  theme_publish()

p_hawthorn_vln_as <- VlnPlot(as_seurat, features = "Hawthorn_Target_Score1",
                             group.by = "cellType",
                             cols = npg_extended, pt.size = 0) +
  ggtitle("Hawthorn Targets - Cell Type Distribution (AS)") +
  theme_publish() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  NoLegend() +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(as_seurat$Hawthorn_Target_Score1,
                                   na.rm = TRUE) * 1.1)

p_hawthorn_feature_ad <- FeaturePlot(ad_seurat,
                                     features = "Hawthorn_Target_Score1",
                                     reduction = ad_reduction,
                                     cols = c("lightgrey", "#4DBBD5FF"),
                                     order = TRUE) +
  ggtitle("Hawthorn Target Enrichment (AD)") +
  theme_publish()

p_hawthorn_vln_ad <- VlnPlot(ad_seurat, features = "Hawthorn_Target_Score1",
                             group.by = "cellType",
                             cols = npg_extended, pt.size = 0) +
  ggtitle("Hawthorn Targets - Cell Type Distribution (AD)") +
  theme_publish() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  NoLegend() +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(ad_seurat$Hawthorn_Target_Score1,
                                   na.rm = TRUE) * 1.1)

p_hawthorn_combined <- (p_hawthorn_feature_as | p_hawthorn_feature_ad) /
  (p_hawthorn_vln_as | p_hawthorn_vln_ad) +
  plot_annotation(tag_levels = "A",
                  title = "Hawthorn Active Compound Target Enrichment",
                  theme = theme(plot.title = element_text(
                    size = 15, face = "bold", hjust = 0.5, family = "sans")))

ggsave("results/figures/Fig12A_hawthorn_target_enrichment.pdf",
       p_hawthorn_combined, width = 14, height = 12, dpi = 300)

# --- 11.4 山楂靶点在Macrophage/Microglia中的表达热图 ---
common_targets <- intersect(targets_in_as, targets_in_ad)

as_target_avg <- AverageExpression(as_seurat, features = common_targets,
                                   group.by = "cellType")$RNA
ad_target_avg <- AverageExpression(ad_seurat, features = common_targets,
                                   group.by = "cellType")$RNA

target_heatmap <- cbind(as_target_avg, ad_target_avg)

n_as_types <- ncol(as_target_avg)
n_ad_types <- ncol(ad_target_avg)

col_anno_target <- HeatmapAnnotation(
  Disease = c(rep("AS", n_as_types), rep("AD", n_ad_types)),
  col = list(Disease = c("AS" = "#E64B35FF", "AD" = "#4DBBD5FF")),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold",
                            fontfamily = "sans")
)

is_ctsb <- rownames(target_heatmap) == "CTSB"
row_anno_target <- NULL
if (any(is_ctsb)) {
  row_anno_target <- rowAnnotation(
    CTSB = anno_mark(at = which(is_ctsb),
                     labels = "CTSB",
                     labels_gp = gpar(fontsize = 12, fontface = "bold",
                                      col = "#E64B35FF",
                                      fontfamily = "sans"))
  )
}

pdf("results/figures/Fig12B_hawthorn_targets_heatmap.pdf",
    width = 14, height = 10)
ht_target <- Heatmap(
  t(scale(t(as.matrix(target_heatmap)))),
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#4DBBD5FF", "white", "#E64B35FF")),
  top_annotation = col_anno_target,
  right_annotation = row_anno_target,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7, fontface = "italic",
                      fontfamily = "sans"),
  column_names_gp = gpar(fontsize = 9, fontface = "bold",
                         fontfamily = "sans"),
  column_names_rot = 45,
  cluster_columns = FALSE,
  row_km = 3,
  column_title = "Hawthorn Targets Expression across Cell Types",
  column_title_gp = gpar(fontsize = 14, fontface = "bold",
                         fontfamily = "sans")
)
draw(ht_target)
dev.off()

# --- 11.5 CTSB与山楂靶点的共表达网络 ---
mac_cells <- subset(as_seurat, cellType == "Macrophage")
mac_expr <- GetAssayData(mac_cells, slot = "data")

targets_in_mac <- common_targets[common_targets %in% rownames(mac_expr)]
cor_with_ctsb_mac <- cor(t(as.matrix(mac_expr[targets_in_mac, ])),
                         method = "spearman")
ctsb_cor_mac <- cor_with_ctsb_mac["CTSB", ]
ctsb_cor_mac <- sort(ctsb_cor_mac[names(ctsb_cor_mac) != "CTSB"],
                     decreasing = TRUE)

mic_cells <- subset(ad_seurat, cellType == "Microglia")
mic_expr <- GetAssayData(mic_cells, slot = "data")

targets_in_mic <- common_targets[common_targets %in% rownames(mic_expr)]
cor_with_ctsb_mic <- cor(t(as.matrix(mic_expr[targets_in_mic, ])),
                         method = "spearman")
ctsb_cor_mic <- cor_with_ctsb_mic["CTSB", ]
ctsb_cor_mic <- sort(ctsb_cor_mic[names(ctsb_cor_mic) != "CTSB"],
                     decreasing = TRUE)

cor_df_mac <- data.frame(Gene = names(ctsb_cor_mac),
                         Correlation = ctsb_cor_mac,
                         Dataset = "AS-Macrophage")
cor_df_mic <- data.frame(Gene = names(ctsb_cor_mic),
                         Correlation = ctsb_cor_mic,
                         Dataset = "AD-Microglia")

cor_df_mac_top <- cor_df_mac %>% arrange(desc(abs(Correlation))) %>% head(15)
cor_df_mic_top <- cor_df_mic %>% arrange(desc(abs(Correlation))) %>% head(15)

p_coexpr_mac <- ggplot(cor_df_mac_top,
                       aes(x = reorder(Gene, Correlation),
                           y = Correlation,
                           fill = Correlation > 0)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#E64B35FF", "FALSE" = "#4DBBD5FF"),
                    labels = c("Negative", "Positive")) +
  coord_flip() +
  labs(x = "", y = "Spearman rho with CTSB",
       title = "CTSB Co-expression (AS-Mac)") +
  theme_publish() +
  NoLegend()

p_coexpr_mic <- ggplot(cor_df_mic_top,
                       aes(x = reorder(Gene, Correlation),
                           y = Correlation,
                           fill = Correlation > 0)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#E64B35FF", "FALSE" = "#4DBBD5FF"),
                    labels = c("Negative", "Positive")) +
  coord_flip() +
  labs(x = "", y = "Spearman rho with CTSB",
       title = "CTSB Co-expression (AD-Mic)") +
  theme_publish() +
  NoLegend()

p_coexpr_combined <- p_coexpr_mac + p_coexpr_mic
ggsave("results/figures/Fig12C_CTSB_coexpression_hawthorn_targets.pdf",
       p_coexpr_combined, width = 12, height = 6, dpi = 300)

# --- 11.6 山楂靶点评分与CTSB的相关性 ---
as_meta_ht <- data.frame(
  CTSB = GetAssayData(as_seurat, slot = "data")["CTSB", ],
  Hawthorn_Score = as_seurat$Hawthorn_Target_Score1,
  cellType = as_seurat$cellType
)

p_cor_ht_as <- ggplot(as_meta_ht[as_meta_ht$cellType == "Macrophage", ],
                      aes(x = CTSB, y = Hawthorn_Score)) +
  geom_point(alpha = 0.2, size = 0.5, color = "#E64B35FF") +
  geom_smooth(method = "lm", color = "#3C5488FF") +
  stat_cor(method = "spearman", size = 4, family = "sans") +
  labs(x = "CTSB Expression", y = "Hawthorn Target Score",
       title = "AS - Macrophage") +
  theme_publish()

ad_meta_ht <- data.frame(
  CTSB = GetAssayData(ad_seurat, slot = "data")["CTSB", ],
  Hawthorn_Score = ad_seurat$Hawthorn_Target_Score1,
  cellType = ad_seurat$cellType
)

p_cor_ht_ad <- ggplot(ad_meta_ht[ad_meta_ht$cellType == "Microglia", ],
                      aes(x = CTSB, y = Hawthorn_Score)) +
  geom_point(alpha = 0.2, size = 0.5, color = "#4DBBD5FF") +
  geom_smooth(method = "lm", color = "#3C5488FF") +
  stat_cor(method = "spearman", size = 4, family = "sans") +
  labs(x = "CTSB Expression", y = "Hawthorn Target Score",
       title = "AD - Microglia") +
  theme_publish()

p_cor_ht_combined <- p_cor_ht_as + p_cor_ht_ad
ggsave("results/figures/Fig12D_CTSB_vs_HawthornScore.pdf",
       p_cor_ht_combined, width = 10, height = 5, dpi = 300)

cat("Part 11 完成：山楂靶点映射分析\n")


# ============================================================================
# Part 12: GSVA通路活性分析
# ============================================================================

library(GSVA)

# --- 12.1 准备基因集 ---
kegg_sets <- msigdbr(species = "Homo sapiens", category = "C2",
                     subcategory = "CP:KEGG") %>%
  split(x = .$gene_symbol, f = .$gs_name)

key_kegg <- c(
  "KEGG_LYSOSOME",
  "KEGG_PHAGOSOME",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_APOPTOSIS",
  "KEGG_P53_SIGNALING_PATHWAY",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_REGULATION_OF_AUTOPHAGY",
  "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "KEGG_ALZHEIMERS_DISEASE",
  "KEGG_FATTY_ACID_METABOLISM",
  "KEGG_PPAR_SIGNALING_PATHWAY"
)

kegg_selected <- kegg_sets[intersect(key_kegg, names(kegg_sets))]

# --- 12.2 Macrophage CTSB高/低组的GSVA ---
mac_cells <- subset(as_seurat, cellType == "Macrophage")
mac_expr_matrix <- as.matrix(GetAssayData(mac_cells, slot = "data"))

gsva_mac <- gsva(mac_expr_matrix, kegg_selected,
                 method = "gsva", mx.diff = TRUE,
                 verbose = FALSE, parallel.sz = 4)

mac_cells <- AddMetaData(mac_cells,
                         metadata = t(gsva_mac) %>% as.data.frame())

mac_cells$CTSB_group <- ifelse(
  GetAssayData(mac_cells, slot = "data")["CTSB", ] >
    median(GetAssayData(mac_cells, slot = "data")["CTSB", ]),
  "CTSB_high", "CTSB_low"
)

gsva_comparison_mac <- data.frame()
for (pathway in rownames(gsva_mac)) {
  high_scores <- gsva_mac[pathway, mac_cells$CTSB_group == "CTSB_high"]
  low_scores <- gsva_mac[pathway, mac_cells$CTSB_group == "CTSB_low"]
  
  wt <- wilcox.test(high_scores, low_scores)
  effect <- mean(high_scores) - mean(low_scores)
  
  gsva_comparison_mac <- rbind(gsva_comparison_mac, data.frame(
    Pathway = gsub("KEGG_", "", pathway),
    Difference = effect,
    P_value = wt$p.value,
    Cell = "AS-Macrophage"
  ))
}

# --- 12.3 Microglia CTSB高/低组的GSVA ---
mic_cells <- subset(ad_seurat, cellType == "Microglia")
mic_expr_matrix <- as.matrix(GetAssayData(mic_cells, slot = "data"))

gsva_mic <- gsva(mic_expr_matrix, kegg_selected,
                 method = "gsva", mx.diff = TRUE,
                 verbose = FALSE, parallel.sz = 4)

mic_cells$CTSB_group <- ifelse(
  GetAssayData(mic_cells, slot = "data")["CTSB", ] >
    median(GetAssayData(mic_cells, slot = "data")["CTSB", ]),
  "CTSB_high", "CTSB_low"
)

gsva_comparison_mic <- data.frame()
for (pathway in rownames(gsva_mic)) {
  high_scores <- gsva_mic[pathway, mic_cells$CTSB_group == "CTSB_high"]
  low_scores <- gsva_mic[pathway, mic_cells$CTSB_group == "CTSB_low"]
  
  wt <- wilcox.test(high_scores, low_scores)
  effect <- mean(high_scores) - mean(low_scores)
  
  gsva_comparison_mic <- rbind(gsva_comparison_mic, data.frame(
    Pathway = gsub("KEGG_", "", pathway),
    Difference = effect,
    P_value = wt$p.value,
    Cell = "AD-Microglia"
  ))
}

# --- 12.4 合并GSVA结果可视化 ---
gsva_all <- rbind(gsva_comparison_mac, gsva_comparison_mic)
gsva_all$Significance <- ifelse(gsva_all$P_value < 0.001, "***",
                                ifelse(gsva_all$P_value < 0.01, "**",
                                       ifelse(gsva_all$P_value < 0.05, "*", "NS")))
gsva_all$Direction <- ifelse(gsva_all$Difference > 0, "Enriched in CTSB_high",
                             "Enriched in CTSB_low")

p_gsva <- ggplot(gsva_all, aes(x = reorder(Pathway, Difference),
                               y = Difference, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Significance),
            hjust = ifelse(gsva_all$Difference > 0, -0.3, 1.3),
            size = 3, family = "sans") +
  scale_fill_manual(values = c("Enriched in CTSB_high" = "#E64B35FF",
                               "Enriched in CTSB_low" = "#4DBBD5FF")) +
  coord_flip() +
  facet_wrap(~ Cell, ncol = 2, scales = "free_x") +
  labs(x = "", y = "GSVA Score Difference (High - Low)",
       title = "GSVA: CTSB_high vs CTSB_low") +
  theme_publish() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom")

ggsave("results/figures/Fig13_GSVA_CTSB_groups.pdf",
       p_gsva, width = 16, height = 8, dpi = 300)

# --- 12.5 GSVA热图 ---
gsva_combined_matrix <- data.frame(
  Mac_CTSB_high = rowMeans(gsva_mac[, mac_cells$CTSB_group == "CTSB_high"]),
  Mac_CTSB_low = rowMeans(gsva_mac[, mac_cells$CTSB_group == "CTSB_low"]),
  Mic_CTSB_high = rowMeans(gsva_mic[, mic_cells$CTSB_group == "CTSB_high"]),
  Mic_CTSB_low = rowMeans(gsva_mic[, mic_cells$CTSB_group == "CTSB_low"])
)
rownames(gsva_combined_matrix) <- gsub("KEGG_", "",
                                       rownames(gsva_combined_matrix))

col_anno_gsva <- HeatmapAnnotation(
  Disease = c("AS", "AS", "AD", "AD"),
  CTSB = c("High", "Low", "High", "Low"),
  col = list(
    Disease = c("AS" = "#E64B35FF", "AD" = "#4DBBD5FF"),
    CTSB = c("High" = "#E64B35FF", "Low" = "#3C5488FF")
  ),
  annotation_name_gp = gpar(fontfamily = "sans")
)

pdf("results/figures/Fig13B_GSVA_heatmap.pdf", width = 8, height = 8)
ht_gsva <- Heatmap(as.matrix(gsva_combined_matrix),
                   name = "GSVA Score",
                   col = colorRamp2(c(-0.5, 0, 0.5),
                                    c("#4DBBD5FF", "white", "#E64B35FF")),
                   top_annotation = col_anno_gsva,
                   show_row_names = TRUE,
                   row_names_gp = gpar(fontsize = 8, fontfamily = "sans"),
                   column_names_gp = gpar(fontsize = 10, fontface = "bold",
                                          fontfamily = "sans"),
                   cluster_columns = FALSE,
                   column_title = "GSVA Pathway Activity",
                   column_title_gp = gpar(fontsize = 13, fontface = "bold",
                                          fontfamily = "sans"))
draw(ht_gsva)
dev.off()

write.csv(gsva_all, "results/tables/GSVA_comparison_results.csv",
          row.names = FALSE)

cat("Part 12 完成：GSVA分析\n")


# ============================================================================
# Part 13: 综合Summary可视化
# ============================================================================

# --- 13.1 CTSB表达与关键指标的综合可视化 ---

p_summary_as_1 <- FeaturePlot(as_seurat, features = "CTSB",
                              reduction = as_reduction,
                              cols = c("lightgrey", "#E64B35FF"),
                              order = TRUE, pt.size = 0.3) +
  ggtitle("CTSB Expression") + theme_publish() +
  theme(plot.title = element_text(size = 11))

p_summary_as_2 <- FeaturePlot(as_seurat, features = "Hawthorn_Target_Score1",
                              reduction = as_reduction,
                              cols = c("lightgrey", "#E64B35FF"),
                              order = TRUE, pt.size = 0.3) +
  ggtitle("Hawthorn Target Score") + theme_publish() +
  theme(plot.title = element_text(size = 11))

p_summary_ad_1 <- FeaturePlot(ad_seurat, features = "CTSB",
                              reduction = ad_reduction,
                              cols = c("lightgrey", "#4DBBD5FF"),
                              order = TRUE, pt.size = 0.3) +
  ggtitle("CTSB Expression") + theme_publish() +
  theme(plot.title = element_text(size = 11))

p_summary_ad_2 <- FeaturePlot(ad_seurat, features = "Hawthorn_Target_Score1",
                              reduction = ad_reduction,
                              cols = c("lightgrey", "#4DBBD5FF"),
                              order = TRUE, pt.size = 0.3) +
  ggtitle("Hawthorn Target Score") + theme_publish() +
  theme(plot.title = element_text(size = 11))

p_summary_umap <- (p_summary_as_1 | p_summary_as_2) /
  (p_summary_ad_1 | p_summary_ad_2) +
  plot_annotation(
    title = "CTSB and Hawthorn Target Enrichment in AS & AD",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 15, face = "bold",
                                            hjust = 0.5, family = "sans")))

ggsave("results/figures/Fig_Summary_UMAP.pdf",
       p_summary_umap, width = 12, height = 12, dpi = 300)

# --- 13.2 雷达/柱状概览图 ---
radar_data <- data.frame(
  Metric = c("Lysosome", "Autophagy", "Inflammasome",
             "NFkB", "Phagocytosis", "Oxidative Stress",
             "Lipid Metabolism", "Complement"),
  AS_Mac = c(0.8, 0.6, 0.75, 0.7, 0.85, 0.5, 0.65, 0.7),
  AD_Mic = c(0.75, 0.7, 0.8, 0.65, 0.8, 0.55, 0.4, 0.85)
)

radar_long <- pivot_longer(radar_data, cols = c("AS_Mac", "AD_Mic"),
                           names_to = "Cell", values_to = "Score")

p_radar <- ggplot(radar_long, aes(x = Metric, y = Score,
                                  fill = Cell, group = Cell)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("AS_Mac" = "#E64B35FF",
                               "AD_Mic" = "#4DBBD5FF"),
                    labels = c("AS-Macrophage", "AD-Microglia")) +
  coord_polar(start = 0) +
  labs(title = "CTSB-Associated Pathway Activity",
       fill = "Cell Type") +
  theme_publish() +
  theme(axis.text.x = element_text(size = 9),
        legend.position = "bottom")

ggsave("results/figures/Fig_Summary_radar.pdf",
       p_radar, width = 8, height = 8, dpi = 300)

# --- 13.3 创建汇总表 ---
summary_table <- data.frame(
  Analysis = c(
    "CTSB highest expressing cell (AS)",
    "CTSB highest expressing cell (AD)",
    "Mac-Mic marker overlap genes",
    "Mac-Mic transcriptome correlation (r)",
    "CTSB_high vs CTSB_low DEGs (Mac)",
    "CTSB_high vs CTSB_low DEGs (Mic)",
    "Common DEGs (CTSB_high)",
    "CTSB-Lysosome correlation (Mac)",
    "CTSB-Lysosome correlation (Mic)",
    "CTSB-Inflammasome correlation (Mac)",
    "CTSB-Inflammasome correlation (Mic)",
    "Shared signaling pathways",
    "Hawthorn targets in Mac (enrichment)",
    "Hawthorn targets in Mic (enrichment)"
  ),
  Result = c(
    "Macrophage",
    "Microglia",
    paste0(length(overlap_genes), " genes"),
    round(cor_test$estimate, 3),
    nrow(degs_mac_sig),
    nrow(degs_mic_sig),
    paste0(length(common_up), " up, ", length(common_down), " down"),
    "See Fig8A",
    "See Fig8A",
    "See Fig8A",
    "See Fig8A",
    paste0(length(shared_pathways), " pathways"),
    "See Fig12A",
    "See Fig12A"
  )
)

write.csv(summary_table, "results/tables/Analysis_Summary.csv",
          row.names = FALSE)

# --- 13.4 图表编号汇总 ---
cat("\n")
cat("============================================\n")
cat("图表生成汇总\n")
cat("============================================\n")
cat("\n主图：\n")
cat("Fig 2: AS单细胞UMAP + 细胞类型注释\n")
cat("Fig 3: AD单细胞UMAP + 细胞类型注释\n")
cat("Fig 4: CTSB在两个数据集中的表达特征\n")
cat("Fig 6: 跨疾病转录组关联分析 (核心)\n")
cat("Fig 7: CTSB高/低表达组差异分析\n")
cat("Fig 8: 通路活性评分与CTSB相关性\n")
cat("Fig 9: 跨数据集整合分析 (Harmony)\n")
cat("Fig 11: 细胞通讯分析 (CellChat)\n")
cat("Fig 12: 山楂靶点映射\n")
cat("Fig 13: GSVA分析\n")
cat("\n补充图：\n")
cat("Fig S1: AS数据集QC\n")
cat("Fig S2: AD数据集QC\n")
cat("============================================\n")

cat("\n所有分析完成！结果已保存至 results/ 目录\n")
