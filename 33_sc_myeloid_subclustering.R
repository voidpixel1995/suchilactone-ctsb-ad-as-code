# =============================================================================
# Step 3.4 - Sub-clustering of macrophages (AS) and microglia (AD)
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

setwd(file.path(WORK_DIR, "singlecell"))   # set by 00_config.R

# ====================================================
# 3.1 AS巨噬细胞亚聚类
# ====================================================

macro <- subset(AS_obj, cellType == "Macrophage")
cat("AS巨噬细胞数量：", ncol(macro), "\n")

macro <- NormalizeData(macro) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(dims = 1:15) %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = 0.5)%>%
  RunUMAP(dims = 1:15)%>%
  RunTSNE(dims = 1:15)

##查找每个聚类的差异基因
logFCfilter=1           #logFC的过滤条件
adjPvalFilter=0.05     #矫正后的pvalue的过滤条件

pbmc.markers <- FindAllMarkers(object = macro,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.table(top10,file="03.top10clusterMarkers.txt",sep="\t",row.names=F,quote=F)


# 巨噬细胞亚型Marker
macro_markers <- c(
  # M1促炎
  "TNF", "IL1B", "IL6", "NOS2", "CXCL10", "CCL2",
  # M2抗炎
  "MRC1", "CD163", "IL10", "TGFB1", "FOLR2",
  # 泡沫细胞
  "ABCA1", "CD36", "OLR1", "MSR1",
  # TREM2+
  "TREM2", "SPP1", "GPNMB", "FABP5", "LPL",
  # 居住型
  "LYVE1", "F13A1", "TIMD4",
  # CTSB和炎症小体
  "CTSB", "NLRP3", "CASP1", "GSDMD", "IL18", "PYCARD"
)
macro_markers <- macro_markers[macro_markers %in% rownames(macro)]

p_macro_dot <- DotPlot(macro, features = macro_markers,
                       cols = c("lightgrey", "#B22222")) +
  RotatedAxis() + ggtitle("Macrophage subtype markers + CTSB/inflammasome")
ggsave("results/03_AS_macro_marker_dotplot.pdf", p_macro_dot, width = 20, height = 6)

# 手动注释亚型（★根据你的结果修改★）
macro_sub_map <- c(
  "0" = "M1_proinflammatory",
  "1" = "TREM2+_Mac",
  "2" = "M2_antiinflammatory",
  "3" = "Foam_cell",
  "4" = "Resident_Mac"
)
macro$macro_subtype <- macro_sub_map[as.character(Idents(macro))]

# 亚型UMAP
p_macro_sub <- DimPlot(macro, group.by = "macro_subtype", label = TRUE, repel = TRUE,
                       cols = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")) +
  ggtitle("Macrophage subtypes (AS)")
ggsave("results/03_AS_macro_subtype_UMAP.pdf", p_macro_sub, width = 8, height = 7)

# CTSB在亚群中的表达
p_ctsb_macro_feature <- FeaturePlot(macro, features = "CTSB",
                                    cols = c("lightgrey", "#B22222"),
                                    order = TRUE) +
  ggtitle("CTSB in macrophage subtypes")

p_ctsb_macro_vln <- VlnPlot(macro, features = "CTSB", group.by = "macro_subtype",
                            pt.size = 0,
                            cols = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")) +
  ggtitle("CTSB across macrophage subtypes") +
  stat_compare_means(method = "kruskal.test", label = "p.format")

p_ctsb_macro_combined <- p_ctsb_macro_feature | p_ctsb_macro_vln
ggsave("results/03_AS_CTSB_macro_subtype.pdf", p_ctsb_macro_combined, width = 14, height = 6)

# CTSB + 炎症小体panel在亚群中
inflam_genes <- c("CTSB", "CTSD", "CTSL", "NLRP3", "PYCARD", "CASP1",
                  "IL1B", "IL18", "GSDMD", "NEK7", "TXNIP")
inflam_genes <- inflam_genes[inflam_genes %in% rownames(macro)]

p_inflam_macro <- DotPlot(macro, features = inflam_genes, group.by = "macro_subtype",
                          cols = c("lightgrey", "#B22222")) +
  RotatedAxis() + ggtitle("CTSB & NLRP3 inflammasome in macrophage subtypes")
ggsave("results/03_AS_inflammasome_macro.pdf", p_inflam_macro, width = 10, height = 5)

# 亚群比例在疾病vs对照中的变化
macro_prop <- macro@meta.data %>%
  group_by(condition, macro_subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(pct = n / sum(n) * 100)

p_macro_prop <- ggplot(macro_prop, aes(x = condition, y = pct, fill = macro_subtype)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")) +
  theme_classic() + labs(title = "Macrophage subtype proportions", y = "%", x = "")
ggsave("results/03_AS_macro_subtype_proportion.pdf", p_macro_prop, width = 7, height = 5)

saveRDS(macro, "results/AS_macrophage_subset.rds")

# ====================================================
# 3.2 AD小胶质细胞亚聚类
# ====================================================
m
micro <- subset(AD_obj, cellType == "Microglia")
cat("AD小胶质细胞数量：", ncol(micro), "\n")

micro <- NormalizeData(micro) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(dims = 1:15) %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = 0.5)%>%
  RunUMAP(dims = 1:15)%>%
  RunTSNE(dims = 1:15)

# 小胶质亚型Marker
micro_markers <- c(
  # 稳态型
  "P2RY12", "CX3CR1", "TMEM119", "CSF1R", "HEXB",
  # DAM
  "TREM2", "SPP1", "APOE", "LPL", "CST7", "GPNMB",
  # 促炎
  "TNF", "IL1B", "CCL2", "CCL3", "CCL4", "CXCL10",
  # IFN响应
  "IFIT1", "ISG15", "MX1",
  # CTSB和炎症小体
  "CTSB", "CTSD", "NLRP3", "CASP1", "GSDMD", "IL18", "PYCARD"
)
micro_markers <- micro_markers[micro_markers %in% rownames(micro)]

p_micro_dot <- DotPlot(micro, features = micro_markers,
                       cols = c("lightgrey", "#8B008B")) +
  RotatedAxis() + ggtitle("Microglia subtype markers + CTSB/inflammasome")
ggsave("results/03_AD_micro_marker_dotplot.pdf", p_micro_dot, width = 18, height = 6)

# 注释（★根据你的结果修改★）
micro_sub_map <- c(
  "0" = "Homeostatic",
  "1" = "DAM",
  "2" = "Proinflammatory",
  "3" = "IFN_responsive",
  "4" = "Proliferating"
)
micro$micro_subtype <- micro_sub_map[as.character(Idents(micro))]

p_micro_sub <- DimPlot(micro, group.by = "micro_subtype", label = TRUE, repel = TRUE,
                       cols = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E")) +
  ggtitle("Microglia subtypes (AD)")
ggsave("results/03_AD_micro_subtype_UMAP.pdf", p_micro_sub, width = 8, height = 7)

# CTSB在小胶质亚群中
p_ctsb_micro_vln <- VlnPlot(micro, features = "CTSB", group.by = "micro_subtype",
                            pt.size = 0,
                            cols = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E")) +
  ggtitle("CTSB across microglia subtypes") +
  stat_compare_means(method = "kruskal.test", label = "p.format")
ggsave("results/03_AD_CTSB_micro_subtype.pdf", p_ctsb_micro_vln, width = 8, height = 5)

# 炎症小体panel
inflam_genes_micro <- inflam_genes[inflam_genes %in% rownames(micro)]
p_inflam_micro <- DotPlot(micro, features = inflam_genes_micro, group.by = "micro_subtype",
                          cols = c("lightgrey", "#8B008B")) +
  RotatedAxis() + ggtitle("CTSB & NLRP3 inflammasome in microglia subtypes")
ggsave("results/03_AD_inflammasome_micro.pdf", p_inflam_micro, width = 10, height = 5)

saveRDS(micro, "results/AD_microglia_subset.rds")
