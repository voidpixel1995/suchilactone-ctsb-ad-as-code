# =============================================================================
# Step 3.3 - Cell-type-resolved expression of CTSB and P2RX7
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

# 加载R包
# ====================================================
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(harmony)
library(pheatmap)
library(ggrepel)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(VennDiagram)
library(CellChat)
library(monocle3)
library(SeuratWrappers)
library(scTenifoldKnk)
library(RColorBrewer)
library(viridis)
library(corrplot)
library(RColorBrewer)

set.seed(12345)
options(future.globals.maxSize = 8000 * 1024^2)

setwd(file.path(WORK_DIR, "singlecell"))   # set by 00_config.R
dir.create("results", showWarnings = FALSE)
load("AD_obj.Rdata")
load("AS_obj.Rdata")

# ====================================================
# 2.1 CTSB在两个数据集中的UMAP表达分布
# ====================================================

p_feature_as <- FeaturePlot(AS_obj, features = "CTSB",
                            cols = c("lightgrey", "#B22222"),
                            order = TRUE, pt.size = 0.3) +
  ggtitle("CTSB expression - AS dataset")

p_feature_ad <- FeaturePlot(AD_obj,reduction = "tsne", features = "CTSB",
                            cols = c("lightgrey", "#B22222"),
                            order = TRUE, pt.size = 0.3) +
  ggtitle("CTSB expression - AD dataset")

p_feature_both <- p_feature_as + p_feature_ad
ggsave("results/02_CTSB_FeaturePlot_both.pdf", p_feature_both, width = 10, height = 8)

# ====================================================
# 2.2 CTSB在各细胞类型中的VlnPlot
# ====================================================

p_vln_as <- VlnPlot(AS_obj, features = "CTSB", group.by = "cellType",
                    pt.size = 0, cols = brewer.pal(11, "Set3")) +
  ggtitle("CTSB - AS dataset") + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

p_vln_ad <- VlnPlot(AD_obj, features = "CTSB", group.by = "cellType",
                    pt.size = 0, cols = brewer.pal(9, "Paired")) +
  ggtitle("CTSB - AD dataset") + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

p_vln_both <- p_vln_as / p_vln_ad
ggsave("results/02_CTSB_VlnPlot_both.pdf", p_vln_both, width = 10, height = 10)

# ====================================================
# 2.3 CTSB在各细胞类型中的DotPlot（两数据集并排）
# ====================================================

p_dot_as <- DotPlot(AS_obj, features = "CTSB", group.by = "cellType",
                    cols = c("lightgrey", "#B22222")) +
  coord_flip() + ggtitle("AS") +
  theme(axis.text.y = element_text(size = 11))

p_dot_ad <- DotPlot(AD_obj, features = "CTSB", group.by = "cellType",
                    cols = c("lightgrey", "#B22222")) +
  coord_flip() + ggtitle("AD") +
  theme(axis.text.y = element_text(size = 11))

p_dot_both <- p_dot_as | p_dot_ad
ggsave("results/02_CTSB_DotPlot_both.pdf", p_dot_both, width = 12, height = 5)

# ====================================================
# 2.4 CTSB疾病组 vs 对照组差异（分细胞类型）
# ====================================================

# AS数据集
metadata <- AS_obj@meta.data
unique(metadata$Sample)
p_split_as <- VlnPlot(AS_obj, features = "CTSB",
                      group.by = "cellType", split.by = "Sample",
                      pt.size = 0,
                      cols = c("Control" ="#D6604D" , "Atherosclerosis" ="#4393C3" )) +
  ggtitle("CTSB: AS vs Control") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/02_AS_CTSB_split_condition.pdf", p_split_as, width = 12, height = 5)

# AD数据集
p_split_ad <- VlnPlot(AD_obj, features = "CTSB",
                      group.by = "cellType", split.by = "Sample",
                      pt.size = 0,
                      cols = c("Control" = "#4393C3", "AD" = "#8B008B")) +
  ggtitle("CTSB: AD vs Control") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/02_AD_CTSB_split_condition.pdf", p_split_ad, width = 12, height = 5)

# ====================================================
# 2.5 统计检验：CTSB在关键细胞类型中的差异
# ====================================================

# 函数：计算每种细胞类型中CTSB的差异
calc_CTSB_diff <- function(obj, disease_label, control_label = "Control") {
  results <- data.frame()
  for(ct in unique(obj$cellType)) {
    sub <- subset(obj, cellType == ct)
    expr_disease <- FetchData(subset(sub, Sample == disease_label), "CTSB")$CTSB
    expr_control <- FetchData(subset(sub, Sample == control_label), "CTSB")$CTSB
    
    if(length(expr_disease) > 10 & length(expr_control) > 10) {
      test <- wilcox.test(expr_disease, expr_control)
      results <- rbind(results, data.frame(
        celltype = ct,
        mean_disease = mean(expr_disease),
        mean_control = mean(expr_control),
        log2FC = log2((mean(expr_disease) + 0.01) / (mean(expr_control) + 0.01)),
        p_value = test$p.value,
        n_disease = length(expr_disease),
        n_control = length(expr_control)
      ))
    }
  }
  results$p_adj <- p.adjust(results$p_value, method = "BH")
  results$sig <- ifelse(results$p_adj < 0.05, "*", "ns")
  return(results)
}

# AS数据集
diff_AS <- calc_CTSB_diff(AS_obj, disease_label = "Atherosclerosis")
cat("\n===== CTSB差异（AS vs Control） =====\n")
print(diff_AS)
write.csv(diff_AS, "results/02_CTSB_diff_AS.csv", row.names = FALSE)

# AD数据集
diff_AD <- calc_CTSB_diff(AD_obj, disease_label = "AD")
cat("\n===== CTSB差异（AD vs Control） =====\n")
print(diff_AD)
write.csv(diff_AD, "results/02_CTSB_diff_AD.csv", row.names = FALSE)

# ====================================================
# 2.6 log2FC柱状图
# ====================================================

p_fc_as <- ggplot(diff_AS, aes(x = reorder(celltype, log2FC), y = log2FC)) +
  geom_bar(stat = "identity", aes(fill = log2FC > 0), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#D6604D", "FALSE" = "#4393C3"), guide = "none") +
  geom_text(aes(label = sig), hjust = ifelse(diff_AS$log2FC > 0, -0.3, 1.3), size = 5) +
  coord_flip() + theme_classic() +
  labs(title = "CTSB log2FC (AS vs Control)", x = "", y = "log2FC") +
  geom_hline(yintercept = 0, linetype = "dashed")

p_fc_ad <- ggplot(diff_AD, aes(x = reorder(celltype, log2FC), y = log2FC)) +
  geom_bar(stat = "identity", aes(fill = log2FC > 0), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#8B008B", "FALSE" = "#4393C3"), guide = "none") +
  geom_text(aes(label = sig), hjust = ifelse(diff_AD$log2FC > 0, -0.3, 1.3), size = 5) +
  coord_flip() + theme_classic() +
  labs(title = "CTSB log2FC (AD vs Control)", x = "", y = "log2FC") +
  geom_hline(yintercept = 0, linetype = "dashed")

p_fc_combined <- p_fc_as | p_fc_ad
ggsave("results/02_CTSB_log2FC_barplot.pdf", p_fc_combined, width = 14, height = 6)

# ====================================================
# 2.7 CTSB表达百分比柱状图
# ====================================================

# 计算各细胞类型中CTSB阳性细胞比例
calc_pct <- function(obj, dataset_name) {
  meta <- obj@meta.data
  meta$CTSB_expr <- FetchData(obj, "CTSB")$CTSB
  meta$CTSB_pos <- meta$CTSB_expr > 0
  
  pct_df <- meta %>%
    group_by(cellType, condition) %>%
    summarise(pct_positive = mean(CTSB_pos) * 100,
              mean_expr = mean(CTSB_expr),
              .groups = "drop") %>%
    mutate(dataset = dataset_name)
  return(pct_df)
}

pct_AS <- calc_pct(AS_obj, "AS")
pct_AD <- calc_pct(AD_obj, "AD")

p_pct_as <- ggplot(pct_AS, aes(x = celltype, y = pct_positive, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Control" = "#4393C3", "AS" = "#D6604D")) +
  theme_classic() +
  labs(title = "CTSB⁺ cell percentage (AS)", y = "% CTSB⁺ cells", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_pct_ad <- ggplot(pct_AD, aes(x = celltype, y = pct_positive, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Control" = "#4393C3", "AD" = "#8B008B")) +
  theme_classic() +
  labs(title = "CTSB⁺ cell percentage (AD)", y = "% CTSB⁺ cells", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_pct_both <- p_pct_as / p_pct_ad
ggsave("results/02_CTSB_positive_percentage.pdf", p_pct_both, width = 10, height = 10)
