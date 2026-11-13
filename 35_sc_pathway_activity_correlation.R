# =============================================================================
# Step 3.6 - Pathway activity scoring and correlation with CTSB expression
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
# 5.1 定义通路基因集
# ====================================================

gene_sets <- list(
  NLRP3_inflammasome = c("NLRP3","PYCARD","CASP1","IL1B","IL18","GSDMD","NEK7","TXNIP"),
  Pyroptosis = c("GSDMD","GSDME","CASP1","CASP4","CASP5","NLRP3","IL1B","IL18"),
  Lysosome = c("CTSB","CTSD","CTSL","CTSS","CTSZ","LAMP1","LAMP2","TFEB","ATP6V1A"),
  Oxidative_stress = c("SOD1","SOD2","CAT","GPX1","GPX4","NFE2L2","NQO1","HMOX1"),
  NFkB = c("NFKB1","NFKB2","RELA","RELB","NFKBIA","IKBKB","CHUK"),
  Apoptosis = c("BAX","BAK1","BCL2","BCL2L1","CASP3","CASP8","CASP9","CYCS"),
  Lipid_metabolism = c("APOE","ABCA1","ABCG1","LPL","LDLR","CD36","FABP4","SREBF1")
)

# ====================================================
# 5.2 AS数据集评分
# ====================================================

gene_sets_AS <- lapply(gene_sets, function(g) g[g %in% rownames(AS_obj)])
gene_sets_AS <- gene_sets_AS[sapply(gene_sets_AS, length) >= 3]

AS_obj <- AddModuleScore(AS_obj, features = gene_sets_AS,
                         name = paste0(names(gene_sets_AS), "_"))

# 重命名列
old_names <- grep("_\\d+$", colnames(AS_obj@meta.data), value = TRUE)
new_names <- sub("_\\d+$", "_score", old_names)
colnames(AS_obj@meta.data)[match(old_names, colnames(AS_obj@meta.data))] <- new_names

# ====================================================
# 5.3 AD数据集评分
# ====================================================

gene_sets_AD <- lapply(gene_sets, function(g) g[g %in% rownames(AD_obj)])
gene_sets_AD <- gene_sets_AD[sapply(gene_sets_AD, length) >= 3]

AD_obj <- AddModuleScore(AD_obj, features = gene_sets_AD,
                         name = paste0(names(gene_sets_AD), "_"))

old_names_ad <- grep("_\\d+$", colnames(AD_obj@meta.data), value = TRUE)
new_names_ad <- sub("_\\d+$", "_score", old_names_ad)
colnames(AD_obj@meta.data)[match(old_names_ad, colnames(AD_obj@meta.data))] <- new_names_ad

# ====================================================
# 5.4 NLRP3通路评分可视化
# ====================================================

# FeaturePlot
p_nlrp3_as <- FeaturePlot(AS_obj, features = "NLRP3_inflammasome_score",
                          cols = c("lightgrey","darkred"), order = TRUE) +
  ggtitle("NLRP3 inflammasome score (AS)")

p_nlrp3_ad <- FeaturePlot(AD_obj,reduction = "tsne", features = "NLRP3_inflammasome_score",
                          cols = c("lightgrey","darkred"), order = TRUE) +
  ggtitle("NLRP3 inflammasome score (AD)")

ggsave("results/05_NLRP3_score_FeaturePlot.pdf", p_nlrp3_as | p_nlrp3_ad, width = 14, height = 6)

# VlnPlot按疾病分组
p_nlrp3_vln_as <- VlnPlot(AS_obj, features = "NLRP3_inflammasome_score",
                          group.by = "cellType", split.by = "Sample",
                          pt.size = 0, cols = c("#4393C3","#D6604D")) +
  ggtitle("NLRP3 score (AS)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_nlrp3_vln_ad <- VlnPlot(AD_obj, features = "NLRP3_inflammasome_score",
                          group.by = "cellType", split.by = "Sample",
                          pt.size = 0, cols = c("#4393C3","#8B008B")) +
  ggtitle("NLRP3 score (AD)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/05_NLRP3_score_VlnPlot.pdf", p_nlrp3_vln_as / p_nlrp3_vln_ad, width = 12, height = 10)

# ====================================================
# 5.5 ★核心图★ CTSB表达与各通路评分的相关性
# ====================================================

# 在AS巨噬细胞中
macro_data <- FetchData(macro, vars = "CTSB")
score_names_AS <- grep("_score$", colnames(AS_obj@meta.data), value = TRUE)
for(s in score_names_AS) {
  macro_data[[s]] <- AS_obj@meta.data[colnames(macro), s]
}

# 计算相关系数
cor_results_macro <- data.frame()
for(s in score_names_AS) {
  ct <- cor.test(macro_data$CTSB, macro_data[[s]], method = "spearman")
  cor_results_macro <- rbind(cor_results_macro, data.frame(
    Pathway = gsub("_score$", "", s),
    rho = ct$estimate,
    pvalue = ct$p.value,
    Dataset = "AS_Macrophage"
  ))
}

# 在AD小胶质细胞中
micro_data <- FetchData(micro, vars = "CTSB")
score_names_AD <- grep("_score$", colnames(AD_obj@meta.data), value = TRUE)
for(s in score_names_AD) {
  micro_data[[s]] <- AD_obj@meta.data[colnames(micro), s]
}

cor_results_micro <- data.frame()
for(s in score_names_AD) {
  ct <- cor.test(micro_data$CTSB, micro_data[[s]], method = "spearman")
  cor_results_micro <- rbind(cor_results_micro, data.frame(
    Pathway = gsub("_score$", "", s),
    rho = ct$estimate,
    pvalue = ct$p.value,
    Dataset = "AD_Microglia"
  ))
}

# 合并
cor_all <- rbind(cor_results_macro, cor_results_micro)
cor_all$sig <- ifelse(cor_all$pvalue < 0.001, "***",
                      ifelse(cor_all$pvalue < 0.01, "**",
                             ifelse(cor_all$pvalue < 0.05, "*", "ns")))

# 相关性柱状图（分面展示）
p_cor_bar <- ggplot(cor_all, aes(x = reorder(Pathway, rho), y = rho, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = sig), position = position_dodge(width = 0.7), 
            vjust = ifelse(cor_all$rho > 0, -0.5, 1.5), size = 4) +
  scale_fill_manual(values = c("AS_Macrophage" = "#D6604D", "AD_Microglia" = "#8B008B")) +
  coord_flip() + theme_classic() +
  labs(title = "CTSB expression correlation with pathway scores",
       x = "", y = "Spearman ρ") +
  geom_hline(yintercept = 0, linetype = "dashed")
ggsave("results/05_CTSB_pathway_correlation_bar.pdf", p_cor_bar, width = 10, height = 6)

# ====================================================
# 5.6 关键散点图：CTSB vs NLRP3通路评分
# ====================================================

p_cor_as <- ggplot(macro_data, aes(x = CTSB, y = NLRP3_inflammasome_score)) +
  geom_point(alpha = 0.15, size = 0.5, color = "#D6604D") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman", size = 4.5) +
  theme_classic() +
  labs(title = "AS Macrophage", x = "CTSB expression", y = "NLRP3 inflammasome score")

p_cor_ad <- ggplot(micro_data, aes(x = CTSB, y = NLRP3_inflammasome_score)) +
  geom_point(alpha = 0.15, size = 0.5, color = "#8B008B") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman", size = 4.5) +
  theme_classic() +
  labs(title = "AD Microglia", x = "CTSB expression", y = "NLRP3 inflammasome score")

ggsave("results/05_CTSB_vs_NLRP3_scatter.pdf", p_cor_as | p_cor_ad, width = 12, height = 5)

# ====================================================
# 5.7 通路评分热图（各细胞类型平均值）
# ====================================================

# AS
pathway_avg_AS <- AS_obj@meta.data %>%
  group_by(cellType) %>%
  summarise(across(all_of(score_names_AS), mean, na.rm = TRUE)) %>%
  column_to_rownames("cellType")
colnames(pathway_avg_AS) <- gsub("_score$", "", colnames(pathway_avg_AS))

pdf("results/05_AS_pathway_heatmap.pdf", width = 9, height = 6)
pheatmap(t(pathway_avg_AS), scale = "row",
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         main = "Pathway activity - AS dataset")
dev.off()

# AD
pathway_avg_AD <- AD_obj@meta.data %>%
  group_by(cellType) %>%
  summarise(across(all_of(score_names_AD), mean, na.rm = TRUE)) %>%
  column_to_rownames("cellType")
colnames(pathway_avg_AD) <- gsub("_score$", "", colnames(pathway_avg_AD))

pdf("results/05_AD_pathway_heatmap.pdf", width = 9, height = 6)
pheatmap(t(pathway_avg_AD), scale = "row",
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         main = "Pathway activity - AD dataset")
dev.off()
