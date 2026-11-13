# =============================================================================
# Step 3.5 - CTSB-high versus CTSB-low differential expression and enrichment
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
# 4.1 AS巨噬细胞：CTSB_high vs CTSB_low 差异分析
# ====================================================

# 按CTSB中位数分组
ctsb_expr_macro <- FetchData(macro, "CTSB")$CTSB
macro$CTSB_group <- ifelse(ctsb_expr_macro > median(ctsb_expr_macro),
                           "CTSB_high", "CTSB_low")

table(macro$CTSB_group)

# 差异基因分析
Idents(macro) <- "CTSB_group"
DEG_CTSB_macro <- FindMarkers(macro, 
                              ident.1 = "CTSB_high", 
                              ident.2 = "CTSB_low",
                              min.pct = 0.1, 
                              logfc.threshold = 0.25)
DEG_CTSB_macro$gene <- rownames(DEG_CTSB_macro)

cat("CTSB_high上调基因：", sum(DEG_CTSB_macro$avg_log2FC > 0.25 & DEG_CTSB_macro$p_val_adj < 0.05), "\n")
cat("CTSB_high下调基因：", sum(DEG_CTSB_macro$avg_log2FC < -0.25 & DEG_CTSB_macro$p_val_adj < 0.05), "\n")

write.csv(DEG_CTSB_macro, "results/04_AS_macro_CTSB_high_vs_low_DEGs.csv", row.names = FALSE)

# 火山图
DEG_CTSB_macro$label <- ifelse(DEG_CTSB_macro$p_val_adj < 0.05 & 
                                 abs(DEG_CTSB_macro$avg_log2FC) > 0.5,
                               DEG_CTSB_macro$gene, "")
DEG_CTSB_macro$color <- case_when(
  DEG_CTSB_macro$avg_log2FC > 0.25 & DEG_CTSB_macro$p_val_adj < 0.05 ~ "Up",
  DEG_CTSB_macro$avg_log2FC < -0.25 & DEG_CTSB_macro$p_val_adj < 0.05 ~ "Down",
  TRUE ~ "NS"
)

p_volcano_macro <- ggplot(DEG_CTSB_macro, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = color), size = 1, alpha = 0.6) +
  scale_color_manual(values = c("Up" = "#D6604D", "Down" = "#4393C3", "NS" = "grey80")) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 25) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  labs(title = "CTSB_high vs CTSB_low in AS macrophages",
       x = "log2FC", y = "-log10(padj)")
ggsave("results/04_AS_CTSB_high_low_volcano.pdf", p_volcano_macro, width = 9, height = 7)

# ====================================================
# 4.2 GO/KEGG富集分析（AS巨噬细胞CTSB_high上调基因）
# ====================================================

up_genes_macro <- DEG_CTSB_macro$gene[DEG_CTSB_macro$avg_log2FC > 0.25 & 
                                        DEG_CTSB_macro$p_val_adj < 0.05]

# GO-BP
ego_macro <- enrichGO(gene = up_genes_macro,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

p_go_macro <- dotplot(ego_macro, showCategory = 15,
                      title = "GO-BP: CTSB_high upregulated genes (AS Macrophage)")
ggsave("results/04_AS_CTSB_high_GO_BP.pdf", p_go_macro, width = 10, height = 8)

# KEGG
gene_entrez_macro <- bitr(up_genes_macro, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
options(timeout = 6000000)  # 时间设置长一点
ekegg_macro <- enrichKEGG(gene = gene_entrez_macro$ENTREZID,
                          organism = "hsa", pvalueCutoff = 0.05)

if(nrow(ekegg_macro@result[ekegg_macro@result$p.adjust < 0.05, ]) > 0) {
  p_kegg_macro <- dotplot(ekegg_macro, showCategory = 15,
                          title = "KEGG: CTSB_high upregulated genes (AS Macrophage)")
  ggsave("results/04_AS_CTSB_high_KEGG.pdf", p_kegg_macro, width = 10, height = 8)
}

write.csv(ego_macro@result, "results/04_AS_CTSB_high_GO_results.csv", row.names = FALSE)

# ====================================================
# 4.3 AD小胶质细胞：CTSB_high vs CTSB_low
# ====================================================

ctsb_expr_micro <- FetchData(micro, "CTSB")$CTSB
micro$CTSB_group <- ifelse(ctsb_expr_micro > median(ctsb_expr_micro),
                           "CTSB_high", "CTSB_low")

Idents(micro) <- "CTSB_group"
DEG_CTSB_micro <- FindMarkers(micro,
                              ident.1 = "CTSB_high",
                              ident.2 = "CTSB_low",
                              min.pct = 0.1,
                              logfc.threshold = 0.25)
DEG_CTSB_micro$gene <- rownames(DEG_CTSB_micro)

write.csv(DEG_CTSB_micro, "results/04_AD_micro_CTSB_high_vs_low_DEGs.csv", row.names = FALSE)

# 火山图
DEG_CTSB_micro$label <- ifelse(DEG_CTSB_micro$p_val_adj < 0.05 & 
                                 abs(DEG_CTSB_micro$avg_log2FC) > 0.5,
                               DEG_CTSB_micro$gene, "")
DEG_CTSB_micro$color <- case_when(
  DEG_CTSB_micro$avg_log2FC > 0.25 & DEG_CTSB_micro$p_val_adj < 0.05 ~ "Up",
  DEG_CTSB_micro$avg_log2FC < -0.25 & DEG_CTSB_micro$p_val_adj < 0.05 ~ "Down",
  TRUE ~ "NS"
)

p_volcano_micro <- ggplot(DEG_CTSB_micro, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = color), size = 1, alpha = 0.6) +
  scale_color_manual(values = c("Up" = "#8B008B", "Down" = "#4393C3", "NS" = "grey80")) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 25) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  labs(title = "CTSB_high vs CTSB_low in AD microglia",
       x = "log2FC", y = "-log10(padj)")
ggsave("results/04_AD_CTSB_high_low_volcano.pdf", p_volcano_micro, width = 9, height = 7)

# GO富集
up_genes_micro <- DEG_CTSB_micro$gene[DEG_CTSB_micro$avg_log2FC > 0.25 & 
                                        DEG_CTSB_micro$p_val_adj < 0.05]

ego_micro <- enrichGO(gene = up_genes_micro,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)

p_go_micro <- dotplot(ego_micro, showCategory = 15,
                      title = "GO-BP: CTSB_high upregulated genes (AD Microglia)")
ggsave("results/04_AD_CTSB_high_GO_BP.pdf", p_go_micro, width = 10, height = 8)

write.csv(ego_micro@result, "results/04_AD_CTSB_high_GO_results.csv", row.names = FALSE)

# KEGG
gene_entrez_micro <- bitr(up_genes_micro, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
options(timeout = 6000000)  # 时间设置长一点
ekegg_micro <- enrichKEGG(gene = gene_entrez_micro$ENTREZID,
                          organism = "hsa", pvalueCutoff = 0.05)

if(nrow(ekegg_micro@result[ekegg_micro@result$p.adjust < 0.05, ]) > 0) {
  p_kegg_micro <- dotplot(ekegg_micro, showCategory = 15,
                          title = "KEGG: CTSB_high upregulated genes (AD micro)")
  ggsave("results/04_AD_CTSB_high_KEGG.pdf", p_kegg_micro, width = 10, height = 8)
}



# ====================================================
# 4.4 ★ 两组CTSB_high上调基因的交集 → 共病通路 ★
# ====================================================

common_CTSB_up <- intersect(up_genes_macro, up_genes_micro)
cat("\n===== CTSB_high共同上调基因（Macro ∩ Micro） =====\n")
cat("巨噬细胞上调：", length(up_genes_macro), "\n")
cat("小胶质上调：", length(up_genes_micro), "\n")
cat("交集：", length(common_CTSB_up), "\n")
print(common_CTSB_up)

# 韦恩图
library(grid)
library(VennDiagram)
#pak::pak("VennDiagram")
pdf("results/04_CTSB_high_common_Venn.pdf", width = 12, height = 6)
grid.newpage()
venn <- venn.diagram(
  x = list("AS Macrophage\nCTSB_high ↑" = up_genes_macro,
           "AD Microglia\nCTSB_high ↑" = up_genes_micro),
  filename = NULL,
  fill = c("#D6604D", "#8B008B"), alpha = 0.5,
  cex = 1.5, cat.cex = 1.1,
  main = "CTSB_high upregulated genes overlap"
)
grid.draw(venn)
dev.off()

# 共同基因的GO富集
if(length(common_CTSB_up) >= 5) {
  ego_common <- enrichGO(gene = common_CTSB_up, OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
  if(nrow(ego_common@result[ego_common@result$p.adjust < 0.05, ]) > 0) {
    p_go_common <- dotplot(ego_common, showCategory = 15,
                           title = "GO: Common CTSB_high upregulated\n(AS Macro ∩ AD Micro)")
    ggsave("results/04_CTSB_high_common_GO.pdf", p_go_common, width = 10, height = 8)
  }
}
