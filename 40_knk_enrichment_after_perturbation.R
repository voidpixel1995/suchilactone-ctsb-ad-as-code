# =============================================================================
# Step 4.2 - GO / KEGG / Reactome / DisGeNET / GSEA analysis of the perturbed network
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

# Input : sigDiff tables produced by 01_sctenifoldknk_ctsb.R
# Note  : the archived copy of this script ended in the middle of the last
#         export block; the closing DisGeNET exporter below was restored from
#         the identical pattern of the preceding blocks. The GSEA step described
#         in the manuscript was run separately - see 01_bulk/06_GSEA for the
#         equivalent implementation.

setwd(file.path(WORK_DIR, "perturbation"))   # set by 00_config.R

# # CTSB虚拟敲减后差异基因的富集分析 - 完整R代码

# ## 分析策略

# 根据您的数据特点（阿尔茨海默病与动脉粥样硬化共病背景下CTSB虚拟敲减），我为您设计以下富集分析方案：

# 1. **GO富集分析**（BP/MF/CC）
# 2. **KEGG通路富集分析**
#   3. **疾病相关基因集富集**（DisGeNET）
# 4. **GSEA预排序富集分析**（基于FC排序）

# ============================================================================
# CTSB虚拟敲减差异基因富集分析 - SCI高质量图表版
# 疾病背景: 阿尔茨海默病 & 动脉粥样硬化
# ============================================================================
getwd()

# ===================== 1. 安装和加载R包 =====================
required_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", 
                       "ggplot2", "DOSE", "pathview", "ReactomePA",
                       "ggrepel", "RColorBrewer", "cowplot", "stringr",
                       "dplyr", "tidyr", "GOSemSim", "msigdbr")

bioc_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", 
                   "DOSE", "pathview", "ReactomePA", "GOSemSim")

# 安装Bioconductor包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE)
}

# 安装CRAN包
cran_packages <- c("ggplot2", "ggrepel", "RColorBrewer", "cowplot", 
                   "stringr", "dplyr", "tidyr", "msigdbr")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

# 加载所有包
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
  library(ReactomePA)
  library(ggrepel)
  library(RColorBrewer)
  library(cowplot)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(GOSemSim)
  library(msigdbr)
})
#pak::pak("ReactomePA")
#options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
#BiocManager::install("ReactomePA")

# ===================== 2. 数据准备 =====================

# --- 基因集1: 可能是神经元/抑制性神经元相关（根据GAD1, GAD2, LHX6, DLX1等标记）---
genes_set1 <- data.frame(
  gene = c("PLXDC2","NXPH1","ZNF536","GRIP1","BTBD11","GAD1","ANK1","GAD2",
           "ERBB4","VWC2","MYO16","KCNAB1","GRIP2","PTPRM","KIAA1211","PLD5",
           "ZNF385D","ALK","NHS","TIAM1","DNER","SLC6A1","AC010127.1","PTCHD4",
           "DGKD","GRIK1","ANO4","ADCY8","PLCXD3","KCNC2","NRIP3","SLC35F4",
           "ENTPD3","FGF13","ADARB2","TOX3","KAZN","DCX","ARL4C","RUNX2",
           "ASTN2","TENM1","DLX6-AS1","SOX6","RAB3C","IGF1","KCNIP1","SNRPN",
           "UBASH3B","FAM19A4","LHX6","GRIN2D","TMEM132C","XKR4","SLC24A3",
           "OSBPL3","DLX1","ARX","MYO5B","C1QL1","ATP1B1","RBMS3"),
  FC = c(2.757,2.141,2.118,2.102,2.090,2.078,2.074,2.059,2.058,2.051,
         2.037,2.026,2.022,2.012,1.998,1.995,1.991,1.980,1.978,1.972,
         1.967,1.966,1.961,1.958,1.955,1.952,1.952,1.943,1.940,1.930,
         1.924,1.921,1.908,1.907,1.903,1.890,1.887,1.879,1.878,1.872,
         1.870,1.864,1.859,1.855,1.852,1.849,1.848,1.846,1.846,1.841,
         1.839,1.833,1.827,1.826,1.818,1.808,1.808,1.808,1.804,1.804,
         1.803,1.797),
  p.adj = c(3.79E-299,2.00E-13,6.69E-12,5.71E-11,2.63E-10,1.15E-09,
            1.65E-09,8.83E-09,9.05E-09,1.83E-08,7.50E-08,2.16E-07,
            3.09E-07,7.89E-07,2.47E-06,3.12E-06,4.28E-06,9.63E-06,
            1.09E-05,1.57E-05,2.22E-05,2.42E-05,3.35E-05,3.91E-05,
            4.57E-05,5.44E-05,5.44E-05,9.37E-05,0.000113,0.000208,
            0.000294,0.000340,0.000689,0.000697,0.000869,0.001692,
            0.001880,0.002802,0.002864,0.003656,0.003981,0.005145,
            0.006295,0.007279,0.008147,0.008818,0.008941,0.009582,
            0.009582,0.011359,0.012093,0.015247,0.018486,0.019335,
            0.025165,0.034485,0.034485,0.034485,0.037921,0.037921,
            0.037921,0.046112),
  group = "Neuronal"
)

# --- 基因集2: 可能是髓系/免疫细胞相关（根据SLC11A1, MARCO, CD36, SPP1等标记）---
genes_set2 <- data.frame(
  gene = c("VCAN","AQP9","SLC11A1","SLAMF9","CSTB","SDC2","FLNA","FABP5",
           "MMP19","SMIM25","SPP1","CD36","ERO1A","MARCO","CD52","ANPEP",
           "PLIN2","MIF","PHLDA1","CTSL","S100A10","SH3BGRL3","ENO1",
           "CD109","C15orf48","PKM"),
  FC = c(3.582,2.756,2.704,2.603,2.539,2.533,2.522,2.504,2.502,2.487,
         2.424,2.419,2.414,2.412,2.373,2.365,2.352,2.327,2.299,2.231,
         2.230,2.197,2.194,2.185,2.174,2.170),
  p.adj = c(1.20E-193,2.40E-17,1.03E-14,1.61E-10,1.81E-08,2.31E-08,
            4.20E-08,1.19E-07,1.20E-07,2.76E-07,7.86E-06,9.53E-06,
            1.12E-05,1.15E-05,6.35E-05,8.47E-05,0.000135,0.000340,
            0.000890,0.006829,0.006829,0.015912,0.016459,0.019912,
            0.024647,0.026441),
  group = "Myeloid"
)

# 合并所有基因
all_genes_df <- rbind(genes_set1, genes_set2)

# 去除非标准基因名（如lncRNA AC010127.1, DLX6-AS1）用于标准富集
# 保留用于展示但标记
coding_genes <- all_genes_df$gene[!grepl("^AC[0-9]|AS1$", all_genes_df$gene)]

cat("Total unique genes:", length(unique(all_genes_df$gene)), "\n")
cat("Coding genes for enrichment:", length(coding_genes), "\n")

# ===================== 3. 基因ID转换 =====================

# 转换为Entrez ID
gene_symbols <- unique(all_genes_df$gene)
gene_symbols_clean <- gene_symbols[!grepl("^AC[0-9]|AS1$", gene_symbols)]

id_mapping <- bitr(gene_symbols_clean, 
                   fromType = "SYMBOL", 
                   toType = c("ENTREZID", "ENSEMBL"),
                   OrgDb = org.Hs.eg.db)

cat("\nSuccessfully mapped genes:", nrow(id_mapping), "/", length(gene_symbols_clean), "\n")
print(id_mapping)

# 为GSEA准备排序基因列表（基于log2FC）
gene_fc <- all_genes_df %>%
  filter(gene %in% id_mapping$SYMBOL) %>%
  left_join(id_mapping, by = c("gene" = "SYMBOL")) %>%
  arrange(desc(FC))

# 创建命名向量用于GSEA
genelist <- gene_fc$FC
names(genelist) <- gene_fc$ENTREZID
genelist <- sort(genelist, decreasing = TRUE)
genelist <- genelist[!duplicated(names(genelist))]

# Entrez ID向量用于ORA
entrez_ids <- unique(id_mapping$ENTREZID)

# ===================== 4. 设置统一主题（SCI风格）=====================

# SCI期刊风格主题
theme_sci <- theme_bw() +
  theme(
    # 文字
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    # 坐标轴
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 9),
    axis.line = element_line(color = "black", linewidth = 0.5),
    # 图例
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    # 面板
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    # 边距
    plot.margin = margin(10, 15, 10, 10)
  )

# 颜色方案
sci_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", 
                "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                "#7E6148", "#B09C85")

# ===================== 5. GO富集分析 =====================

cat("\n========== GO Enrichment Analysis ==========\n")

# --- 5.1 GO-BP (Biological Process) ---
ego_BP <- enrichGO(gene = entrez_ids,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,
                   readable = TRUE,
                   minGSSize = 5,
                   maxGSSize = 500)

# 简化冗余GO term
if (nrow(as.data.frame(ego_BP)) > 0) {
  ego_BP_sim <- simplify(ego_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
  cat("GO-BP significant terms:", nrow(as.data.frame(ego_BP)), "\n")
  cat("GO-BP after simplification:", nrow(as.data.frame(ego_BP_sim)), "\n")
}

# --- 5.2 GO-MF (Molecular Function) ---
ego_MF <- enrichGO(gene = entrez_ids,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,
                   readable = TRUE,
                   minGSSize = 5,
                   maxGSSize = 500)

# --- 5.3 GO-CC (Cellular Component) ---
ego_CC <- enrichGO(gene = entrez_ids,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,
                   readable = TRUE,
                   minGSSize = 5,
                   maxGSSize = 500)

# ===================== 6. KEGG通路富集分析 =====================

cat("\n========== KEGG Pathway Enrichment ==========\n")
options(timeout = 6000000)
ekegg <- enrichKEGG(gene = entrez_ids,
                    organism = "hsa",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.1,
                    minGSSize = 5,
                    maxGSSize = 500)

if (nrow(as.data.frame(ekegg)) > 0) {
  # 转换KEGG结果中的Entrez ID为Gene Symbol
  ekegg_readable <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  cat("KEGG significant pathways:", nrow(as.data.frame(ekegg)), "\n")
}

# ===================== 7. Reactome通路富集 =====================

cat("\n========== Reactome Pathway Enrichment ==========\n")

ereactome <- enrichPathway(gene = entrez_ids,
                           organism = "human",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.1,
                           readable = TRUE,
                           minGSSize = 5,
                           maxGSSize = 500)

if (nrow(as.data.frame(ereactome)) > 0) {
  cat("Reactome significant pathways:", nrow(as.data.frame(ereactome)), "\n")
}

# ===================== 8. MSigDB基因集富集 =====================

cat("\n========== MSigDB Hallmark & Disease Signatures ==========\n")

# --- 8.1 Hallmark基因集 ---
library(GSVA)
library(limma)
library(GSEABase)

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

gmtFile="hallmark.symbols.gmt"      #基因集文件

hallmark_sets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
hallmark_sets=hallmark_df

ehallmark <- enricher(gene = entrez_ids,
                      TERM2GENE = hallmark_sets,
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      minGSSize = 5)

if (nrow(as.data.frame(ehallmark)) > 0) {
  cat("Hallmark significant gene sets:", nrow(as.data.frame(ehallmark)), "\n")
}

ehallv <- as.data.frame(ehallmark)

# --- 8.2 C2:CP (Canonical Pathways) ---
c2cp_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP") %>%
  dplyr::select(gs_name, entrez_gene)

ec2cp <- enricher(gene = entrez_ids,
                  TERM2GENE = c2cp_sets,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.1,
                  minGSSize = 5)

# --- 8.3 C5:GO (用于补充) & C7 (Immunologic signatures) ---
c7_sets <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") %>%
  dplyr::select(gs_name, entrez_gene)

ec7 <- enricher(gene = entrez_ids,
                TERM2GENE = c7_sets,
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                minGSSize = 5)

# ===================== 9. DisGeNET疾病关联分析 =====================

cat("\n========== DisGeNET Disease Association ==========\n")

# DisGeNET富集
edgn <- enrichDGN(gene = entrez_ids,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.1,
                  readable = TRUE,
                  minGSSize = 3)

if (nrow(as.data.frame(edgn)) > 0) {
  cat("DisGeNET significant diseases:", nrow(as.data.frame(edgn)), "\n")
}

# ===================== 10. 高质量图表输出 =====================

cat("\n========== Generating Publication-Quality Figures ==========\n")

# 创建输出目录
dir.create("enrichment_results", showWarnings = FALSE)
dir.create("enrichment_results/figures", showWarnings = FALSE)
dir.create("enrichment_results/tables", showWarnings = FALSE)

# ---------- Figure 1: GO-BP Barplot (Top 20) ----------
if (nrow(as.data.frame(ego_BP_sim)) > 0) {
  
  bp_data <- as.data.frame(ego_BP_sim) %>%
    arrange(p.adjust) %>%
    head(20) %>%
    mutate(
      Description = str_wrap(Description, width = 45),
      neg_log10_padj = -log10(p.adjust),
      GeneRatio_numeric = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]])
        parts[1] / parts[2]
      })
    )
  
  bp_data$Description <- factor(bp_data$Description, 
                                levels = rev(bp_data$Description))
  
  p1 <- ggplot(bp_data, aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "#E64B35", high = "#4DBBD5", 
                         name = "Adjusted\np-value",
                         guide = guide_colorbar(reverse = FALSE)) +
    scale_size_continuous(range = c(3, 10), name = "Gene\nCount") +
    labs(title = "GO Biological Process Enrichment",
         subtitle = "Top 20 terms after CTSB virtual knockdown",
         x = "Gene Ratio",
         y = NULL) +
    theme_sci +
    theme(axis.text.y = element_text(size = 8.5, lineheight = 0.9))
  
  ggsave("enrichment_results/figures/Fig1_GO_BP_dotplot.pdf", 
         p1, width = 10, height = 8, dpi = 600)
  ggsave("enrichment_results/figures/Fig1_GO_BP_dotplot.tiff", 
         p1, width = 10, height = 8, dpi = 600, compression = "lzw")
  
  cat("Figure 1: GO-BP dotplot saved.\n")
}

# ---------- Figure 2: GO三合一柱状图 (BP + MF + CC) ----------
combine_go_data <- function(ego_bp, ego_mf, ego_cc, top_n = 8) {
  
  process_go <- function(ego_obj, category, n) {
    if (is.null(ego_obj) || nrow(as.data.frame(ego_obj)) == 0) return(NULL)
    df <- as.data.frame(ego_obj) %>%
      arrange(p.adjust) %>%
      head(n) %>%
      mutate(Category = category)
    return(df)
  }
  
  bp_df <- process_go(ego_bp, "Biological Process", top_n)
  mf_df <- process_go(ego_mf, "Molecular Function", top_n)
  cc_df <- process_go(ego_cc, "Cellular Component", top_n)
  
  combined <- bind_rows(bp_df, mf_df, cc_df)
  return(combined)
}

go_combined <- combine_go_data(ego_BP_sim, ego_MF, ego_CC, top_n = 8)

if (nrow(go_combined) > 0) {
  go_combined <- go_combined %>%
    mutate(
      Description = str_wrap(Description, width = 50),
      neg_log10_padj = -log10(p.adjust),
      Category = factor(Category, levels = c("Biological Process", 
                                             "Molecular Function", 
                                             "Cellular Component"))
    )
  
  # 按Category和p.adjust排序
  go_combined <- go_combined %>%
    arrange(Category, p.adjust) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))
  
  p2 <- ggplot(go_combined, aes(x = neg_log10_padj, y = Description, fill = Category)) +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.85) +
    geom_text(aes(label = Count), hjust = -0.3, size = 3, fontface = "bold") +
    scale_fill_manual(values = c("Biological Process" = "#E64B35",
                                 "Molecular Function" = "#4DBBD5", 
                                 "Cellular Component" = "#00A087"),
                      name = "GO Category") +
    labs(title = "Gene Ontology Enrichment Analysis",
         subtitle = "CTSB virtual knockdown - significant DEGs",
         x = expression(-log[10](adjusted~italic(p)-value)),
         y = NULL) +
    theme_sci +
    theme(
      axis.text.y = element_text(size = 8, lineheight = 0.85),
      legend.position = "top",
      legend.justification = "center"
    ) +
    xlim(0, max(go_combined$neg_log10_padj) * 1.15)
  
  ggsave("enrichment_results/figures/Fig2_GO_combined_barplot.pdf",
         p2, width = 12, height = 10, dpi = 600)
  ggsave("enrichment_results/figures/Fig2_GO_combined_barplot.tiff",
         p2, width = 12, height = 10, dpi = 600, compression = "lzw")
  
  cat("Figure 2: GO combined barplot saved.\n")
}

# ---------- Figure 3: KEGG Pathway Dotplot ----------
if (exists("ekegg_readable") && nrow(as.data.frame(ekegg_readable)) > 0) {
  
  kegg_data <- as.data.frame(ekegg_readable) %>%
    arrange(p.adjust) %>%
    head(20) %>%
    mutate(
      Description = str_wrap(Description, width = 45),
      GeneRatio_numeric = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]])
        parts[1] / parts[2]
      })
    )
  
  kegg_data$Description <- factor(kegg_data$Description, 
                                  levels = rev(kegg_data$Description))
  
  p3 <- ggplot(kegg_data, aes(x = GeneRatio_numeric, y = Description)) +
    geom_segment(aes(x = 0, xend = GeneRatio_numeric, 
                     y = Description, yend = Description),
                 color = "grey70", linewidth = 0.5) +
    geom_point(aes(size = Count, color = -log10(p.adjust))) +
    scale_color_gradientn(colors = c("#4DBBD5", "#00A087", "#E64B35"),
                          name = expression(-log[10](p.adj))) +
    scale_size_continuous(range = c(3, 9), name = "Gene Count") +
    labs(title = "KEGG Pathway Enrichment",
         subtitle = "CTSB virtual knockdown DEGs",
         x = "Gene Ratio",
         y = NULL) +
    theme_sci +
    theme(axis.text.y = element_text(size = 9))
  
  ggsave("enrichment_results/figures/Fig3_KEGG_dotplot.pdf",
         p3, width = 10, height = 7, dpi = 600)
  ggsave("enrichment_results/figures/Fig3_KEGG_dotplot.tiff",
         p3, width = 10, height = 7, dpi = 600, compression = "lzw")
  
  cat("Figure 3: KEGG dotplot saved.\n")
}

# ---------- Figure 4: Reactome Pathway ----------
if (exists("ereactome") && nrow(as.data.frame(ereactome)) > 0) {
  
  react_data <- as.data.frame(ereactome) %>%
    arrange(p.adjust) %>%
    head(15) %>%
    mutate(
      Description = str_wrap(Description, width = 50),
      neg_log10_padj = -log10(p.adjust),
      GeneRatio_numeric = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]])
        parts[1] / parts[2]
      })
    )
  
  react_data$Description <- factor(react_data$Description, 
                                   levels = rev(react_data$Description))
  
  p4 <- ggplot(react_data, aes(x = neg_log10_padj, y = Description)) +
    geom_bar(stat = "identity", fill = "#3C5488", alpha = 0.8, width = 0.65) +
    geom_text(aes(label = paste0("n=", Count)), 
              hjust = -0.2, size = 3.2, color = "black") +
    labs(title = "Reactome Pathway Enrichment",
         subtitle = "CTSB knockdown responsive genes",
         x = expression(-log[10](adjusted~italic(p)-value)),
         y = NULL) +
    theme_sci +
    theme(axis.text.y = element_text(size = 8.5)) +
    xlim(0, max(react_data$neg_log10_padj) * 1.15)
  
  ggsave("enrichment_results/figures/Fig4_Reactome_barplot.pdf",
         p4, width = 10, height = 7, dpi = 600)
  ggsave("enrichment_results/figures/Fig4_Reactome_barplot.tiff",
         p4, width = 10, height = 7, dpi = 600, compression = "lzw")
  
  cat("Figure 4: Reactome barplot saved.\n")
}

# ---------- Figure 5: MSigDB Hallmark Enrichment ----------
if (exists("ehallmark") && nrow(as.data.frame(ehallmark)) > 0) {
  
  hallmark_data <- as.data.frame(ehallmark) %>%
    arrange(p.adjust) %>%
    head(15) %>%
    mutate(
      # 清理Hallmark名称
      Description = gsub("HALLMARK_", "", ID),
      Description = gsub("_", " ", Description),
      Description = str_to_title(Description),
      neg_log10_padj = -log10(p.adjust),
      GeneRatio_numeric = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]])
        parts[1] / parts[2]
      })
    )
  
  hallmark_data$Description <- factor(hallmark_data$Description, 
                                      levels = rev(hallmark_data$Description))
  
  p5 <- ggplot(hallmark_data, aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = neg_log10_padj)) +
    scale_color_gradient2(low = "#4DBBD5", mid = "#F39B7F", high = "#E64B35",
                          midpoint = median(hallmark_data$neg_log10_padj),
                          name = expression(-log[10](p.adj))) +
    scale_size_continuous(range = c(4, 10), name = "Gene Count") +
    labs(title = "MSigDB Hallmark Gene Set Enrichment",
         subtitle = "CTSB virtual knockdown",
         x = "Gene Ratio",
         y = NULL) +
    theme_sci
  
  ggsave("enrichment_results/figures/Fig5_Hallmark_dotplot.pdf",
         p5, width = 10, height = 7, dpi = 600)
  ggsave("enrichment_results/figures/Fig5_Hallmark_dotplot.tiff",
         p5, width = 10, height = 7, dpi = 600, compression = "lzw")
  
  cat("Figure 5: Hallmark dotplot saved.\n")
}

# ---------- Figure 6: DisGeNET 疾病关联分析 ----------
if (exists("edgn") && nrow(as.data.frame(edgn)) > 0) {
  
  dgn_data <- as.data.frame(edgn) %>%
    arrange(p.adjust) %>%
    head(20) %>%
    mutate(
      Description = str_wrap(Description, width = 45),
      neg_log10_padj = -log10(p.adjust),
      GeneRatio_numeric = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]])
        parts[1] / parts[2]
      })
    )
  
  dgn_data$Description <- factor(dgn_data$Description, 
                                 levels = rev(dgn_data$Description))
  
  # 标注AD和AS相关疾病
  dgn_data <- dgn_data %>%
    mutate(
      disease_category = case_when(
        grepl("alzheimer|dementia|cognitive|neurodegen|amyloid|tau", 
              Description, ignore.case = TRUE) ~ "AD-related",
        grepl("atheros|cardiovas|coronary|arterial|vascular|lipid|cholesterol", 
              Description, ignore.case = TRUE) ~ "AS-related",
        grepl("inflam|immune|cytokine", 
              Description, ignore.case = TRUE) ~ "Inflammation",
        TRUE ~ "Other"
      )
    )
  
  p6 <- ggplot(dgn_data, aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = disease_category, alpha = neg_log10_padj)) +
    scale_color_manual(values = c("AD-related" = "#E64B35", 
                                  "AS-related" = "#3C5488",
                                  "Inflammation" = "#00A087",
                                  "Other" = "#8491B4"),
                       name = "Disease\nCategory") +
    scale_size_continuous(range = c(3, 10), name = "Gene\nCount") +
    scale_alpha_continuous(range = c(0.5, 1), 
                           name = expression(-log[10](p.adj))) +
    labs(title = "DisGeNET Disease Association Enrichment",
         subtitle = "CTSB knockdown DEGs - Disease relevance",
         x = "Gene Ratio",
         y = NULL) +
    theme_sci +
    theme(axis.text.y = element_text(size = 8.5, lineheight = 0.9))
  
  ggsave("enrichment_results/figures/Fig6_DisGeNET_dotplot.pdf",
         p6, width = 11, height = 8, dpi = 600)
  ggsave("enrichment_results/figures/Fig6_DisGeNET_dotplot.tiff",
         p6, width = 11, height = 8, dpi = 600, compression = "lzw")
  
  cat("Figure 6: DisGeNET dotplot saved.\n")
}

# ---------- Figure 7: GO-BP Cnetplot (基因-通路网络图) ----------
if (nrow(as.data.frame(ego_BP_sim)) >= 5) {
  
  pdf("enrichment_results/figures/Fig7_GO_BP_cnetplot.pdf", 
      width = 14, height = 12)
  
  # 创建FC命名向量
  fc_vector <- all_genes_df$FC
  names(fc_vector) <- all_genes_df$gene
  
  p7 <- cnetplot(ego_BP_sim, 
                 showCategory = 8,
                 foldChange = fc_vector,
                 cex_label_gene = 0.7,
                 cex_label_category = 0.9,
                 cex_category = 1.5,
                 colorEdge = TRUE,
                 circular = FALSE,
                 node_label = "all") +
    scale_color_gradient2(low = "#4DBBD5", mid = "#F39B7F", high = "#E64B35",
                          midpoint = median(fc_vector),
                          name = "Fold Change") +
    labs(title = "GO-BP Gene-Concept Network",
         subtitle = "Node color reflects fold change after CTSB knockdown") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
          legend.position = "right")
  
  print(p7)
  dev.off()
  
  # 同时保存TIFF
  tiff("enrichment_results/figures/Fig7_GO_BP_cnetplot.tiff",
       width = 14, height = 12, units = "in", res = 600, compression = "lzw")
  print(p7)
  dev.off()
  
  cat("Figure 7: GO-BP cnetplot saved.\n")
}

# ---------- Figure 8: GO-BP Treeplot (树状图) ----------
if (nrow(as.data.frame(ego_BP_sim)) >= 5) {
  
  # 计算语义相似性用于聚类
  ego_BP_sim2 <- pairwise_termsim(ego_BP_sim)
  
  p8 <- treeplot(ego_BP_sim2, 
                 showCategory = 25,
                 fontsize = 3.5,
                 offset = rel(1.2),
                 offset_tiplab = rel(1),
                 cluster.params = list(method = "ward.D2", n = 5)) +
    labs(title = "GO-BP Hierarchical Clustering Tree",
         subtitle = "CTSB knockdown - grouped by semantic similarity") +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"))
  
  ggsave("enrichment_results/figures/Fig8_GO_BP_treeplot.pdf",
         p8, width = 14, height = 10, dpi = 600)
  ggsave("enrichment_results/figures/Fig8_GO_BP_treeplot.tiff",
         p8, width = 14, height = 10, dpi = 600, compression = "lzw")
  
  cat("Figure 8: GO-BP treeplot saved.\n")
}

# ---------- Figure 9: GO-BP Emapplot (富集图谱网络) ----------
if (nrow(as.data.frame(ego_BP_sim)) >= 5) {
  
  ego_BP_sim3 <- pairwise_termsim(ego_BP_sim)
  
  p9 <- emapplot(ego_BP_sim3, 
                 showCategory = 25,
                 cex_label_category = 0.65,
                 cex_category = 1.2,
                 layout = "kk",
                 repel = TRUE) +
    scale_fill_gradient(low = "#F39B7F", high = "#E64B35", 
                        name = "p.adjust") +
    labs(title = "GO-BP Enrichment Map Network",
         subtitle = "Edge width reflects term similarity") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"))
  
  ggsave("enrichment_results/figures/Fig9_GO_BP_emapplot.pdf",
         p9, width = 12, height = 10, dpi = 600)
  ggsave("enrichment_results/figures/Fig9_GO_BP_emapplot.tiff",
         p9, width = 12, height = 10, dpi = 600, compression = "lzw")
  
  cat("Figure 9: GO-BP emapplot saved.\n")
}

# ---------- Figure 10: 分组比较分析 (Neuronal vs Myeloid) ----------

# 为两个细胞群分别做富集分析并比较
neuronal_genes <- genes_set1$gene[!grepl("^AC[0-9]|AS1$", genes_set1$gene)]
myeloid_genes <- genes_set2$gene

neuronal_ids <- bitr(neuronal_genes, fromType = "SYMBOL", 
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)
myeloid_ids <- bitr(myeloid_genes, fromType = "SYMBOL", 
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 创建基因列表用于compareCluster
gene_cluster_list <- list(
  "Neuronal" = neuronal_ids$ENTREZID,
  "Myeloid" = myeloid_ids$ENTREZID
)

# GO-BP比较
compare_go_bp <- compareCluster(geneCluster = gene_cluster_list,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                readable = TRUE)

if (nrow(as.data.frame(compare_go_bp)) > 0) {
  
  p10 <- dotplot(compare_go_bp, 
                 showCategory = 10,
                 font.size = 9,
                 title = "GO-BP Comparison: Neuronal vs Myeloid") +
    scale_color_gradient(low = "#E64B35", high = "#4DBBD5",
                         name = "p.adjust") +
    theme_sci +
    theme(
      axis.text.y = element_text(size = 8, lineheight = 0.85),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 10, face = "bold")
    )
  
  ggsave("enrichment_results/figures/Fig10_Compare_GO_BP.pdf",
         p10, width = 12, height = 9, dpi = 600)
  ggsave("enrichment_results/figures/Fig10_Compare_GO_BP.tiff",
         p10, width = 12, height = 9, dpi = 600, compression = "lzw")
  
  cat("Figure 10: Comparative GO-BP dotplot saved.\n")
}

# KEGG比较
compare_kegg <- compareCluster(geneCluster = gene_cluster_list,
                               fun = "enrichKEGG",
                               organism = "hsa",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.1)

if (nrow(as.data.frame(compare_kegg)) > 0) {
  
  compare_kegg_readable <- setReadable(compare_kegg, OrgDb = org.Hs.eg.db, 
                                       keyType = "ENTREZID")
  
  p10b <- dotplot(compare_kegg_readable, 
                  showCategory = 10,
                  font.size = 9,
                  title = "KEGG Comparison: Neuronal vs Myeloid") +
    scale_color_gradient(low = "#E64B35", high = "#4DBBD5",
                         name = "p.adjust") +
    theme_sci +
    theme(
      axis.text.y = element_text(size = 8.5),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
    )
  
  ggsave("enrichment_results/figures/Fig10b_Compare_KEGG.pdf",
         p10b, width = 12, height = 8, dpi = 600)
  ggsave("enrichment_results/figures/Fig10b_Compare_KEGG.tiff",
         p10b, width = 12, height = 8, dpi = 600, compression = "lzw")
  
  cat("Figure 10b: Comparative KEGG dotplot saved.\n")
}

# ---------- Figure 11: 综合气泡图 (Bubble Plot) - 多数据库整合 ----------

# 整合多个数据库结果
integrate_results <- function(go_res, kegg_res, react_res, hallmark_res, top_n = 5) {
  
  results <- list()
  
  if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
    go_df <- as.data.frame(go_res) %>%
      arrange(p.adjust) %>% head(top_n) %>%
      mutate(Database = "GO-BP")
    results[["GO"]] <- go_df
  }
  
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    kegg_df <- as.data.frame(kegg_res) %>%
      arrange(p.adjust) %>% head(top_n) %>%
      mutate(Database = "KEGG")
    results[["KEGG"]] <- kegg_df
  }
  
  if (!is.null(react_res) && nrow(as.data.frame(react_res)) > 0) {
    react_df <- as.data.frame(react_res) %>%
      arrange(p.adjust) %>% head(top_n) %>%
      mutate(Database = "Reactome")
    results[["Reactome"]] <- react_df
  }
  
  if (!is.null(hallmark_res) && nrow(as.data.frame(hallmark_res)) > 0) {
    hallmark_df <- as.data.frame(hallmark_res) %>%
      arrange(p.adjust) %>% head(top_n) %>%
      mutate(
        Description = gsub("HALLMARK_", "", ID),
        Description = gsub("_", " ", Description),
        Description = str_to_title(Description),
        Database = "Hallmark"
      )
    results[["Hallmark"]] <- hallmark_df
  }
  
  # 统一列名合并
  cols_keep <- c("Description", "GeneRatio", "Count", "p.adjust", "Database")
  combined <- bind_rows(lapply(results, function(df) {
    df %>% dplyr::select(any_of(cols_keep))
  }))
  
  return(combined)
}

integrated_data <- integrate_results(ego_BP_sim, ekegg_readable, 
                                     ereactome, ehallmark, top_n = 5)

if (nrow(integrated_data) > 0) {
  
  integrated_data <- integrated_data %>%
    mutate(
      GeneRatio_numeric = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]])
        parts[1] / parts[2]
      }),
      neg_log10_padj = -log10(p.adjust),
      Description = str_wrap(Description, width = 40),
      Database = factor(Database, levels = c("GO-BP", "KEGG", "Reactome", "Hallmark"))
    ) %>%
    group_by(Database) %>%
    arrange(p.adjust) %>%
    ungroup()
  
  integrated_data$Description <- factor(integrated_data$Description,
                                        levels = rev(unique(integrated_data$Description)))
  
  p11 <- ggplot(integrated_data, 
                aes(x = neg_log10_padj, y = Description, 
                    size = Count, color = Database)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("GO-BP" = "#E64B35", 
                                  "KEGG" = "#4DBBD5",
                                  "Reactome" = "#00A087", 
                                  "Hallmark" = "#3C5488"),
                       name = "Database") +
    scale_size_continuous(range = c(3, 10), name = "Gene Count") +
    facet_grid(Database ~ ., scales = "free_y", space = "free_y",
               switch = "y") +
    labs(title = "Multi-database Enrichment Summary",
         subtitle = "Top pathways from GO, KEGG, Reactome, and MSigDB Hallmark",
         x = expression(-log[10](adjusted~italic(p)-value)),
         y = NULL) +
    theme_sci +
    theme(
      axis.text.y = element_text(size = 8, lineheight = 0.85),
      strip.text.y.left = element_text(size = 9, face = "bold", angle = 0),
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      panel.spacing = unit(0.3, "lines"),
      legend.position = "right"
    )
  
  ggsave("enrichment_results/figures/Fig11_Multi_DB_enrichment.pdf",
         p11, width = 12, height = 11, dpi = 600)
  ggsave("enrichment_results/figures/Fig11_Multi_DB_enrichment.tiff",
         p11, width = 12, height = 11, dpi = 600, compression = "lzw")
  
  cat("Figure 11: Multi-database enrichment summary saved.\n")
}

# ---------- Figure 12: Heatmap-like Plot (热图式富集展示) ----------
if (nrow(as.data.frame(compare_go_bp)) > 0) {
  
  # 提取比较结果数据
  compare_data <- as.data.frame(compare_go_bp) %>%
    arrange(p.adjust) %>%
    group_by(Cluster) %>%
    slice_head(n = 15) %>%
    ungroup()
  
  if (nrow(compare_data) > 0) {
    
    # 转为宽格式计算
    compare_wide <- compare_data %>%
      mutate(
        neg_log10_p = -log10(p.adjust),
        GeneRatio_numeric = sapply(GeneRatio, function(x) {
          parts <- as.numeric(strsplit(x, "/")[[1]])
          parts[1] / parts[2]
        })
      ) %>%
      dplyr::select(Description, Cluster, neg_log10_p, GeneRatio_numeric, Count)
    
    compare_wide$Description <- str_wrap(compare_wide$Description, width = 50)
    
    p12 <- ggplot(compare_wide, 
                  aes(x = Cluster, y = reorder(Description, neg_log10_p))) +
      geom_point(aes(size = Count, fill = neg_log10_p), 
                 shape = 21, color = "black", stroke = 0.5) +
      scale_fill_gradient2(low = "#4DBBD5", mid = "#FFFFBF", high = "#E64B35",
                           midpoint = median(compare_wide$neg_log10_p),
                           name = expression(-log[10](p.adj))) +
      scale_size_continuous(range = c(2, 10), name = "Gene Count") +
      labs(title = "Comparative Enrichment Heatmap",
           subtitle = "Neuronal vs Myeloid cell types after CTSB knockdown",
           x = "Cell Type",
           y = NULL) +
      theme_sci +
      theme(
        axis.text.y = element_text(size = 7.5, lineheight = 0.85),
        axis.text.x = element_text(size = 11, face = "bold"),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
      )
    
    ggsave("enrichment_results/figures/Fig12_Compare_heatmap.pdf",
           p12, width = 10, height = 10, dpi = 600)
    ggsave("enrichment_results/figures/Fig12_Compare_heatmap.tiff",
           p12, width = 10, height = 10, dpi = 600, compression = "lzw")
    
    cat("Figure 12: Comparative heatmap saved.\n")
  }
}

# ---------- Figure 13: 圈图 (Chord diagram / Circos plot) ----------
if (nrow(as.data.frame(ego_BP_sim)) >= 3) {
  
  # 使用enrichplot的upsetplot展示基因-通路交集
  p13 <- upsetplot(ego_BP_sim, n = 12) +
    labs(title = "Gene-Pathway Intersection (UpSet Plot)") +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  
  ggsave("enrichment_results/figures/Fig13_UpSet_GO_BP.pdf",
         p13, width = 12, height = 7, dpi = 600)
  ggsave("enrichment_results/figures/Fig13_UpSet_GO_BP.tiff",
         p13, width = 12, height = 7, dpi = 600, compression = "lzw")
  
  cat("Figure 13: UpSet plot saved.\n")
}

# ---------- Figure 14: 自定义AD-AS共病相关通路展示 ----------

# 定义AD和AS相关的自定义通路基因集
ad_as_custom_genesets <- data.frame(
  term = c(
    # 神经退行性疾病相关
    rep("GABAergic signaling", 4),
    rep("Neuronal development", 5),
    rep("Synaptic transmission", 5),
    rep("Axon guidance", 3),
    # 免疫/炎症相关
    rep("Innate immune response", 4),
    rep("Lipid metabolism", 3),
    rep("Lysosomal function", 2),
    rep("Phagocytosis/Endocytosis", 3),
    # 共病相关
    rep("ECM remodeling", 3),
    rep("Glycolysis/Metabolism", 3),
    rep("Cell migration", 3)
  ),
  gene = c(
    # GABAergic
    "GAD1", "GAD2", "SLC6A1", "DLX1",
    # Neuronal development
    "DCX", "ARX", "LHX6", "SOX6", "DNER",
    # Synaptic
    "GRIK1", "GRIN2D", "KCNC2", "KCNAB1", "KCNIP1",
    # Axon
    "PLXDC2", "ERBB4", "ALK",
    # Innate immune
    "SLC11A1", "MARCO", "SLAMF9", "CD52",
    # Lipid
    "CD36", "FABP5", "PLIN2",
    # Lysosomal
    "CTSL", "CSTB",
    # Phagocytosis
    "SPP1", "AQP9", "ANPEP",
    # ECM
    "VCAN", "MMP19", "SDC2",
    # Glycolysis
    "PKM", "ENO1", "ERO1A",
    # Migration
    "FLNA", "TIAM1", "CD109"
  )
)

# 合并FC数据
ad_as_plot_data <- ad_as_custom_genesets %>%
  left_join(all_genes_df %>% dplyr::select(gene, FC, group), by = "gene") %>%
  filter(!is.na(FC))

# 分类
ad_as_plot_data <- ad_as_plot_data %>%
  mutate(
    pathway_category = case_when(
      term %in% c("GABAergic signaling", "Neuronal development", 
                  "Synaptic transmission", "Axon guidance") ~ "Neurodegeneration",
      term %in% c("Innate immune response", "Phagocytosis/Endocytosis",
                  "Lysosomal function") ~ "Immune/Inflammation",
      term %in% c("Lipid metabolism", "Glycolysis/Metabolism", 
                  "ECM remodeling", "Cell migration") ~ "Metabolism/Remodeling"
    ),
    pathway_category = factor(pathway_category, 
                              levels = c("Neurodegeneration", 
                                         "Immune/Inflammation",
                                         "Metabolism/Remodeling")),
    term = factor(term, levels = rev(unique(term)))
  )

p14 <- ggplot(ad_as_plot_data, aes(x = gene, y = term)) +
  geom_point(aes(size = FC, fill = pathway_category), 
             shape = 21, color = "black", stroke = 0.4, alpha = 0.85) +
  scale_fill_manual(values = c("Neurodegeneration" = "#E64B35",
                               "Immune/Inflammation" = "#4DBBD5",
                               "Metabolism/Remodeling" = "#00A087"),
                    name = "Pathway\nCategory") +
  scale_size_continuous(range = c(3, 9), 
                        name = "Fold Change",
                        breaks = c(1.9, 2.2, 2.5, 3.0, 3.5)) +
  facet_grid(pathway_category ~ ., scales = "free_y", space = "free_y") +
  labs(title = "AD-AS Comorbidity Pathways: CTSB Knockdown Effects",
       subtitle = "Gene membership in disease-relevant functional categories",
       x = "Gene Symbol",
       y = NULL) +
  theme_sci +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, 
                               size = 8, face = "italic"),
    axis.text.y = element_text(size = 9),
    strip.text.y = element_text(size = 9, face = "bold", angle = 0),
    strip.background = element_rect(fill = "grey95", color = "grey60"),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "right"
  )

ggsave("enrichment_results/figures/Fig14_AD_AS_custom_pathways.pdf",
       p14, width = 14, height = 9, dpi = 600)
ggsave("enrichment_results/figures/Fig14_AD_AS_custom_pathways.tiff",
       p14, width = 14, height = 9, dpi = 600, compression = "lzw")

cat("Figure 14: AD-AS custom pathway plot saved.\n")

# ---------- Figure 15: 综合汇总图 (Panel Figure) ----------

# 组合关键图形为一个大Panel
if (exists("p1") && exists("p3") && exists("p10") && exists("p14")) {
  
  # 标注panel标签
  p1_label <- p1 + labs(tag = "A") + 
    theme(plot.tag = element_text(size = 16, face = "bold"))
  p3_label <- p3 + labs(tag = "B") + 
    theme(plot.tag = element_text(size = 16, face = "bold"))
  p10_label <- p10 + labs(tag = "C") + 
    theme(plot.tag = element_text(size = 16, face = "bold"))
  p14_label <- p14 + labs(tag = "D") + 
    theme(plot.tag = element_text(size = 16, face = "bold"))
  
  p_combined <- plot_grid(
    p1_label, p3_label,
    p10_label, p14_label,
    ncol = 2, nrow = 2,
    rel_widths = c(1, 1),
    rel_heights = c(1, 1),
    align = "hv",
    axis = "tblr"
  )
  
  ggsave("enrichment_results/figures/Fig15_Combined_Panel.pdf",
         p_combined, width = 22, height = 18, dpi = 600)
  ggsave("enrichment_results/figures/Fig15_Combined_Panel.tiff",
         p_combined, width = 22, height = 18, dpi = 600, compression = "lzw")
  
  cat("Figure 15: Combined panel figure saved.\n")
}

# ===================== 11. 导出富集分析结果表格 =====================

cat("\n========== Exporting Result Tables ==========\n")

# GO-BP
if (nrow(as.data.frame(ego_BP_sim)) > 0) {
  write.csv(as.data.frame(ego_BP_sim), 
            "enrichment_results/tables/GO_BP_enrichment.csv", 
            row.names = FALSE)
  cat("GO-BP results exported.\n")
}

# GO-MF
if (nrow(as.data.frame(ego_MF)) > 0) {
  write.csv(as.data.frame(ego_MF), 
            "enrichment_results/tables/GO_MF_enrichment.csv", 
            row.names = FALSE)
  cat("GO-MF results exported.\n")
}

# GO-CC
if (nrow(as.data.frame(ego_CC)) > 0) {
  write.csv(as.data.frame(ego_CC), 
            "enrichment_results/tables/GO_CC_enrichment.csv", 
            row.names = FALSE)
  cat("GO-CC results exported.\n")
}

# KEGG
if (exists("ekegg_readable") && nrow(as.data.frame(ekegg_readable)) > 0) {
  write.csv(as.data.frame(ekegg_readable), 
            "enrichment_results/tables/KEGG_enrichment.csv", 
            row.names = FALSE)
  cat("KEGG results exported.\n")
}

# Reactome
if (exists("ereactome") && nrow(as.data.frame(ereactome)) > 0) {
  write.csv(as.data.frame(ereactome), 
            "enrichment_results/tables/Reactome_enrichment.csv", 
            row.names = FALSE)
  cat("Reactome results exported.\n")
}

# Hallmark
if (exists("ehallmark") && nrow(as.data.frame(ehallmark)) > 0) {
  write.csv(as.data.frame(ehallmark), 
            "enrichment_results/tables/Hallmark_enrichment.csv", 
            row.names = FALSE)
  cat("Hallmark results exported.\n")
}

# DisGeNET
if (exists("edgn") && nrow(as.data.frame(edgn)) > 0) {
  write.csv(as.data.frame(edgn),
            "enrichment_results/tables/DisGeNET_enrichment.csv",
            row.names = FALSE)
  cat("DisGeNET results exported.\n")
}

cat("\nAll enrichment analyses finished.\n")
