# =============================================================================
# Step 1.6 - KEGG enrichment of the shared targets (clusterProfiler)
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

# Input : interGenes.txt from 04_venn/

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("circlize")
#install.packages("RColorBrewer")
#install.packages("ggpubr")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#BiocManager::install("ComplexHeatmap")


#引用包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)

pvalueFilter=0.05         #p值过滤条件
p.adjustFilter=0.05      #矫正后的p值过滤条件

setwd(file.path(WORK_DIR, "bulk", "05_enrichment"))   # set by 00_config.R
rt=read.table("interGenes.txt", header=F, sep="\t", check.names=F)     #读取输入文件

#提取交集基因的名称,将基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#KEGG富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<p.adjustFilter),]
KEGG=KEGG[KEGG$category != "Human Diseases",]
KEGG=na.omit(KEGG)
#保存显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#设置展示通路的数目
showNum = ifelse(nrow(KEGG)>30, 30, nrow(KEGG))
showKEGG = KEGG[1:showNum,]
showKEGG$pathway = factor(showKEGG$Description, levels = rev(showKEGG$Description))

#绘制柱状图
bar <- ggplot(data = showKEGG, aes(x = Count, y = pathway, fill = pvalue)) +
	geom_bar(stat = "identity", width = 0.75) +
	scale_fill_distiller(palette = "Spectral", direction = 1) +
	labs(x = "Gene count", y = "", title = "Enriched KEGG Pathways") +
	theme_bw()+ theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
		 axis.title = element_text(size = 11), axis.text = element_text(size = 10),
          legend.title = element_text(size = 11), legend.text = element_text(size = 10))
#输出图形
pdf(file="barplot.pdf", width=7, height=5.5)
print(bar)
dev.off()

#绘制气泡图
showKEGG$GeneRatio=sapply(showKEGG$GeneRatio, function(x) eval(parse(text = x)))
showKEGG=showKEGG[order(showKEGG$GeneRatio),]
showKEGG$pathway = factor(showKEGG$pathway, levels = showKEGG$pathway)
bubble <- ggplot(data = showKEGG, aes(x = GeneRatio, y = pathway)) +
	geom_point(aes(size = Count, color = pvalue)) + 
	scale_color_distiller(palette = "Spectral", direction = 1) +
	labs(x = "Gene ratio", y = "", title = "Enriched KEGG Pathways", size = "Count") + 
	theme_bw()+ theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
		axis.title = element_text(size = 11), axis.text = element_text(size = 10),
          legend.title = element_text(size = 11), legend.text = element_text(size = 10))
#输出图形
pdf(file="bubble.pdf", width=7, height=5.5)
print(bubble)
dev.off()
