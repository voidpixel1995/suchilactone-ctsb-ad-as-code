# =============================================================================
# Step 3.1 - Quality control, clustering and SingleR annotation of the scRNA-seq data
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
#BiocManager::install("scrapper")
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

logFCfilter=0.585       #logFC的过滤条件(logFC=0.585,差异倍数1.5倍;logFC=1,差异2倍;logFC=2,差异4倍)
adjPvalFilter=0.05      #矫正后的pvalue的过滤条件
diseaseName="GC"        #疾病组的名称

#设置工作目录
workDir="E:\\1SX\\202604网络药理+心脑血管\\AD\\05.Seurat"
setwd(file.path(WORK_DIR, "singlecell"))   # set by 00_config.R

#读取数据
dirs=list.dirs(workDir)
dirs_sample=dirs[-1]
names(dirs_sample)=gsub(".+\\/(.+)", "\\1", dirs_sample)
counts <- Read10X(data.dir = dirs_sample)
pbmc = CreateSeuratObject(counts, min.cells=3, min.features=100)

#使用PercentageFeatureSet函数计算线粒体基因的百分比
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
#绘制基因特征的小提琴图
pdf(file="01.featureViolin.pdf", width=10, height=6.5)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
dev.off()
pbmc=subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 15)    #对数据进行过滤

#绘制测序深度相关性的图形
pdf(file="01.featureCor.pdf", width=13, height=7)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#对数据进行标准化
pbmc <- NormalizeData(object=pbmc, normalization.method="LogNormalize", scale.factor=10000)
#提取细胞间变异系数较大的基因
pbmc <- FindVariableFeatures(object=pbmc, selection.method="vst", nfeatures=1500)
#输出特征方差图
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="01.featureVar.pdf", width=10, height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()



##PCA分析
pbmc=ScaleData(pbmc)        #PCA降维之前的标准预处理步骤
pbmc=RunPCA(object= pbmc, npcs=20, pc.genes=VariableFeatures(object=pbmc))     #PCA分析
pbmc=RunHarmony(pbmc, "orig.ident")

#绘制每个PCA成分的特征基因
pdf(file="02.pcaGene.pdf", width=10, height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca", nfeatures=20)
dev.off()

#绘制PCA图形
pdf(file="02.PCA.pdf", width=7.5, height=5)
DimPlot(object=pbmc, reduction="pca")
dev.off()

#PCA的热图
pdf(file="02.pcaHeatmap.pdf", width=10, height=8)
DimHeatmap(object=pbmc, dims=1:4, cells=500, balanced=TRUE, nfeatures=30, ncol=2)
dev.off()

#得到每个PC的p值分布
pbmc <- JackStraw(object=pbmc, num.replicate=100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="02.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object=pbmc, dims=1:20)
dev.off()



##聚类分析
pcSelect=20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)     #计算邻接距离

#对细胞分组,对细胞标准模块化
pbmc <- FindClusters(pbmc, resolution=seq(0.5, 1.2, by=0.1))
pbmc <- FindClusters(object = pbmc, resolution=0.6)
#输出聚类图形
pdf(file="03.cluster.pdf", width=7, height=6)
#pbmc <-RunUMAP(object = pbmc, dims = 1:pcSelect)        #UMAP聚类
#DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE)      #UMAP可视化
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)             #TSNE聚类
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)     #TSNE可视化
dev.off()
write.table(pbmc$seurat_clusters,file="03.Cluster.txt",quote=F,sep="\t",col.names=F)

##查找每个聚类的差异基因
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#绘制marker在每个聚类的热图
pdf(file="03.clusterHeatmap.pdf",width=15, height=15)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()



pbmc_for_SingleR <- GetAssayData(pbmc, layer="data")
clusters<-pbmc@meta.data$seurat_clusters
ref1=get(load("ref_Human_all.RData"))
ref2=get(load("ref_Hematopoietic.RData"))
ref3=get(load("DatabaseImmuneCellExpressionData.Rdata"))
ref4=get(load("BlueprintEncode_bpe.se_human.RData"))
ref5=get(load("HumanPrimaryCellAtlas_hpca.se_human.RData"))
ref6=get(load("MonacoImmuneData.Rdata"))
ref7=get(load("NovershternHematopoieticData.Rdata"))
singler=SingleR(test=pbmc_for_SingleR, ref =list(ref1, ref2, ref3, ref4, ref5, ref6, ref7),
              labels=list(ref1$label.main,ref2$label.main,ref3$label.main,ref4$label.main,ref5$label.main,ref6$label.main,ref7$label.main), clusters = clusters)

singler$labels=gsub("_|-", " ", singler$labels)
singler$labels=gsub("T cells, CD4\\+", "CD4+ T cells", singler$labels)
singler$labels=gsub("T cells, CD8\\+", "CD8+ T cells", singler$labels)
singler$labels=gsub("Macrophages", "Macrophage", singler$labels)
clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)

#输出细胞注释的结果
cellAnn <- clusterAnn[match(pbmc$seurat_clusters, clusterAnn[,1]), 2]
cellAnnOut=cbind(names(pbmc$seurat_clusters), cellAnn)
colnames(cellAnnOut)=c("id", "labels")
write.table(cellAnnOut, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

#细胞注释后的可视化
newLabels=singler$labels
names(newLabels)=levels(pbmc)
pbmc=RenameIdents(pbmc, newLabels)
pbmc$cellType=Idents(pbmc)
pdf(file="04.cellAnn.pdf", width=7.5, height=6)
#DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE)     #UMAP可视化
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)           #TSNE可视化
dev.off()
#分组可视化
Type=gsub("(.*?)\\..*", "\\1", colnames(pbmc))
names(Type)=colnames(pbmc)
pbmc=AddMetaData(object=pbmc, metadata=Type, col.name="Type")
pbmc$Sample=ifelse(pbmc$Type=="Control", "Control", diseaseName)
pdf(file="04.group.cellAnn.pdf", width=11, height=6)
#DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE, split.by="Sample")     #UMAP可视化
TSNEPlot(object = pbmc, pt.size = 1, label = TRUE, split.by="Sample")     #TSNE可视化
dev.off()

#细胞类型的差异分析
cell.markers=FindAllMarkers(object = pbmc,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = logFCfilter)
sig.cellMarkers=cell.markers[(abs(as.numeric(as.vector(cell.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(cell.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.cellMarkers,file="04.cellMarkers.txt",sep="\t",row.names=F,quote=F)

#细胞组间差异分析
pbmc.markers=FindMarkers(pbmc, ident.1="Treat", ident.2="Control", group.by='Type', logfc.threshold=0.1)
sig.markersGroup=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
sig.markersGroup=cbind(Gene=row.names(sig.markersGroup), sig.markersGroup)
write.table(sig.markersGroup,file=paste0("05.",diseaseName,"_","all.diffGene.txt"),sep="\t",row.names=F,quote=F)

#细胞组间差异分析(每种细胞单独进行差异分析)
groups=gsub("(.*?)\\..*", "\\1", colnames(pbmc))
groups=paste0(groups, "_", cellAnn)
names(groups)=colnames(pbmc)
pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
allGroupMarkers=data.frame()
for(cellName in unique(cellAnn)){
	conName=paste0("Control_", cellName)
	treatName=paste0("Treat_", cellName)
	if( (length(groups[groups==conName])>10) & (length(groups[groups==treatName])>10) ){
		pbmc.markers=FindMarkers(pbmc, ident.1=treatName, ident.2=conName, group.by='group', logfc.threshold=0.1)
		pbmc.markers$gene=row.names(pbmc.markers)
		pbmc.markers$cluster=cellName
		allGroupMarkers=rbind(allGroupMarkers, pbmc.markers)
		sig.markersGroup=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
		sig.markersGroup=cbind(Gene=row.names(sig.markersGroup), sig.markersGroup)
		write.table(sig.markersGroup,file=paste0("05.",diseaseName,"_",cellName,".diffGene.txt"),sep="\t",row.names=F,quote=F)
	}
}

#保存单细胞数据的对象
save(pbmc, cellAnn, clusterAnn, sig.markers, sig.cellMarkers, allGroupMarkers, file="Seurat.Rdata")
