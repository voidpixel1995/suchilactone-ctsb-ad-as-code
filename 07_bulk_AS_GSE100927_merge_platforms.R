# =============================================================================
# Step 1.1b - Merge GSE100927 sub-series and map probes to gene symbols
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

####下载geo的数据####

#先解压GSE42872_series_matrix.txt.gz，注意解压到当前目录，再读入
#comment.char="!" 意思是！后面的内容不要读取，可以打开文件看一下

setwd(file.path(WORK_DIR, "bulk", "AS", "00_normalize", "GSE100927"))   # set by 00_config.R

exprSet <- read.table("GSE100927_series_matrix.txt",
                      comment.char="!",
                      stringsAsFactors=F,
                      header=T)

#高频操作来了：第一列变成行名
#分两步操作
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]

####自动log化与校正####
#注意此时数据的格式是矩阵
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

## 开始判断
if (LogC) { 
  # 将小于或等于0的值替换为NaN
  ex[ex <= 0] <- NA  # 使用NA而不是NaN
  # 取log2
  exprSet <- log2(ex)
  print("log2 transform finished")
} else {
  print("log2 transform not needed")
}

library(limma) 
boxplot(exprSet,outline=FALSE, notch=T, las=2)
### 该函数默认使用quntile 矫正差异 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)
## 这步把矩阵转换为数据框很重要
class(exprSet)
exprSet <- as.data.frame(exprSet)


####探针基因名转换####
#platformMap 中有常见的平台个R注释包的对应关系
platformMap <- data.table::fread("platformMap.txt",data.table = F)

#平台的名称如何知道?
index <- "GPL6244"
#数据储存在bioc_package这一列中
paste0(platformMap$bioc_package[grep(index,platformMap$gpl)],".db")

#安装R包,可以直接安装，这里用了判断
if(!requireNamespace("hugene10sttranscriptcluster.db")){
  options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
  BiocManager::install("hugene10sttranscriptcluster.db",update = F,ask = F)
  } 

#加载R包
library(hugene10sttranscriptcluster.db)
#获取探针和基因的对应关系：这是探针注释的关键步骤
probe2symbol_df <- toTable(get("hugene10sttranscriptclusterSYMBOL"))
#探针有多少个？
length(unique(probe2symbol_df$probe_id))
#这么多行中，基因名称有重复的么？
length(unique(probe2symbol_df$symbol))

### 探针转换以及去重，获得最终的表达矩阵
### 拆分体会,要学会分步探索排查，#号的使用
library(dplyr)
library(tibble)
exprSet <- exprSet %>% 
  ## 行名转列名,因为只有变成数据框的列,才可以用inner_join
  rownames_to_column("probe_id") %>% 
  ## 合并探针的信息
  inner_join(probe2symbol_df,by="probe_id") %>% 
  ## 去掉多余信息
  select(-probe_id) %>%  
  ## 重新排列
  select(symbol,everything()) %>%  
  ## rowMeans求出行的平均数(这边的.代表上面传入的数据)
  ## .[,-1]表示去掉出入数据的第一列，然后求行的平均值
  mutate(rowMean =rowMeans(.[,-1])) %>% 
  ## 把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  ## 去重，symbol留下第一个
  distinct(symbol,.keep_all = T) %>% 
  ## 反向选择去除rowMean这一列
  select(-rowMean) %>% 
  ## 列名转行名
  column_to_rownames("symbol")

exprSet1 <- select(exprSet,-c(1:46))
exprSet <- exprSet1
### 保存数据Rdata格式
save(exprSet,file = "exprSet_rmdup.Rdata")

### 保存数据csv格式
write.csv(exprSet1,file = "GSE25097.csv")

### 保存为txt格式
# 将行名与数据框结合，行名成为第一列
exprSet <- cbind(Rownames = rownames(exprSet), exprSet)
#修改第一列名称
colnames(exprSet)[1] <- "id"
#输出txt格式
write.table(exprSet, file="GSE25097geneMatrix.txt", 
            sep="\t", quote=F, col.names=T, row.names=F)















### 使用limma来做芯片的差异分析
#加载limma包，用于校正和比较差异
rm(list = ls())
library(limma)
### 加载数据，注意解决报错
load(file = "output/exprSet_rmdup.Rdata")

### 1.创建分组
### 这一步根据样本来就行，原则就是: 跟样本匹配，取决于样本的排序
group <- c(rep("treat",8),rep("con",20)) 
##group <- c("treat","con","treat","con","treat","treat") 
### 分组变成向量，并且限定leves的顺序
### levels里面，把对照组放在前面
group <- factor(group,levels = c("con","treat"))

### 主成分分析PCA：提前预测结果
### 行是样本列是基因
res.pca <- prcomp(t(exprSet), scale = TRUE)

###如果报错用下面这段代码
if(T){
  # 1.转置
  data <- t(exprSet)
  ## 2. NAN值变成0
  data[is.nan(data)] <- 0
  ##3.删除没有变异的基因，使用函数 var
  data <- data[,which(apply(data,2,var) != 0)]
  res.pca <- prcomp(data,scale = TRUE)
}

library(factoextra)
fviz_pca_ind(res.pca,col.ind = group)

### 构建比较矩阵
design <- model.matrix(~group)
### 比较矩阵命名
colnames(design) <- levels(group)
design

### 2.线性模型拟合
fit <- lmFit(exprSet,design)
### 3.贝叶斯检验
fit2 <- eBayes(fit)
### 4.输出差异分析结果,其中coef的数目不能操过design的列数
### 此处的2代表的是design中第二列和第一列的比较
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
### 这个数据很重要需要保存一下
save(allDiff,file = "output/allDiff.Rdata")
### 定义差异基因：差异倍数2倍，矫正后的p值小于0.05

library(dplyr)
diffgene <- allDiff %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >1)

### 如果出现行名丢失的情况，需要先把行名变成列，处理完毕后再把列变成行名
### 这个工作是由tibble这个包里面的rownames_to_column()和column_to_rownames()完成的
library(tibble)
diffgene <- allDiff %>% 
  rownames_to_column() %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >1) %>% 
  column_to_rownames()

### 可选方案:使用subset直接获取,&是and的意思
diffgene <- subset(allDiff,abs(logFC) >1 & adj.P.Val < 0.05)
test <- allDiff[allDiff$adj.P.Val < 0.05 & abs(allDiff$logFC)>1,]
### 该数据也需要保存，此处一次性保存两个数据，如果是多个，一次写入变量名称即可。
save(diffgene,group,file = "output/diffgene.Rdata")
### 到此差异基因的分析就结束了
## 作图环节
## 1.把现在数据调整成可以作图的格式
### 这个技能是data wrangling部分重点掌握的技能
### 复习一下流程：输入数据是表达量，经过三步
### 1.探针ID转换，2.行列转置，3，添加分组信息。最终获得的是数据框

### 行列转置
exprSet <- as.data.frame(t(exprSet))
### 添加分组信息
dd <- cbind(group=group,exprSet)
### 截取部分展示,这就是清洁数据
test = dd[,1:10]

## 2.作图展示
library(ggplot2)
ggplot(data = dd,aes(x=group,y=CD36,fill=group))+
  geom_boxplot()+
  geom_point()+
  theme_bw()

## 3.steal plot
my_comparisons <- list(
  c("con", "treat")
)
library(ggpubr)
ggboxplot(
  dd, x = "group", y = "CD36",
  color = "group", palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

## 改写成函数
diffplot <- function(gene){
  my_comparisons <- list(
    c("treat", "con")
  )
  library(ggpubr)
  ggboxplot(
    dd, x = "group", y = gene,
    color = "group", palette = c("#00AFBB", "#E7B800"),
    add = "jitter"
  )+
    stat_compare_means(comparisons = my_comparisons, method = "t.test")
}

diffplot("CD36")
diffplot("MOXD1")
diffplot("NOSIP")
## 4.多个基因作图查看
## 先把基因提取出来
genelist <- rownames(diffgene)[1:6]
## 再提取表达量，使用名称选取行
data <- dd[,c("group",genelist)]
## 用pivot_longer调整数据，数据变长，增加的是行
library(tidyr)
data <- data %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "expression")
## 多基因作图
## 作图
ggplot(data = data,aes(x=gene,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  stat_compare_means(aes(group=group), label = "p.signif", method = "t.test")

## 尝试更清晰的展示
pdf(file = "output/diffgenboxplot.pdf", width = 8, height = 6)
ggplot(data = data,aes(x=group,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter(size =0.5)+
  theme_bw()+
  facet_grid(.~gene)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")
dev.off()

## 图片导出
library(export)
## 导成PPT可编辑的格式
graph2ppt(file="output/diffgenboxplot.pptx")
## 其他自己想要的格式
graph2pdf(file="output/diffgenboxplot.pdf")
graph2tif(file="output/diffgenboxplot.tif")
## 导成AI可以编辑的状态
graph2eps(file="output/diffgenboxplot.eps")

####heatmap热图
### 热图就是表达量数据，行是差异基因，列是样本
### 想一想该如何提取

### 用行名提取数据
rm(list = ls())
## 1.加载表达数据
load(file = "output/exprSet_rmdup.Rdata")
## 2.加载差异列表
load(file = "output/allDiff.Rdata")
library(dplyr)
library(tibble)
diffgene <- allDiff %>% 
  rownames_to_column() %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >1) %>% 
  column_to_rownames()

## 3.用名称提取部分数据用作热图绘制
heatdata <- exprSet[rownames(diffgene),]
## 4.制作一个分组信息用于注释
group <- c(rep("treat",8),rep("con",20)) 
annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(heatdata)

## 加载热图的R包
library(pheatmap)
## 颜色包 viridisLite
library(viridisLite)
## 直接作图
pheatmap(heatdata)

### 经过修饰的图
pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = F,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color = viridis(10, alpha = 1, begin = 0.5, end = 1, direction = 1),#调色
         cellwidth = 60, # 格子宽度
         cellheight = 0.4,# 格子高度
         fontsize = 10 # 字体大小
)

### 可以重新筛选，阈值设大一点
diffgene <- allDiff %>% 
  rownames_to_column() %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >7) %>% 
  column_to_rownames()

## 用名称提取部分数据用作热图绘制
heatdata <- exprSet[rownames(diffgene),]
## 制作一个分组信息用于注释
group <- c(rep("treat",8),rep("con",20)) 
annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(heatdata)

#如果注释出界，可以通过调整格子比例和字体修正
pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color = viridis(10, alpha = 1, begin = 0.5, end = 1, direction = 1),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 10, # 格子宽度
         cellheight = 14,# 格子高度
         fontsize = 10 # 字体大小
)

### 导出图
## 加载R包
library(export)
## 导成PPT可编辑的格式
graph2ppt(file="output/heatmap_modified.pptx")
graph2pdf(file="output/heatmap_modified.pdf")

### 画火山图，订制火山图
rm(list = ls())
##用ggplot2
library(ggplot2)
library(ggrepel)
library(dplyr)
load(file = "output/allDiff.Rdata")
### 无论是什么名称，都改为data
data <- allDiff
### 如果没有gene这一列就需要添加
data$gene <- rownames(data)
## 仔细观察data数据
## 如果是你自己的数据，至少有三列
## logFC，P.Value，gene
ggplot(data=data, aes(x=logFC, y =-log10(P.Value))) +
  ## 三个部分分别画点
  geom_point(data=subset(data,abs(data$logFC) <= 1),aes(size=abs(logFC)),color="black",alpha=0.1) +
  geom_point(data=subset(data,data$P.Value<0.05 & data$logFC > 1),aes(size=abs(logFC)),color="red",alpha=0.2) +
  geom_point(data=subset(data,data$P.Value<0.05 & data$logFC < -1),aes(size=abs(logFC)),color="green",alpha=0.2) +
  ## 画线
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  ## 主题
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
  labs(x="log2 (fold change)",y="-log10 (q-value)")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position='none')+
  ## 标签
  geom_text_repel(data=subset(data, abs(logFC) > 3), aes(label=gene),col="black",alpha = 0.8)


##换一种风格
library(ggplot2)
library(ggrepel)
data <- allDiff
data$gene <- rownames(data)
data$significant <- as.factor(data$adj.P.Val<0.05 & abs(data$logFC) > 1)
data$gene <- rownames(data)
ggplot(data=data, aes(x=logFC, y =-log10(adj.P.Val),color=significant)) +
  geom_point(alpha=0.8, size=1.2,col="black")+
  geom_point(data=subset(data, logFC > 1),alpha=0.8, size=1.2,col="red")+
  geom_point(data=subset(data, logFC < -1),alpha=0.6, size=1.2,col="blue")+
  labs(x="log2 (fold change)",y="-log10 (adj.P.Val)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black")) +
  geom_point(data=subset(data, abs(logFC) >= 3),alpha=0.8, size=3,col="green")+
  geom_text_repel(data=subset(data, abs(logFC) > 3), 
                  aes(label=gene),col="black",alpha = 0.8)
## 加载R包
library(export)
## 导成PPT可编辑的格式
#graph2ppt(file="output/volcano.pptx")
graph2pdf(file="output/volcano.pdf")

ggplot(data=data, aes(x=logFC, y =-log10(P.Value))) +
  ## 三个部分分别画点
  geom_point(data=subset(data,abs(data$logFC) <= 1),aes(size=abs(logFC)),color="black",alpha=0.1) +
  geom_point(data=subset(data,data$P.Value<0.05 & data$logFC > 1),aes(size=abs(logFC)),color="red",alpha=0.2) +
  geom_point(data=subset(data,data$P.Value<0.05 & data$logFC < -1),aes(size=abs(logFC)),color="green",alpha=0.2) +
  ## 画线
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  ## 主题
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
  labs(x="log2 (fold change)",y="-log10 (q-value)")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position='none')+
  ## 标签
  geom_text_repel(data=subset(data, abs(logFC) > 3), aes(label=gene),col="black",alpha = 0.8)


##换一种风格
library(ggplot2)
library(ggrepel)
data <- allDiff
data$gene <- rownames(data)
data$significant <- as.factor(data$adj.P.Val<0.05 & abs(data$logFC) > 1)
data$gene <- rownames(data)
ggplot(data=data, aes(x=logFC, y =-log10(adj.P.Val),color=significant)) +
  geom_point(alpha=0.8, size=1.2,col="black")+
  geom_point(data=subset(data, logFC > 1),alpha=0.8, size=1.2,col="red")+
  geom_point(data=subset(data, logFC < -1),alpha=0.6, size=1.2,col="blue")+
  labs(x="log2 (fold change)",y="-log10 (adj.P.Val)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black")) +
  geom_point(data=subset(data, abs(logFC) >= 3),alpha=0.8, size=3,col="green")+
  geom_text_repel(data=subset(data, abs(logFC) > 3), 
                  aes(label=gene),col="black",alpha = 0.8)

###再换一种风格火山图（可用）
logFCfilter=2
logFCcolor=3
###标记上下调
index = data$adj.P.Val<0.05 & abs(data$logFC) > logFCfilter
data$group <- 0
data$group[index & data$logFC>0]=1
data$group[index & data$logFC<0]=-1
data$group <- factor(data$group,levels = c(1,0,-1), labels=c("Up","NS","Down"))
###正式画图
ggplot(data = data, aes(x = logFC, y = -log10(adj.P.Val), color = group)) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_color_manual(values = c("red", "grey50", "blue4")) +
  labs(x = "log2 (fold change)", y = "-log10 (adj.P.Val)") +  
  theme(plot.title = element_text(hjust = 0.4)) +
  geom_hline(yintercept = -log10(0.05), lty = 4, lwd = 0.6, alpha = 0.8) +
  geom_vline(xintercept = c(-logFCfilter, logFCfilter), lty = 4, lwd = 0.6, alpha = 0.8) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme(legend.position = "top") +  
  geom_point(
    data = subset(data, abs(logFC) >= logFCfilter & adj.P.Val < 0.05),
    alpha = 0.8
  ) +
  ggrepel::geom_text_repel(  
    data = subset(data, abs(logFC) >= logFCfilter & adj.P.Val < 0.05),
    aes(label = gene),
    color = "black",
    alpha = 0.8
  )

## 加载R包
library(export)
## 导成PPT可编辑的格式
graph2ppt(file="output/volcano.pptx")
graph2pdf(file="output/volcano.pdf")

### KEGG分析
rm(list = ls())
library(clusterProfiler)
load(file = "output/diffgene.Rdata")
### 这个分析需要什么数据？
### 获得基因列表
gene <- rownames(diffgene)
#基因名称转换，返回的是数据框
gene = bitr(gene, 
            fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")
head(gene)
#**KEGG分析**
#install.packages("D:\\1SX\\R packages\\KEGG.db_3.2.3.tar.gz",
#repos = NULL, type = "source")
library(KEGG.db)
### organism = 'hsa'
### http://rest.kegg.jp/list/organism
### https://mp.weixin.qq.com/s/PwrdQAkG3pTlwMB6Mj8wXQ
### library(createKEGGdb)
### species <-c("mmu","hsa")
### create_kegg_db(species)

EGG <- enrichKEGG(gene = gene$ENTREZID,
                  organism = 'hsa',
                  pvalueCutoff = 0.05,
                  use_internal_data =F)#F使用在线数据

## 画图
barplot(EGG)
dotplot(EGG)

#barplot(EGG,
#       x = "Count",
#       color = "p.adjust",
#       showCategory = 20,
#       font.size = 12,
#       title = "",
#       mlabel_format = 30)
### KEGG的富集分析比较特殊，他的背后是个网站
KEGG_df <- as.data.frame(EGG)
symboldata <- setReadable(EGG, OrgDb="org.Hs.eg.db", keyType = "ENTREZID")
symboldata  <- as.data.frame(symboldata )

browseKEGG(EGG, 'hsa04110')

save(EGG,file = "output/EGG.Rdata")


### 本节任务：使用pathview 这个R包可视化KEGG的结果

rm(list = ls())
load(file = "output/EGG.Rdata")
KEGG_df <- as.data.frame(EGG)

## 安装这个R包
if(!requireNamespace("pathview")){
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
  BiocManager::install("pathview",update = F,ask = F)
} 

### 这个分析需要什么啊？？
### 很重要的部分，制作geneList
### 什么是geneList

library(clusterProfiler)
load(file = "output/diffgene.Rdata")
gene <- rownames(diffgene)
## 基因名称转换，从symbol 到ENTREZID
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
## 新的数据框
gene_df <- data.frame(logFC=diffgene$logFC,SYMBOL = rownames(diffgene))
## 根据基因名称合并数据
gene_df <- merge(gene_df,gene,by="SYMBOL")

## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$ENTREZID
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

head(geneList)

### pathview可视化
library(pathview)
pathway.id = "hsa04110"
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "hsa",
                   limit      = list(gene=max(abs(geneList)), cpd=1))
## 改变倍数的颜色
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "hsa"#图片对比度更好
                   #limit      = list(gene=max(abs(geneList)), cpd=1)
)
## 改变构图
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "hsa",
                   #limit      = list(gene=max(abs(geneList)), cpd=1),
                   #改变构图,保存为PDF格式
                   kegg.native=F)

### pathview批量画图
### 新建文件夹，
dir.create("output/pathview_out")
### 设置工作目录到想要的地方
getwd()
#remotes::install_github("datapplab/pathview")
#library(pathview)
# 重新运行生成代码

### 然后循环绘图
for (pathway.id in KEGG_df$ID ){
  print(pathway.id)
  pathview(gene.data  = geneList,
           pathway.id = pathway.id,
           species    = "hsa"
  )
}

### 结束后切换回原来的工作目录
getwd()

### 本节任务：基因集富集分析GSEA的原理以及用法
## 解释原理
rm(list = ls())
library(clusterProfiler)
load(file = "output/allDiff.Rdata")
#获得基因列表
gene <- rownames(allDiff)
## 转换
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=allDiff$logFC,
                      SYMBOL = rownames(allDiff))
gene_df <- merge(gene_df,gene,by="SYMBOL")

## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$SYMBOL
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

head(geneList)

### GSEA 变化在于gene set
### 1.hallmarks gene set
## 读入hallmarks gene set，从哪来？
hallmarks <- read.gmt("resource/h.all.v2024.1.Hs.symbols (1).gmt")

### 主程序GSEA
y <- GSEA(geneList,TERM2GENE =hallmarks)

yd <- as.data.frame(y)
### 看整体分布
library(ggplot2)
dotplot(y,showCategory=30,
        split=".sign",
        font.size = 8,
        label_format = 60)+facet_grid(~.sign)

#df <- ggplot2::fortify(y, showCategory = 30, split=".sign")
### https://github.com/YuLab-SMU/enrichplot/blob/master/R/method-fortify.R

### 可以修改标签长度
library(stringi)
library(ggplot2)
dotplot(y,showCategory=12,split=".sign")+
  facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) stri_sub(x,10))

### 选择需要呈现的来作图
library(enrichplot)
gseaplot2(y,"HALLMARK_ALLOGRAFT_REJECTION",color = "red",pvalue_table = T)
gseaplot2(y,10,color = "red",pvalue_table = T)
gseaplot2(y, geneSetID = 1:3)

### cutting edge作图
if(!requireNamespace("ggnewscale",quietly = TRUE)){
  options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  install.packages("ggnewscale",update = F,ask = F)
}
cnetplot(y,showCategory = 4,foldChange = geneList,colorEdge = T)

### 2.kegg 通路
## 读入kegg gene set
kegg <- read.gmt("resource/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt")
y <- GSEA(geneList,TERM2GENE =kegg)
yd <- as.data.frame(y)
### 看整体分布
dotplot(y,showCategory=12,split=".sign")+facet_grid(~.sign)

### 3.转录因子
## 读入转录因子,高能，需要专业知识消化
ENCODE_TF <- read.gmt("resource/ENCODE_TF_ChIP-seq_2015.txt")
tfbs <- GSEA(geneList,TERM2GENE = ENCODE_TF)
tfbsd <- as.data.frame(tfbs)
### 看整体分布 
dotplot(tfbs,showCategory=30,split=".sign")+facet_grid(~.sign)

gseaplot2(tfbs,"E2F4 MEL cell line mm9",color = "red",pvalue_table = T)
### 这些gene set从哪来的
### 1.board 官网
### 2.https://amp.pharm.mssm.edu/Enrichr/#stats
### 该部分所有数据集已经下载，在附赠教程中
### 3.自己自作或者收集

## GZ06_必备技能GSEA，富集分析的神器，量化一切通路。 (课程赠送，微店领取) 
## https://weidian.com/item.html?itemID=3749911422
## 原理文字版
## https://mp.weixin.qq.com/s/FyXMarhQY6r-oymTd491Sw
## 还可以选择火山图来呈现
## https://mp.weixin.qq.com/s/jolWmKoLic5m_M5F1E8K2g
## 其他物种的GSEA，比如老鼠该怎么做呢
### 使用msigdbr包
### https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
## GEO部分，还有个多芯片联合分析，在附赠部分！
## 一定要完成的家庭作业！！！
## https://mp.weixin.qq.com/s/kTFNOrLfwz1pND877V8F5A
## 下游分析特别好的教程！！！！
## https://hbctraining.github.io/DGE_workshop_salmon/lessons/functional_analysis_2019.html
### PPi 网络构建，参考附赠视频
### 以一款大杀器，结束GEO的教程
### https://amp.pharm.mssm.edu/Enrichr/

## GEO教程长期更新的链接是这个:
## https://codingsoeasy.com/archives/geo
