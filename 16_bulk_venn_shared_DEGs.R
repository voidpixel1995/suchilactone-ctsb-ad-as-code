# =============================================================================
# Step 1.5 - Shared up- and down-regulated DEGs between AD and AS (Venn)
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

# Input : one *.diff.txt table per disease (copy of 03_differential/*/diff.xls,
#         e.g. AD.diff.txt and AS.diff.txt) in the working directory
# Output: up.pdf, down.pdf, interGenes.Type.txt, interGenes.List.txt

#install.packages("ggvenn")


library(ggvenn)      #引用包
setwd(file.path(WORK_DIR, "bulk", "04_venn"))   # set by 00_config.R

upList=list()        #上调基因的列表
downList=list()      #下调基因的列表
files=list.files(pattern = "\\.diff.txt$")     #提取目录下diff.txt结尾的文件

#对差异分析的结果文件进行循环
for(i in 1:length(files)){
    inputFile=files[i]
	#读取输入文件
    rt=read.table(inputFile, header=T, sep="\t", check.names=F)
    upRT=rt[rt$logFC>0,]       #提取上调基因的数据
    downRT=rt[rt$logFC<0,]     #提取下调基因的数据
    upGenes=as.vector(upRT[,1])          #提取上调基因名称
    downGenes=as.vector(downRT[,1])      #提取下调基因名称
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    upList[[paste0(header[1],"_up")]]=upGenes
    downList[[paste0(header[1],"_down")]]=downGenes
}

#绘制上调基因的venn图
pdf(file="up.pdf", width=6, height=6)
ggvenn(upList, show_percentage = T,
	stroke_color = "white", stroke_size = 0.5,
	fill_color = c("#E41A1C","#1E90FF"),
	set_name_color =c("#E41A1C","#1E90FF"),
	set_name_size=6, text_size=4.5)
dev.off()

#绘制下调基因的venn图
pdf(file="down.pdf", width=6, height=6)
ggvenn(downList, show_percentage = T,
	stroke_color = "white", stroke_size = 0.5,
	fill_color = c("#E41A1C","#1E90FF"),
	set_name_color =c("#E41A1C","#1E90FF"),
	set_name_size=6, text_size=4.5)
dev.off()

#提取交集基因
upInterGenes=Reduce(intersect, upList)        #交集的上调基因
upTab=cbind(upInterGenes, "up")
downInterGenes=Reduce(intersect, downList)    #交集的下调基因
downTab=cbind(downInterGenes, "down")
interTab=data.frame()
if(length(upInterGenes)>0 & length(downInterGenes)>0){
	interTab=rbind(upTab, downTab)
}else if(length(upInterGenes)>0){
	interTab=upTab
}else if(length(downInterGenes)>0){
	interTab=downTab
}
colnames(interTab)=c("Gene", "Type")
#输出交集基因的属性文件
write.table(file="interGenes.Type.txt", interTab, sep="\t", quote=F, col.names=T, row.names=F)
#输出交集基因的列表文件
write.table(file="interGenes.List.txt", interTab[,1], sep="\t", quote=F, col.names=F, row.names=F)
