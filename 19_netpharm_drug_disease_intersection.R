# =============================================================================
# Step 2.3 - Intersect hawthorn targets with the shared AD-AS DEGs
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

# Input : drug target list + interGenes.txt from 01_bulk/04_venn/

#install.packages("ggvenn")


#引用包
library(ggvenn)

setwd(file.path(WORK_DIR, "network"))   # set by 00_config.R

#获取目录下txt结尾的文件
files=list.files(pattern="*.txt$") 
geneList=list()

#读取所有txt文件中的基因信息，保存到geneList
for(inputFile in files){
	if(inputFile=="Disease.txt"){next}
    rt=read.table(inputFile, header=F, sep="\t", check.names=F)      #读取输入文件
    geneNames=unlist(strsplit(as.vector(rt[,1]), " "))       #提取基因名称
    geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
    uniqGene=unique(geneNames)                 #基因取unique
    header=unlist(strsplit(inputFile,"\\.|\\-"))
    geneList[[header[1]]]=uniqGene
}
#绘制venn图
pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList,show_percentage = T,
  stroke_color = "white", stroke_size = 0.5,
  fill_color = c("#E41A1C","#1E90FF","#FF8C00","#31A354FF"),
  set_name_color =c("#E41A1C","#1E90FF","#FF8C00","#31A354FF"),
  set_name_size=6, text_size=4.5)
dev.off()

#输出四个数据库的并集基因
unionGenes=Reduce(union, geneList)
write.table(file="Disease.txt", unionGenes, sep="\t", quote=F, col.names=F, row.names=F)
