# =============================================================================
# Step 2.2 - Merge the compound-target tables of the three prediction platforms
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

# SwissTargetPrediction, SEA and TCMSP target lists are merged and
# standardised to HGNC symbols; duplicated edges are removed.

probFilter=0.05      #Probability的过滤条件
setwd(file.path(WORK_DIR, "network"))   # set by 00_config.R

#获取目录下csv结尾的文件
files=list.files(pattern="*.csv$")
files=setdiff(files, "swissadme.filter.csv")

#读取SwissADME数据库过滤后的结果文件
admeRT=read.csv("swissadme.filter.csv", header=T, sep=",", check.names=F, row.names=1)

#读取所有csv文件中的基因信息
outTab=data.frame()
for(inputFile in files){
	rt=read.csv(inputFile, header=T, sep=",",check.names=FALSE)
	rt=subset(rt, `Probability*` > probFilter)
	geneNames=unlist(strsplit(as.vector(rt[,2]), " "))       #提取基因名称
	geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
	uniqGene=unique(geneNames)                 #基因取unique
	header=gsub(".csv", "", inputFile)
	Ingredient=admeRT[header, "Ingredient name"]
	if(length(uniqGene)>0){
		outTab=rbind(outTab, cbind(id=header, Ingredient=Ingredient, Gene=uniqGene))
	}
}

#输出过滤后的结果
write.table(outTab, file="Drug.txt", sep="\t", quote=F, row.names=F)
