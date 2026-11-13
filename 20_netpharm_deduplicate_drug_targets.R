# =============================================================================
# Step 2.4 - Remove duplicated compound-target pairs before network construction
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

library(ggvenn)

drugName = "Pueraria"
diseaseName = "Alzheimer"
setwd(file.path(WORK_DIR, "network"))   # set by 00_config.R

# 读取药物基因（带格式清洗）
drugRT = read.table("Drug.txt", header=TRUE, sep="\t", comment.char="", quote="")
drugGenes = unique(toupper(trimws(as.vector(drugRT[,3]))))

# 读取疾病基因（带格式清洗）
disRT = read.table("Disease.txt", header=FALSE, sep="\t", comment.char="", quote="")
diseaseGenes = unique(toupper(trimws(as.vector(disRT[,1]))))

# 打印实际基因数量验证
cat("实际药物基因:", length(drugGenes), "\n",
    "实际疾病基因:", length(diseaseGenes), "\n",
    "实际交集:", length(intersect(drugGenes, diseaseGenes)), "\n")

# 创建绘图数据
geneList = list(Pueraria = drugGenes, Alzheimer = diseaseGenes)

# 绘制韦恩图
pdf("venn.pdf", 6, 6)
ggvenn(geneList, 
       show_percentage = TRUE,
       stroke_size = 0.5,
       fill_color = c("#E41A1C","#1E90FF"),
       set_name_color = c("#E41A1C","#1E90FF"))
dev.off()
