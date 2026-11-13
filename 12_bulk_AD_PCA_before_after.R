# =============================================================================
# Step 1.3 - PCA before and after batch correction (AD)
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

# 安装包 (如果尚未安装)
# install.packages("ggplot2")
# install.packages("ggpubr")

# 引用包
library(ggplot2)
library(ggpubr)

# 设置工作目录 (请确保路径存在)
setwd(file.path(WORK_DIR, "bulk", "AD", "02_PCA"))   # set by 00_config.R

# 定义PCA分析的函数
bioPCA = function(inputFile = null, outFile = null, titleName = null){
  # 读取输入文件,提取数据
  # check.names=F 防止基因名中的特殊字符被修改
  rt = read.table(inputFile, header = T, sep = "\t", check.names = F, row.names = 1)
  data = t(rt)
  
  # 获取分组信息 (GEO数据库研究的id)
  Project = gsub("(.*?)\\_.*", "\\1", rownames(data))    
  
  # PCA分析
  data.pca = prcomp(data)
  pcaPredict = predict(data.pca)
  PCA = data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2], Type = Project)
  
  # 绘制图形
  pdf(file = outFile, width = 5.5, height = 4.25)
  
  # 修改处：添加 palette = "npg" 参数以应用 Nature 配色
  # ggpubr 支持的主题包括: "npg", "aaas", "lancet", "jco", "ucscgb", "d3", "gse", "jama", "nejm"
  p1 = ggscatter(data = PCA, 
                 x = "PC1", y = "PC2", 
                 color = "Type", 
                 shape = "Type", 
                 palette = "npg",          # <--- 这里启用了 NPG 配色
                 ellipse = T, 
                 ellipse.type = "norm", 
                 ellipse.border.remove = F, 
                 ellipse.alpha = 0.1,
                 size = 2, 
                 main = titleName, 
                 legend = "right") +
    theme(plot.margin = unit(rep(1.5, 4), 'lines'), 
          plot.title = element_text(hjust = 0.5))
  
  print(p1)
  dev.off()
  
  # 可选：打印提示
  cat("Plot saved to:", outFile, "\n")
}

# 调用函数, 绘制批次矫正前的图形
bioPCA(inputFile = "merge.preNorm.txt", 
       outFile = "PCA.preNorm.pdf", 
       titleName = "Before batch correction")

# 调用函数, 绘制批次矫正后的图形
# 注意：修正了原代码中 'normalzie' 的拼写错误为 'normalize'
bioPCA(inputFile = "merge.normalize.txt", 
       outFile = "PCA.normalize.pdf", 
       titleName = "After batch correction")
