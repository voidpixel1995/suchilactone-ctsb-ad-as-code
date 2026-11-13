# =============================================================================
# Step 1.4 - limma differential expression analysis (AS)
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

# Input : merge.normalize.txt from 01_combat/
# Output: all.txt, diff.xls, diffGeneExp.txt, volcano/heatmap figures

# 安装包 (如果尚未安装)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("limma")
# install.packages("pheatmap")
# install.packages("ggplot2")
# install.packages("ggpubr") # 需要此包来获取 NPG 配色

# 引用包
library(limma)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggpubr) # 引入 ggpubr 以使用 get_palette 和 scale_color_npg

logFCfilter = 0.585           # logFC的过滤条件
adj.P.Val.Filter = 0.05      # 矫正后p值的过滤条件
inputFile = "merge.normalize.txt"      # 表达数据文件
setwd(file.path(WORK_DIR, "bulk", "AS", "03_differential"))   # set by 00_config.R

# --- 数据读取与预处理 ---
rt = read.table(inputFile, header=T, sep="\t", check.names=F)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[,2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data = avereps(data)

# --- 样品分组信息处理 ---
Type = gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data = data[, order(Type)]       # 根据样品的分组信息对样品进行排序
Project = gsub("(.+)\\_(.+)\\_(.+)", "\\1", colnames(data))    # 获得GEO数据研究的id
Type = gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
colnames(data) = gsub("(.+)\\_(.+)\\_(.+)", "\\2", colnames(data))

# --- 差异分析 (limma) ---
design <- model.matrix(~0 + factor(Type))
colnames(design) <- c("Control", "Treat")
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(Treat - Control, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# --- 输出结果 ---
# 输出所有基因的差异情况
allDiff = topTable(fit2, adjust='fdr', number=200000)
allDiffOut = rbind(id=colnames(allDiff), allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

# 输出显著的差异基因
diffSig = allDiff[with(allDiff, (abs(logFC) > logFCfilter & adj.P.Val < adj.P.Val.Filter)), ]
diffSigOut = rbind(id=colnames(diffSig), diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

# 输出差异基因表达量
if(nrow(diffSig) > 0) {
  diffGeneExp = data[row.names(diffSig), ]
  diffGeneExpOut = rbind(id=paste0(colnames(diffGeneExp), "_", Type), diffGeneExp)
  write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)
} else {
  message("No significant differential genes found. Skipping expression output.")
}

# --- 绘制差异基因热图 (修改点 1: NPG 配色) ---
geneNum = 25     # 定义展示基因的数目
if(nrow(diffSig) > 0) {
  diffUp = diffSig[diffSig$logFC > 0, ]
  diffDown = diffSig[diffSig$logFC < 0, ]
  
  geneUp = if(nrow(diffUp) > 0) row.names(diffUp) else character(0)
  geneDown = if(nrow(diffDown) > 0) row.names(diffDown) else character(0)
  
  if(length(geneUp) > geneNum) { geneUp = geneUp[1:geneNum] }
  if(length(geneDown) > geneNum) { geneDown = geneDown[1:geneNum] }
  
  selectedGenes = c(geneUp, geneDown)
  
  if(length(selectedGenes) > 0) {
    hmExp = data[selectedGenes, ]
    
    # 准备注释文件
    names(Type) = colnames(data)
    Type_df = as.data.frame(Type)
    Type_df = cbind(Project, Type_df)
    colnames(Type_df) = c("Project", "Type") # 确保列名清晰
    
    # 【修改处】获取 NPG 配色方案用于热图
    # get_palette 返回十六进制颜色码，我们取前两个主要颜色作为高低表达的代表色，或者生成渐变
    # 这里使用 NPG 调色板生成一个从 低->中->高 的渐变色
    npg_colors = get_palette("npg", 10) 
    # 使用 NPG 风格的颜色生成渐变 (例如取调色板中的蓝色调和红色调部分，或者直接用其风格)
    # 为了保持热图习惯 (蓝-白-红)，我们从 NPG 调色板中选取类似的深色作为端点
    # 或者更简单地，直接使用 NPG 的特定颜色组合。
    # 下面这行代码生成了一个基于 NPG 风格的 50 色渐变
    my_col = colorRampPalette(c(npg_colors[4], "white", npg_colors[1]))(50) 
    # 注：npg 调色板顺序通常是：红, 蓝, 绿, 橙... 索引可能因版本微调，这里假设 1 是红/橙，4 是蓝/紫
    # 如果不确定具体索引，可以使用标准的 NPG 颜色手动指定端点，例如：
    # my_col = colorRampPalette(c("#0073C2", "white", "#EFC000"))(50) # 示例，但最好动态获取
    # 最稳妥的方式：直接使用 ggpubr 提供的颜色向量构建渐变
    # 这里我们采用一种通用策略：用 NPG 的第一个颜色(通常红)和第二个颜色(通常蓝)做两端
    # 确保低表达是冷色，高表达是暖色
    low_col = npg_colors[2] # 通常是蓝色系
    high_col = npg_colors[1] # 通常是红色/橙色系
    my_col = colorRampPalette(c(low_col, "white", high_col))(50)
    
    # 输出图形
    pdf(file="heatmap.pdf", width=10, height=7)
    pheatmap(hmExp, 
             annotation_col = Type_df, 
             color = my_col,          # 使用 NPG 风格渐变
             cluster_cols = F,
             show_colnames = F,
             scale = "row",
             fontsize = 8,
             fontsize_row = 5.5,
             fontsize_col = 8,
             legend = TRUE)
    dev.off()
  } else {
    message("Not enough genes to plot heatmap.")
  }
} else {
  message("No significant genes to plot heatmap.")
}

# --- 绘制火山图 (修改点 2: NPG 配色) ---
# 重新读取或整理数据用于绘图
if(file.exists("all.txt")) {
  rt_vol = read.table("all.txt", header=T, sep="\t", check.names=F, row.names=1) # 读取时去掉第一行作为列名的干扰，如果all.txt第一行是id
  
  # 注意：之前的代码写入 all.txt 时加了 "id" 行，读取时需要处理
  # 如果 read.table 读入后第一列是基因名（因为 row.names=1），则不需要再处理
  # 但原代码写的是: allDiffOut=rbind(id=colnames(allDiff),allDiff) -> 这意味着第一行是列名，第一列是 "id"
  # 让我们重新构建一下用于绘图的数据框，避免文件读取的歧义
  # 直接使用内存中的 allDiff 对象更安全
  
  rt_vol <- allDiff
  rt_vol$Gene <- rownames(allDiff) # 保存基因名
  
  # 定义显著性
  Sig = ifelse((rt_vol$adj.P.Val < adj.P.Val.Filter) & (abs(rt_vol$logFC) > logFCfilter), 
               ifelse(rt_vol$logFC > logFCfilter, "Up", "Down"), 
               "Not")
  
  rt_vol = mutate(rt_vol, Sig = Sig)
  
  # 将 Sig 转换为因子，并指定顺序，确保颜色映射正确
  rt_vol$Sig <- factor(rt_vol$Sig, levels = c("Down", "Not", "Up"))
  
  # 【修改处】获取 NPG 配色
  # ggpubr::get_palette("npg", 3) 会返回 3 个颜色，通常顺序是 红, 蓝, 绿 (取决于具体版本实现)
  # 为了严格控制：Down=蓝, Not=灰, Up=红 (符合常规且美观)
  # 我们从 NPG 调色板取红色和蓝色，灰色保持中性
  npg_pal <- get_palette("npg", 10)
  # 假设 npg_pal[1] 是红色系, npg_pal[2] 是蓝色系
  # 构造自定义颜色向量：对应因子水平 c("Down", "Not", "Up")
  custom_colors <- c("Down" = npg_pal[2], "Not" = "grey", "Up" = npg_pal[1])
  
  p = ggplot(rt_vol, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(col = Sig), alpha = 0.7, size = 1.5) +
    scale_color_manual(values = custom_colors) + # 使用混合了 NPG 色的手动映射
    # 或者如果你只想完全自动应用 NPG 且不关心具体的上下调颜色对应关系，可以用:
    # scale_color_npg() 
    # 但为了科学绘图的严谨性 (上红下蓝)，推荐上面的 manual 方法结合 NPG 色值
    
    labs(title = "Volcano Plot", x = "logFC", y = "-log10(adj.P.Val)") +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.title = element_blank(),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank())
  
  # 输出图形
  pdf(file="vol.pdf", width = 5.5, height = 4.5)
  print(p)
  dev.off()
} else {
  message("all.txt not found, skipping volcano plot.")
}
