# =============================================================================
# Step 3.7 - Cell-cell communication analysis with CellChat
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

setwd(file.path(WORK_DIR, "singlecell"))   # set by 00_config.R

# ====================================================
# 6.1 AS数据集CellChat
# ====================================================
library(CellChat)

CellChatDB <- CellChatDB.human

# AS疾病组
AS_d <- subset(AS_obj, Sample == "Atherosclerosis")
cc_AS_d <- createCellChat(AS_d, group.by = "cellType")
cc_AS_d@DB <- CellChatDB
cc_AS_d <- subsetData(cc_AS_d)
cc_AS_d <- identifyOverExpressedGenes(cc_AS_d)
cc_AS_d <- identifyOverExpressedInteractions(cc_AS_d)
cc_AS_d <- computeCommunProb(cc_AS_d, type = "triMean")
cc_AS_d <- filterCommunication(cc_AS_d, min.cells = 10)
cc_AS_d <- computeCommunProbPathway(cc_AS_d)
cc_AS_d <- aggregateNet(cc_AS_d)

# AS对照组
AS_c <- subset(AS_obj, condition == "Control")
cc_AS_c <- createCellChat(AS_c, group.by = "celltype")
cc_AS_c@DB <- CellChatDB
cc_AS_c <- subsetData(cc_AS_c)
cc_AS_c <- identifyOverExpressedGenes(cc_AS_c)
cc_AS_c <- identifyOverExpressedInteractions(cc_AS_c)
cc_AS_c <- computeCommunProb(cc_AS_c, type = "triMean")
cc_AS_c <- filterCommunication(cc_AS_c, min.cells = 10)
cc_AS_c <- computeCommunProbPathway(cc_AS_c)
cc_AS_c <- aggregateNet(cc_AS_c)

# 比较
cc_AS_list <- list(Control = cc_AS_c, AS = cc_AS_d)
cc_AS_merged <- mergeCellChat(cc_AS_list, add.names = names(cc_AS_list))

# 巨噬细胞发出的信号比较
pdf("results/06_AS_cellchat_macro_bubble.pdf", width = 12, height = 8)
netVisual_bubble(cc_AS_merged, sources.use = "Macrophage",
                 comparison = c(1, 2), angle.x = 45,
                 title.name = "Macrophage outgoing signaling: Control vs AS")
dev.off()

# 信号通路强度比较
pdf("results/06_AS_pathway_strength.pdf", width = 10, height = 7)
rankNet(cc_AS_merged, mode = "comparison", comparison = c(1, 2),
        stacked = TRUE, do.stat = TRUE)
dev.off()

# 整体通讯网络
pdf("results/06_AS_network.pdf", width = 12, height = 5)
par(mfrow = c(1, 2))
netVisual_circle(cc_AS_c@net$count, weight.scale = TRUE, title.name = "Control")
netVisual_circle(cc_AS_d@net$count, weight.scale = TRUE, title.name = "AS")
dev.off()

# ====================================================
# 6.2 AD数据集CellChat
# ====================================================

AD_d <- subset(AD_obj, condition == "AD")
cc_AD_d <- createCellChat(AD_d, group.by = "celltype")
cc_AD_d@DB <- CellChatDB
cc_AD_d <- subsetData(cc_AD_d)
cc_AD_d <- identifyOverExpressedGenes(cc_AD_d)
cc_AD_d <- identifyOverExpressedInteractions(cc_AD_d)
cc_AD_d <- computeCommunProb(cc_AD_d, type = "triMean")
cc_AD_d <- filterCommunication(cc_AD_d, min.cells = 10)
cc_AD_d <- computeCommunProbPathway(cc_AD_d)
cc_AD_d <- aggregateNet(cc_AD_d)

AD_c <- subset(AD_obj, condition == "Control")
cc_AD_c <- createCellChat(AD_c, group.by = "celltype")
cc_AD_c@DB <- CellChatDB
cc_AD_c <- subsetData(cc_AD_c)
cc_AD_c <- identifyOverExpressedGenes(cc_AD_c)
cc_AD_c <- identifyOverExpressedInteractions(cc_AD_c)
cc_AD_c <- computeCommunProb(cc_AD_c, type = "triMean")
cc_AD_c <- filterCommunication(cc_AD_c, min.cells = 10)
cc_AD_c <- computeCommunProbPathway(cc_AD_c)
cc_AD_c <- aggregateNet(cc_AD_c)

cc_AD_list <- list(Control = cc_AD_c, AD = cc_AD_d)
cc_AD_merged <- mergeCellChat(cc_AD_list, add.names = names(cc_AD_list))

# 小胶质细胞信号比较
pdf("results/06_AD_cellchat_micro_bubble.pdf", width = 12, height = 8)
netVisual_bubble(cc_AD_merged, sources.use = "Microglia",
                 comparison = c(1, 2), angle.x = 45,
                 title.name = "Microglia outgoing signaling: Control vs AD")
dev.off()

pdf("results/06_AD_pathway_strength.pdf", width = 10, height = 7)
rankNet(cc_AD_merged, mode = "comparison", comparison = c(1, 2),
        stacked = TRUE, do.stat = TRUE)
dev.off()

saveRDS(cc_AS_merged, "results/06_cellchat_AS_merged.rds")
saveRDS(cc_AD_merged, "results/06_cellchat_AD_merged.rds")
