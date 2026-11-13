# =============================================================================
# Step 2.1 - SwissADME pharmacokinetic filtering of the hawthorn compounds
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

# Input : one SwissADME CSV export per compound (data/reference/swissadme/)
# Output: swissadme.filter.csv

inputFile="swissadme.csv"       #输入文件
setwd(file.path(WORK_DIR, "network"))   # set by 00_config.R

data=read.csv(inputFile, header=T, sep=",", check.names=FALSE)    #读取输入文件

#对药物成分进行过滤
outTab=subset(data, MW<500 &                        #分子量小于500
			 `#H-bond donors` < 5 &                #氢键供体数目小于5
			 `#H-bond acceptors` < 10 &           #氢键受体数目小于10
			 `#Rotatable bonds` < 10 &            #可旋转键的数量不超过10个
			 (`iLOGP` > -2 & `iLOGP` < 5) &       #脂水分配系数小于5
			 (`XLOGP3` > -2 & `XLOGP3` < 5) &   
			 (`Consensus Log P` > -2 & `Consensus Log P` < 5) )
#保留胃肠道吸收高的药物成分
outTab=subset(outTab, `GI absorption`=="High")
#保留具有血脑屏障通透性的药物成分
#outTab=subset(outTab, `BBB permeant`=="Yes")

#输出过滤后的结果
write.csv(outTab, file="swissadme.filter.csv", row.names=F)
