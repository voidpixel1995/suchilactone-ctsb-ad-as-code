# =============================================================================
# Step 1.10 - SHAP interpretation of the best-performing model
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

# Input : normalize.txt, interGenes.txt, refer.methodLists.txt

#install.packages("caret")
#install.packages("DALEX")
#install.packages("ggplot2")
#install.packages("randomForest")
#install.packages("kernlab")
#install.packages("kernelshap")
#install.packages("pROC")
#install.packages("shapviz")
#install.packages("xgboost")
#install.packages("klaR")


#引用包
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(kernelshap)
library(pROC)
library(shapviz)
library(xgboost)
library(klaR)
#names(getModelInfo())

#set.seed(12345)      #设置种子
inputFile="merge.normalize.txt"      #表达数据文件
geneFile="interGenes.txt"             #基因列表文件
setwd(file.path(WORK_DIR, "bulk", "07_machine_learning"))   # set by 00_config.R

#读取表达数据文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#读取基因列表文件,提取交集特征基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
row.names(data)=gsub("-", "_", row.names(data))

#获取样品分组信息(对照组和实验组)
data=t(data)
group=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))
data=as.data.frame(data)
data$Type=group

#对数据进行分组(训练组和测试组)
inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

#提取用于预测的数据
yTestClass=test$Type
yTest=ifelse(yTestClass=="Control", 0, 1)
test<-test[,-ncol(test)]

#读取机器学习方法的文件
control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
methodRT <- read.table("refer.methodLists.txt", header=T, sep="\t", check.names=F)

#对机器学习方法进行循环,构建机器学习的模型
modelList=list()
AUC=c()
ROCcolor=rainbow(nrow(methodRT))
for(i in 1:nrow(methodRT)){
	name=methodRT[i,"Name"]
	method=methodRT[i,"Method"]
	#构建机器学习的模型
	model=train(Type ~ ., data = train, method=method, trControl = control)
	if(name=="SVM"){
		model=train(Type ~ ., data = train, method=method, prob.model=TRUE, trControl = control)
	}
	#得到每个样品预测的结果
	pred=predict(model, newdata=test, type="prob")
	#绘制ROC曲线
	roc=roc(yTest, as.numeric(pred[,2]))
	AUC=c(AUC, paste0(name, ': ', sprintf("%.03f",roc$auc)))
	modelList[[method]]=as.numeric(roc$auc)
	if(i==1){
		pdf(file="ROC.pdf", width=5.5, height=5)
		plot(roc, print.auc=F, legacy.axes=T, main="", col=ROCcolor[i], lwd=3)
	}else{
		plot(roc, print.auc=F, legacy.axes=T, main="", col=ROCcolor[i], lwd=3, add=T)
	}
}
legend('bottomright', AUC, col=ROCcolor, lwd=3, bty = 'n', cex=0.9)
dev.off()

#根据ROC曲线下的面积获取最优的模型
aucValue=unlist(modelList)
bestMethod=names(which(aucValue==max(aucValue)))
bestMethod
train$Type=ifelse(train$Type=="Control", 0, 1)
model=train(Type ~., data = train, method = bestMethod, trControl=control)

#计算SHAP值
fit=permshap(model, train[,-ncol(train)])
#fit=kernelshap(model, train[,-ncol(train)])      #基因比较多的时候运行这命令
shp <- shapviz(fit, X_pred = train[,-ncol(train)], X=train[,-ncol(train)], interactions=T)
#根据贡献程度对基因进行排序
important=sort(colMeans(abs(shp$S)), decreasing=T)
showVars=names(important)

#绘制柱状图
pdf(file="barplot.pdf", width=6, height=6)
sv_importance(shp, kind="bar", show_numbers=TRUE)+theme_bw()
dev.off()

#绘制蜂群图
pdf(file="bee.pdf", width=7, height=6)
sv_importance(shp, kind = "bee", show_numbers=TRUE)+theme_bw()
dev.off()

#散点图(依赖图)
pdf(file="dependence.pdf", width=9, height=6)
sv_dependence(shp, v = showVars)+theme_bw()
dev.off()

#瀑布图
pdf(file="waterfall.pdf", width=7, height=5)
sv_waterfall(shp, row_id = 1)
dev.off()

#单样品力图
pdf(file="force.pdf", width=9, height=5)
sv_force(shp, row_id = 1)
dev.off()
