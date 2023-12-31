library(pROC)                  
inputFile=""    
setwd("")  
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))
y=ifelse(y=="con", 0, 1)
roc1=roc(y, as.numeric(rt[,2]))
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)
pdf(file="", width=5, height=5)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main="Train group")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()




