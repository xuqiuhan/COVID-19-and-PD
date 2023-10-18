
library(randomForest)
set.seed(123456)
inputFile=""    
setwd("")      
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)
importance=importance(x=rf2)
pdf(file="", width=6.2, height=5.8)
varImpPlot(rf2, main="")
dev.off()
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>2])     
rfGenes=names(rfGenes[1:30])        
write.table(rfGenes, file="", sep="\t", quote=F, col.names=F, row.names=F)
sigExp=t(data[,rfGenes])
sigExpOut=rbind(ID=colnames(sigExp),sigExp)
write.table(sigExpOut, file="", sep="\t", quote=F, col.names=F)




