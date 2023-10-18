


library(limma)             
expFile=""     
diffFile=""       
setwd("")  
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
diffRT=read.table(diffFile, header=T, sep="\t", check.names=F, row.names=1)
diffRT=diffRT[row.names(data),]
Up=data[diffRT[,"logFC"]>0,]
dataDown=data[diffRT[,"logFC"]<0,]
dataUp2=t(apply(dataUp,1,function(x)ifelse(x>median(x),1,0)))
dataDown2=t(apply(dataDown,1,function(x)ifelse(x>median(x),0,1)))

#??bind(dataUp2, dataDown2)
outTab=rbind(id=colnames(outTab), outTab)
write.table(outTab, file="geneScor\t", quote=F, col.names=F)


######Vi