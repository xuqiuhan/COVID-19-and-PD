library(limma)
library(pheatmap)
inputFile=""
setwd("")    
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="", width=8, height=6)
pheatmap(data, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =T,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
dev.off()


