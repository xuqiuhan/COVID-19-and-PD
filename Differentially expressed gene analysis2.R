rm(list=ls())
library(pacman)
p_load(limma,pheatmap,DESeq2)
setwd("")
inputfile=""
conNum= 3   
treatNum= 6   
pfilter=0.05   
logfcfilter=2    
yifeng<-read.table(inputfile,header = T,sep = "\t",check.names = F)
yifeng=as.matrix(yifeng)
rownames(yifeng)=yifeng[,1]
GeneExp=yifeng[,2:ncol(yifeng)]
yifeng=matrix(as.numeric(as.matrix(GeneExp)),nrow=nrow(GeneExp),dimnames=list(rownames(GeneExp),colnames(GeneExp)))
yifeng=avereps(yifeng)
yifeng=yifeng[rowMeans(yifeng)>1,]
yifeng=round(yifeng,0)
coldata <- data.frame(condition = factor(c(rep('control',conNum), rep('treat', treatNum)),
                                         levels = c('control', 'treat')))
desdsf<- DESeqDataSetFromMatrix(countData =yifeng ,colData = coldata ,design= ~condition)
desdsf2 <- DESeq(desdsf, parallel = T)
myres<- results(desdsf2, contrast = c('condition', 'treat', 'control'))
all_genes<- data.frame(myres, stringsAsFactors = FALSE, check.names = FALSE)
all_genes=na.omit(all_genes)
write.table(all_genes,"",sep="\t",quote = F)
diff_genes = all_genes[(all_genes$padj < pfilter & (all_genes$log2FoldChange>=logfcfilter | all_genes$log2FoldChange<=(-logfcfilter))),]
write.table(diff_genes, "",sep="\t",quote=F)
up_genes = all_genes[(all_genes$padj < pfilter & (all_genes$log2FoldChange>=logfcfilter)),]
write.table(up_genes, "",sep="\t",quote=F)
down_genes = all_genes[(all_genes$padj < pfilter & (all_genes$log2FoldChange<=(-logfcfilter))),]
write.table(down_genes, "",sep="\t",quote=F)
Normalizegeneexp=as.data.frame(counts(desdsf2 , normalized=TRUE)) 
Normalizegeneexp1=rbind(id=colnames(Normalizegeneexp),Normalizegeneexp)
write.table(Normalizegeneexp1,"",sep="\t",quote=F,col.names=F)   
diff_genesexp=rbind(id=colnames(Normalizegeneexp1),Normalizegeneexp1[rownames(diff_genes),])
write.table(diff_genesexp,"",sep="\t",quote=F,col.names=F)
inputheatmap<-Normalizegeneexp[rownames(diff_genes),]
inputheatmap=log2(inputheatmap+1)
inputheatmap=inputheatmap[1:20,]
Type=c(rep("Normal",conNum),rep("COVID-19",treatNum))   
names(Type)=colnames(inputheatmap)
Type=as.data.frame(Type)
ann_colors = list(
  group = c(control="#008080", treat="#DC143C"))
pdf("heatmap.pdf",10,8)
pheatmap(inputheatmap, annotation_col=Type, 
         #scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),cluster_cols =F,
         fontsize = 10,fontsize_row=10,fontsize_col=5,
         show_colnames = F,
         annotation_legend = T,
         annotation_names_col = T,
         annotation_colors =ann_colors[1])
dev.off()
nosig<-all_genes[abs(all_genes$log2FoldChange)< logfcfilter | all_genes$padj>=pfilter,]
xmax<-max(all_genes$log2FoldChange)
ymax<-max(-log10(all_genes$padj))
down_genes<- transform(down_genes,padj=-log10(down_genes$padj))
up_genes<- transform(up_genes,padj=-log10(up_genes$padj))
nosig<- transform(nosig,padj=-log10(nosig$padj))
pdf("")
plot(nosig$log2FoldChange,nosig$padj,xlim = c(-xmax,xmax),ylim=c(0,30),col="black",
     pch=16,cex=0.9,main = "Volcano",xlab = "log2FoldChange",ylab="-log10(padj)")
points(up_genes$log2FoldChange,up_genes$padj,col="#DC143C",pch=16,cex=0.9)
points(down_genes$log2FoldChange,down_genes$padj,col="#008080",pch=16,cex=0.9)
abline(v=0,lwd=3,lty=2)
dev.off()

