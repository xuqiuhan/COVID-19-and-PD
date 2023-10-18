

rm(list=ls())

library(pacman)
p_load(edgeR,pheatmap)
setwd("")
inputfile=""
conNum= 86
treatNum= 106
pfilter=0.05   
logfcfilter=0.585  
yifeng<-read.table(inputfile,header = T,sep = "\t",check.names = F)
yifeng=as.matrix(yifeng)
rownames(yifeng)=yifeng[,1]
GeneExp=yifeng[,2:ncol(yifeng)]
yifeng=matrix(as.numeric(as.matrix(GeneExp)),nrow=nrow(GeneExp),dimnames=list(rownames(GeneExp),colnames(GeneExp)))
yifeng=avereps(yifeng)
yifeng=yifeng[rowMeans(yifeng)>1,]
design=c(rep("control",conNum),rep("treat",treatNum))                       
mydesign <- model.matrix(~design)
mydgelist <- DGEList(counts=yifeng,group=design)
mydgelist<- calcNormFactors(mydgelist)
mydgelist<- estimateCommonDisp(mydgelist)
mydgelist <- estimateTagwiseDisp(mydgelist,trend = "movingave")
mytest <- exactTest(mydgelist,pair = c("control","treat"))
all_genes<-topTags(mytest,n=20000000)
all_genes=all_genes$table
iddata<-mydgelist$pseudo.counts
write.table(all_genes,"",sep="\t",quote = F)
diff_genes = all_genes[(all_genes$FDR < pfilter & (all_genes$logFC>=logfcfilter | all_genes$logFC<=(-logfcfilter))),]
write.table(diff_genes, "",sep="\t",quote=F)
up_genes = all_genes[(all_genes$FDR < pfilter & (all_genes$logFC>=logfcfilter)),]
write.table(up_genes, "",sep="\t",quote=F)
down_genes = all_genes[(all_genes$FDR < pfilter & (all_genes$logFC<=(-logfcfilter))),]
write.table(down_genes, "",sep="\t",quote=F)
Normalizegeneexp=rbind(id=colnames(iddata),iddata)
write.table(Normalizegeneexp,"",sep="\t",quote=F,col.names=F)   
diff_genesexp=rbind(id=colnames(iddata),iddata[rownames(diff_genes),])
write.table(diff_genesexp,"",sep="\t",quote=F,col.names=F)
inputheatmap<-iddata[rownames(diff_genes),]
inputheatmap=log2(inputheatmap+0.01)
inputheatmap=inputheatmap[1:20,]
Type=c(rep("Normal",conNum),rep("Parkinson disease",treatNum))   
names(Type)=colnames(inputheatmap)
Type=as.data.frame(Type)
ann_colors = list(
group = c(control="#DC143C", treat="#008080"))
pdf("",10,8)
pheatmap(inputheatmap, 
         scale = "row",
         annotation_col=Type, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(150),cluster_cols =F,
         fontsize = 10,fontsize_row=10,fontsize_col=5,
         show_colnames = F,
         annotation_legend = T,
         annotation_names_col = T,
         annotation_colors =ann_colors[1]
)
dev.off()
nosig<-all_genes[abs(all_genes$logFC)< logfcfilter | all_genes$FDR>=pfilter,]
xmax<-max(all_genes$logFC)
ymax<-max(-log10(all_genes$FDR))
down_genes<- transform(down_genes,FDR=-log10(down_genes$FDR))
up_genes<- transform(up_genes,FDR=-log10(up_genes$FDR))
nosig<- transform(nosig,FDR=-log10(nosig$FDR))
pdf("")
plot(nosig$logFC,nosig$FDR,xlim = c(-xmax,xmax),ylim=c(0,ymax),col="black",
     pch=16,cex=0.9,main = "Volcano",xlab = "logFC",ylab="-log10(FDR)")
points(up_genes$logFC,up_genes$FDR,col="#DC143C",pch=16,cex=0.9)
points(down_genes$logFC,down_genes$FDR,col="#008080",pch=16,cex=0.9)
abline(v=0,lwd=3,lty=2)
dev.off()

