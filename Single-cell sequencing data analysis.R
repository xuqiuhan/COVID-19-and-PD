setwd("")
rm(list = ls())
source("/mnt/8t/gse/代码库/单细胞/library.R")
DA_scRNA=qs::qread("DA_scRNA_recluster.qs")
DA_scRNA$group
DA_scRNA$subcluster=DA_scRNA$RNA_snn_res.0.3
DA_scRNA$subcluster=gsub("0","10",DA_scRNA$subcluster)
DA_scRNA$subcluster=paste0("Cluster ",DA_scRNA$subcluster)
DimPlot(DA_scRNA,group.by = "group",label = T)+ggtitle("UMAP")
ggsave("1.UMAP.pdf",width = 10,height = 10)
Idents(DA_scRNA)=DA_scRNA$group
plan(multisession,workers=10)
ALL_MARKER=FindAllMarkers(DA_scRNA,logfc.threshold = 0.1,min.pct = 0.1,only.pos = T,test.use = "LR")
plan(multisession,workers=1)
top5=ALL_MARKER %>% group_by(cluster) %>% top_n(5,avg_log2FC)
pdf("2.AverageHeatmap.pdf",width = 10,height = 15)
scRNAtoolVis::AverageHeatmap(DA_scRNA,markerGene = top5$gene)
dev.off()
VlnPlot(DA_scRNA,features = c("","","","","","","","","","","",""),
        pt.size = 0,ncol = 4)
ggsave("3.VlnPlot.pdf",width = 20,height = 15)
Cellratio <- prop.table(table(Idents(DA_scRNA), DA_scRNA$disease__ontology_label), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("Cluster","Var2","Freq")
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#DEB887")
library(ggplot2)
p <- ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Cluster),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Cluster',y = 'Cell Fraction')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+RotatedAxis()+ggplot2::coord_flip()
ggsave(plot=p,filename = "4.1.pdf",width = 10,height = 6)
DA_scRNA$orig.ident=as.character.factor(DA_scRNA$orig.ident)
Cellratio <- prop.table(table(Idents(DA_scRNA), DA_scRNA$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("Cluster","Var2","Freq")
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#DEB887")
library(ggplot2)
p <- ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Cluster),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Cluster',y = 'Cell Fraction')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+RotatedAxis()+ggplot2::coord_flip()
ggsave(plot=p,filename = "4.2.pdf",width = 10,height = 15)
DA_scRNA$subcluster_pc=paste0(DA_scRNA$group,"_",DA_scRNA$disease__ontology_label)
VlnPlot(DA_scRNA,features = "",group.by = "new_celltypes")
ggsave(filename = "5.1.pdf",width = 20,height = 10)














    
    
    
    

  
  
  
  
  
  
  
 
  
 
    
    
 
    
    
    
 