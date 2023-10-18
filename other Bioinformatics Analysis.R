
library(ReactomePA) 
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
setwd("") 
genelist_input <- fread(file="", header = T, sep='\t', data.table = F)
genename <- as.character(genelist_input[,1])
gene_map <- select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID")) 
non_duplicates_idx <- which(duplicated(gene_map$SYMBOL) == FALSE)
gene_map <- gene_map[non_duplicates_idx, ]
colnames(gene_map)[1]<-"Gene"
temp<-inner_join(gene_map,genelist_input,by = "Gene")
temp<-temp[,-1]
temp<-na.omit(temp)
temp$logFC<-sort(temp$logFC,decreasing = T)
geneList = temp[,2]
names(geneList) = as.character(temp[,1])
geneList
KEGG_gseresult <- gseKEGG(geneList, pvalueCutoff=0.05) 
go_results<-as.data.frame(Go_gseresult)
kegg_results<-as.data.frame(KEGG_gseresult)
write.csv (go_results, file ="Go_gseresult.csv")
write.csv (kegg_results, file ="KEGG_gseresult.csv")
gseaplot2(Go_gseresult, 1:5, title = "", pvalue_table = FALSE)
