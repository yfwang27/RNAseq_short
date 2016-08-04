library(edgeR)
library(DESeq2)
library(goseq)
library(RColorBrewer)
library(ggplot2)
library(KEGG.db)
library(org.Mm.eg.db)
library(biomaRt)

targets <- read.table("targets.txt",sep="\t",header=TRUE) 

AllCounts<-read.csv(file="AllCounts.csv",row.names = 1)

head(AllCounts)

class(AllCounts)

cData<-data.frame(name=targets$Sample,                                                                          Group=targets$Group,Batch=targets$Batch)

rownames(cData)<-cData[,1]

cData$Group<-relevel(cData$Group,ref="Viv")

dds<-DESeqDataSetFromMatrix(
  countData= AllCounts,colData=cData,
  design=~Group)

dds<-DESeq(dds)

sizeFactors(dds)

head(dispersions(dds))

plotDispEsts(dds)

res<-results(dds) 

# Order results by adjusted p value 
resOrdered<-res[order(res$padj),]

?results      

head(resOrdered)

library(biomaRt)

mart=useMart('ENSEMBL_MART_ENSEMBL',
             dataset='mmusculus_gene_ensembl',
             host="may2012.archive.ensembl.org")

bm<-getBM(attributes=c('ensembl_gene_id',                                                                                         'mgi_symbol'),
          filters ='ensembl_gene_id',
          values=rownames(resOrdered), mart=mart)

# see the first few rows of "bm" object
head(bm)     

resAnnotated <- merge(as.data.frame(resOrdered),bm,by.x=0,by.y=1)
head(resAnnotated)
# change the column name
colnames(resAnnotated)[1]<-"ensembl_gene_id"

resAnnotated<-resAnnotated[order(resAnnotated$pvalue,                                                                                                   decreasing=F),]

# show the result with gene symbol annotation
head(resAnnotated)

write.table(resAnnotated,file="DESeq_result.txt",sep="\t")
write.csv(resAnnotated,file="DESeq_result.csv",  row.names=F)

summary(res)

# How many adjusted p-values were less than 0.05?
sum(res$padj < 0.05, na.rm=TRUE)

plotMA(res, main="DESeq2", ylim=c(-4,4))

idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

plotCounts(dds,gene=which.min(res$padj),intgroup="Group")

rld<-rlog(dds)  

vsd <- varianceStabilizingTransformation(dds)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Group")])

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)


pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


rlogcount <- assay(rld)    
sampleDists <- as.matrix(dist(t(rlogcount)))


png(file="sample_dis_map.png")
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=showcols[colData(dds)$Group], 
          RowSideColors=showcols[colData(dds)$Group],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()


plotPCA(rld, intgroup="Group")
# save the plot

library(ggplot2)
ggsave(file="PCA_plot_version1.png")

library(RColorBrewer)

# Creates nice looking color palettes
showcols <- brewer.pal(8, "Set1")[1:length(unique(colData(dds)$Group))]

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2,label=colData(dds)$name))+
  geom_text(col=showcols[colData(dds)$Group],                                                                                                alpha=0.8,size=4)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
ggsave(file="PCA_plot_version2.png")


res_contrast<-results(dds,contrast=c("Group","Hfd","Viv")) 


summary(res_contrast)


ddsMF<-DESeqDataSetFromMatrix(countData= AllCounts,colData= cData,design=~ Batch + Group)

ddsMF <- DESeq(ddsMF)

resMF <- results(ddsMF)

resMForder<-resMF[order(resMF$padj),]
head(resMForder)

summary(resMForder)

ddsLRT<-DESeqDataSetFromMatrix(countData= AllCounts,colData= cData,design=~ Batch + Group)

# We would like to see the Group effect hence the reduced=~Batch     
ddsLRT <- DESeq(ddsLRT, test="LRT", 
                full=~Batch+ Group, 
                reduced=~Batch)

resLRT<-results(ddsLRT)

resLRTorder<-resLRT[order(resLRT$padj),]

head(resLRTorder)


summary(resLRTorder)

library(KEGG.db)
library(goseq)

# remove the NAs

resdat<- res[complete.cases(res$padj),]

degenes<-as.integer(resdat$padj<0.05)
names(degenes)<-rownames(resdat)

# remove duplicate gene names
degenes<-degenes[match(unique(names(degenes)),                                                                                              names(degenes))]

table(degenes)


pwf=nullp(degenes,genome="mm9",'ensGene', plot.fit=FALSE)

head(pwf)

plotPWF(pwf)

xx <- as.list(KEGGPATHID2NAME)
temp <- cbind(names(xx),unlist(xx))

addKeggTogoseq <- function(JX,temp){
  for(l in 1:nrow(JX)){
    if(JX[l,1] %in% temp[,1]){
      JX[l,"term"] <- temp[temp[,1] %in% JX[l,1],2]
      JX[l,"ontology"] <- "KEGG"
    }
  }
  return(JX)
}

go<-goseq(pwf,genome="mm9",'ensGene', test.cats=c("GO:BP","GO:MF","KEGG"))


head(go)


restemp<-addKeggTogoseq(go,temp)   

head(restemp)

write.table(restemp,file="GO_Kegg_Wallenius.txt", row.names=F,sep="\t")

write.csv(restemp,file="GO_Kegg_Wallenius.csv", row.names=F)

