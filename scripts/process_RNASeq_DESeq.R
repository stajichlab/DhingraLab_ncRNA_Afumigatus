#!/usr/bin/env Rscript

library(DESeq2)
library(tximport)
library(dplyr)
library(ggplot2)
library(magrittr)
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(fdrtool)
library(geneplotter)
library(EDASeq)
library(tidyverse)
my_pal2 = mypal2 <- colorRampPalette(brewer.pal(6, "YlOrRd"))

countdata <- read.table("results/read_count_A1163.tsv", header=TRUE, row.names=1)

countdata <- countdata[ ,6:ncol(countdata)]
countdata <- as.matrix(countdata)
head(countdata)

samples <- read_csv("samples.csv",col_names=TRUE) %>% arrange(SAMPLE)
exprnames <- samples$SAMPLE
#exprnames <- sub("","",exprnames,perl=FALSE)

# check that experimental columns match in order
all(exprnames %in% colnames(countdata))
all(exprnames == colnames(countdata))
# reorder the columns
#countdata <- countdata[,exprnames]
all(exprnames == colnames(countdata))

# DEseq2 analyses
geno = factor( samples$Genotype)
rep = factor( samples$Replicate)
treatment = factor (samples$Condition)

sampleTable <- data.frame(genotype = geno,
                          replicate = rep,
                          condition = treatment)
rownames(sampleTable) = exprnames

# check this is right order!
#rownames(sampleTable) = colnames(countdata)

ddsAll <- DESeqDataSetFromMatrix(countData = countdata,
                                colData   = sampleTable, 
                                design    = genotype ~ condition )
nrow(ddsAll)
dds <- ddsAll[ rowSums(counts(ddsAll)) > 1, ]
nrow(ddsAll)
ddsAll <- estimateSizeFactors(ddsAll)
ddsAll <- estimateDispersions(ddsAll)

rld <- rlog(ddsAll, blind=FALSE)
vsd <- vst(ddsAll, blind=FALSE)
df <- bind_rows(
  as.data.frame(log2(counts(ddsAll, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as.data.frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as.data.frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

pdf("plots/RNASeq_sumplots_All.pdf")

plotDispEsts(ddsAll)

multidensity( counts(ddsAll, normalized = T),
              xlab="mean counts", xlim=c(0, 1000))
multiecdf( counts(ddsAll, normalized = T),
           xlab="mean counts", xlim=c(0, 1000))

MA.idx = t(combn(1:4, 2))
for( i in  seq_along( MA.idx[,1])){
  MDPlot(counts(ddsAll, normalized = T),
         c(MA.idx[i,1],MA.idx[i,2]),
         main = paste( colnames(ddsAll)[MA.idx[i,1]], " vs ",
                       colnames(ddsAll)[MA.idx[i,2]] ), ylim = c(-3,3))
}
select <- order(rowMeans(counts(ddsAll,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df2 <- as.data.frame(colData(ddsAll)[,c("condition","genotype")])
rownames(df2) = colnames(countdata)
colnames(df2) = c("Treatment", "Genotype")

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD Top Expression")

topVar <- head(order(rowVars(assay(vsd)),
                     decreasing=TRUE),60)
mat  <- assay(vsd)[ topVar, ]
pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD Most different")

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD Top Expression")

pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD Top Expression")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$genotype, vsd$condition,sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pcaData <- plotPCA(vsd, intgroup=c("genotype","condition","replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

dev.off()

p<- ggplot(pcaData, aes(PC1, PC2, color=genotype,shape=treatment,label=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()

ggsave("plots/PCA_expression.pdf",p)

# only look at one treatment (growth condition) at a time contrast genotypes
calcPlotPairwiseCondition <- function(c) {
  c=unique(sampleTable$condition)[1]
  sampleTablePair = sampleTable %>% dplyr::filter(condition == c)
  print(sampleTablePair)
  pdf(sprintf("plots/RNASeq_%s_contrastGenotype.pdf",c))
  genes = rownames(countdata)
  countdataPair <- as.data.frame(as_tibble(countdata) %>% dplyr::select(contains(sprintf('_%s_',c))))
  
  rownames(countdataPair) <- genes
  
  dds <- DESeqDataSetFromMatrix(countData = countdataPair,
                                colData   = sampleTablePair, 
                                design    = ~ genotype )

  nrow(dds)
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  nrow(dds)

  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)

  vsd <- vst(dds, blind=FALSE)
  rld <- rlog(dds, blind=FALSE)
  head(assay(vsd), 3)
  head(assay(rld), 3)

  df <- bind_rows(
      as.data.frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
        mutate(transformation = "log2(x + 1)"),
      as.data.frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
      as.data.frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst")
      )
    

  colnames(df)[1:2] <- c("x", "y")
  
  select <- order(rowMeans(counts(ddsAll,normalized=TRUE)),
                  decreasing=TRUE)[1:50]

  MA.idx = t(combn(1:4, 2))
  for( i in  seq_along( MA.idx[,1])){ 
    MDPlot(counts(dds, normalized = T), 
          c(MA.idx[i,1],MA.idx[i,2]), 
          main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                         colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3))
  }
  
  ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation)
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:50]
  
  df2 <- as.data.frame(colData(dds)[,c("condition","genotype")])
  rownames(df2) = colnames(countdataPair)
  colnames(df2) = c("Treatment", "Genotype")
  pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
           fontsize_row = 7,fontsize_col = 7,
           cluster_cols=FALSE, annotation_col=df2,main="VSD Top Expression")
  
  topVar <- head(order(rowVars(assay(vsd)),
                       decreasing=TRUE),60)
  mat  <- assay(vsd)[ topVar, ]
  
  pheatmap(mat, show_rownames=TRUE,
           fontsize_row = 7,fontsize_col = 7,
           cluster_cols=FALSE, annotation_col=df2,main="VSD Most different")
  
  pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
           fontsize_row = 7,fontsize_col = 7,
           cluster_cols=FALSE, annotation_col=df2,main="RLD Top Expression")
  
  pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
           fontsize_row = 7,fontsize_col = 7,
           cluster_cols=FALSE, annotation_col=df2,main="RLD Top Expression")
  
  topVar <- head(order(rowVars(assay(rld)),
                       decreasing=TRUE),60)
  mat  <- assay(rld)[ topVar, ]
  mat  <- mat - rowMeans(mat)
  pheatmap(mat, show_rownames=TRUE,
           fontsize_row = 7,fontsize_col = 7,
           cluster_cols=FALSE, annotation_col=df2,main="RLD Most different")
  
  norm.counts <- counts(dds, normalized=TRUE)
  log.norm.counts <- log2(norm.counts + 1)
  
  topVarGenes <- order(-rowVars(log.norm.counts)[0:100])
  mat<-log.norm.counts[topVarGenes,]
  mat<-mat -rowMeans(mat)
  
  pheatmap(mat,method="complete",main = "TopVar normalized", show_rownames = F,
          annotation_legend = FALSE, legend=T, cluster_rows=TRUE, cexRow=0.55 )
  
  topVarGenes1 <- order(-rowVars(assay(rld)))[0:1000]
  mat1 <- assay(rld)[ topVarGenes1, ]
  mat1<- mat1 - rowMeans(mat1)
  
  pheatmap(mat1, method="complete",
           main = "Unsupervised 1000 genes ",
           show_rownames = F,annotation_legend = FALSE, 
           legend=T, cluster_cols=TRUE)
  
  dds <- DESeq(dds)
  res <- results(dds)
  res
  print(resultsNames(dds))
  maxCompare = length(resultsNames(dds))
  for (coefcompare in resultsNames(dds)[2:maxCompare] ) {
    print(coefcompare)
    
    resLFC <- lfcShrink(dds, coef=coefcompare, type="apeglm")
    summary(resLFC)
    res05 <- results(dds,alpha=0.05)
    summary(res05)
#    sum(res$padj < 0.01, na.rm=TRUE)
    
    resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) >= 1 & baseMean >= 5)
    resSig <- resSig[order(resSig$pvalue),]
    pdf(sprintf("plots/plotMA_%s_%s.pdf",c,coefcompare))
    plotMA(res05, ylim=c(-3,3),main=sprintf("%s",coefcompare))
    plotMA(resLFC, ylim=c(-3,3),main=sprintf("%s",coefcompare))
    plotMA(resSig, ylim=c(-3,3),main=sprintf("signif only %s",coefcompare))
    dev.off()
    write.csv(resSig,sprintf("reports/%s_%s.csv",c,coefcompare))
    write.csv(fpm(dds),sprintf("reports/%s_allGeno_FPM.csv",c))
  }
}

# only look at one genotype (strain) at a time contrast conditions
calcPlotPairwiseGenotype <- function(g) {
  sampleTablePair = sampleTable %>% dplyr::filter(genotype == g)
  print(sampleTablePair)
  pdf(sprintf("plots/RNASeq_%s_contrastCondition.pdf",g))
  genes = rownames(countdata)
  countdataPair <- as.data.frame(as_tibble(countdata) %>% dplyr::select(starts_with(sprintf('%s_',g))))
  
  rownames(countdataPair) <- genes
  
  dds <- DESeqDataSetFromMatrix(countData = countdataPair,
                                colData   = sampleTablePair, 
                                design    = ~ condition )

  nrow(dds)
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  nrow(dds)

  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)

  vsd <- vst(dds, blind=FALSE)
  rld <- rlog(dds, blind=FALSE)
  head(assay(vsd), 3)
  head(assay(rld), 3)

  df <- bind_rows(
      as.data.frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
        mutate(transformation = "log2(x + 1)"),
      as.data.frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
      as.data.frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst")
      )
    

  colnames(df)[1:2] <- c("x", "y")
  
  select <- order(rowMeans(counts(ddsAll,normalized=TRUE)),
                  decreasing=TRUE)[1:50]

  MA.idx = t(combn(1:4, 2))
  for( i in  seq_along( MA.idx[,1])){ 
    MDPlot(counts(dds, normalized = T), 
          c(MA.idx[i,1],MA.idx[i,2]), 
          main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                         colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3))
  }
  
  ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation)
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:50]
  
  #datCollapsed <- collapseReplicates(dds, groupby=dds$genotype,run=dds$replicate,renameCols=TRUE)
  
  df2 <- as.data.frame(colData(dds)[,c("condition","genotype")])
  rownames(df2) = colnames(countdataPair)
  colnames(df2) = c("Treatment", "Genotype")
  pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
           fontsize_row = 7,fontsize_col = 7,
           cluster_cols=FALSE, annotation_col=df2,main="VSD Top Expression")
  
  topVar <- head(order(rowVars(assay(vsd)),
                       decreasing=TRUE),60)
  mat  <- assay(vsd)[ topVar, ]
  
  pheatmap(mat, show_rownames=TRUE,
           fontsize_row = 7,fontsize_col = 7,
           cluster_cols=FALSE, annotation_col=df2,main="VSD Most different")
  
  pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
           fontsize_row = 7,fontsize_col = 7,
           cluster_cols=FALSE, annotation_col=df2,main="RLD Top Expression")
  
  pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
           fontsize_row = 7,fontsize_col = 7,
           cluster_cols=FALSE, annotation_col=df2,main="RLD Top Expression")
  
  topVar <- head(order(rowVars(assay(rld)),
                       decreasing=TRUE),60)
  mat  <- assay(rld)[ topVar, ]
  mat  <- mat - rowMeans(mat)
  pheatmap(mat, show_rownames=TRUE,
           fontsize_row = 7,fontsize_col = 7,
           cluster_cols=FALSE, annotation_col=df2,main="RLD Most different")
  
 
  # FIX ME HERE
  norm.counts <- counts(dds, normalized=TRUE)
  log.norm.counts <- log2(norm.counts + 1)
  
  topVarGenes <- order(-rowVars(log.norm.counts)[0:100])
  mat<-log.norm.counts[topVarGenes,]
  mat<-mat -rowMeans(mat)
  
  pheatmap(mat,method="complete",main = "TopVar normalized", show_rownames = F,
          annotation_legend = FALSE, legend=T, cluster_rows=TRUE, cexRow=0.55 )
  
  topVarGenes1 <- order(-rowVars(assay(rld)))[0:1000]
  mat1 <- assay(rld)[ topVarGenes1, ]
  mat1<- mat1 - rowMeans(mat1)
  
  pheatmap(mat1, method="complete",
           main = "Unsupervised 1000 genes ",
           show_rownames = F,annotation_legend = FALSE, 
           legend=T, cluster_cols=TRUE)
  dev.off()
  
  dds <- DESeq(dds)
  res <- results(dds)
  res    
  #res <- results(dds, contrast=c("treatment","AF293","CEA10"))
  #res
  
  resultsNames(dds)
  resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
  resLFC
  summary(resLFC)
  summary(res)
  res05 <- results(dds,alpha=0.05)
  summary(res05)
  sum(res$padj < 0.01, na.rm=TRUE)

  pdf(sprintf("plots/plotMA_%s_%s.pdf",g,resultsNames(dds)[2]))
  plotMA(res, ylim=c(-2,2),main=sprintf("%s",resultsNames(dds)[2]))
  plotMA(resLFC, ylim=c(-2,2),main=sprintf("%s",resultsNames(dds)[2]))
  
  resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) >= 1 & baseMean >= 5)
  resSig <- resSig[order(resSig$pvalue),]
  write.csv(resSig,sprintf("reports/%s_%s.csv",g,resultsNames(dds)[2]))
  write.csv(fpm(dds),sprintf("reports/%s_allConditions_FPM.csv",g))
  dev.off()
}
lapply(unique(sampleTable$condition),calcPlotPairwiseCondition)
lapply(unique(sampleTable$genotype),calcPlotPairwiseGenotype)



