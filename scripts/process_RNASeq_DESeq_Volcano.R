#!/usr/bin/env Rscript

# Here is a list of volcano plots:

# 1. WT vs de
# 2. Wt vs comp
# 3. Wt vs OE
# Al in GMM (No Drug)

# 4. WT GMM vs WT drug
# 5. delta GMM vs delta drug
# 5. WT drug vs del drug
# 6. Wt drug vs comp drug
# 7 WT drug vs OE drug

# Colored dots to focus on:
#  
#  Green ones are significant in RNAseq WT GMM vs WT drug
#First 7 in yellow are significant WT GMM vs WT drug
# Our focus is on afub_000050 and afub_000060. They are negatively regulated by delta 182 so they go up in WT GMM vs WT drug and go up in WT drug vs delta drug

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
library(EnhancedVolcano)
my_pal2 = mypal2 <- colorRampPalette(brewer.pal(6, "YlOrRd"))

pdf("plots/misc_plots.pdf")
countdata <- read.table("results/read_count_A1163.tsv", header=TRUE, row.names=1)

countdata <- countdata[ ,6:ncol(countdata)]
countdata <- as.matrix(countdata)
head(countdata)

samplesFile <- read_csv("samples.csv",col_names=TRUE) %>% arrange(SAMPLE)
exprnames <- samplesFile$SAMPLE
exprnamesGMM <-  (samplesFile %>% dplyr::filter(Condition == "GMM"))$SAMPLE
#exprnames <- sub("","",exprnames,perl=FALSE)

# check that experimental columns match in order
all(exprnames %in% colnames(countdata))
all(exprnames == colnames(countdata))
# reorder the columns
#countdata <- countdata[,exprnames]
all(exprnames == colnames(countdata))

# DEseq2 analyses
geno = factor( samplesFile$Genotype)
rep = factor( samplesFile$Replicate)
treatment = factor (samplesFile$Condition)

sampleTable <- data.frame(genotype = geno,
                          replicate = rep,
                          condition = treatment)
# GMM only
sampleTableGMM = subset(sampleTable,condition == "GMM")
rownames(sampleTableGMM) = exprnamesGMM

# check this is right order!
#rownames(sampleTable) = colnames(countdata)

countdataGMM = countdata[,exprnamesGMM]

ddsGMM <- DESeqDataSetFromMatrix(countData = countdataGMM,
                                 colData   = sampleTableGMM, 
                                 design    = ~ genotype)
nrow(ddsGMM)
ddsGMM <- ddsGMM[ rowSums(counts(ddsGMM)) > 1, ]
nrow(ddsGMM)
ddsGMM <- estimateSizeFactors(ddsGMM)
ddsGMM <- estimateDispersions(ddsGMM)
rld <- rlogTransformation(ddsGMM)
ddsGMM <- DESeq(ddsGMM)
sizeFactors(ddsGMM)

res <- results(ddsGMM)

res <- results(ddsGMM,
               contrast = c('genotype', 'WT','Del'))
res <- lfcShrink(ddsGMM,
                 contrast = c('genotype','WT','Del'), res=res, type = 'ashr')
res_tbl <- as_tibble(res, rownames = "ENSEMBL")

p<-EnhancedVolcano(res,
                   FCcutoff = 0.5,
                   pointSize = 2.0,
                   labSize = 6.0,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  title = 'WT versus Delta',
                  colAlpha = 1,
                  y = 'pvalue'
                  )
p
ggsave("plots/Volcano.GMM.WT-vs-Delta.pdf",p,width=16,height=12)

best_genes <- res_tbl %>%
  arrange(padj)  %>%
  head(9)
write_tsv(res_tbl,"results/GMM.tsv")

p <- as_tibble(counts(ddsGMM[best_genes$ENSEMBL, ], normalize = TRUE),
               rownames = 'ENSEMBL') %>% 
  pivot_longer(names_to = "sample", values_to = "counts", -ENSEMBL) %>%  
  left_join(as_tibble(sampleTableGMM, rownames = "sample")) %>% 
  ggplot(aes(x = sample, y = counts, fill = genotype)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ ENSEMBL, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))

p
ggsave("plots/barPlot.GMM.top_9.pdf",p,width=16,height=12)

p <- as_tibble(counts(ddsGMM[best_genes$ENSEMBL, ], normalize = TRUE),
          rownames = 'ENSEMBL') %>% 
  pivot_longer(names_to = "sample", values_to = "counts", -ENSEMBL) %>%  
  left_join(as_tibble(sampleTableGMM, rownames = "sample")) %>% 
  filter(genotype %in% c('WT', 'Del')) %>%
  ggplot(aes(x = sample, y = counts, fill = genotype)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ ENSEMBL, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))

p
ggsave("plots/barPlot.GMM.top_9_WT_v_Delta.pdf",p,width=16,height=12)


res <- results(ddsGMM,
               contrast = c('genotype', 'WT','Comp'))
res <- lfcShrink(ddsGMM,
                 contrast = c('genotype','WT','Comp'), res=res, type = 'ashr')
res_tbl <- as_tibble(res, rownames = "ENSEMBL")
p<-EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                title = 'WT versus Comp',
                y = 'pvalue')

ggsave("plots/Volcano.GMM.WT-vs-Comp.pdf",p,width=6,height=8)
p <- as_tibble(counts(ddsGMM[best_genes$ENSEMBL, ], normalize = TRUE),
               rownames = 'ENSEMBL') %>% 
  pivot_longer(names_to = "sample", values_to = "counts", -ENSEMBL) %>%  
  left_join(as_tibble(sampleTableGMM, rownames = "sample")) %>% 
  filter(genotype %in% c('WT', 'Comp')) %>%
  ggplot(aes(x = sample, y = counts, fill = genotype)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ ENSEMBL, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))

p
ggsave("plots/barPlot.GMM.top_9_WT_v_Comp.pdf",p,width=16,height=12)


res <- results(ddsGMM,
               contrast = c('genotype', 'WT','OE'))
res <- lfcShrink(ddsGMM,
                 contrast = c('genotype','WT','OE'), res=res, type = 'ashr')
res_tbl <- as_tibble(res, rownames = "ENSEMBL")
p <- EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                title = 'WT versus OverExpression',
                y = 'pvalue')

ggsave("plots/Volcano.GMM.WT-vs-OE.pdf",p,width=6,height=8)
p <- as_tibble(counts(ddsGMM[best_genes$ENSEMBL, ], normalize = TRUE),
               rownames = 'ENSEMBL') %>% 
  pivot_longer(names_to = "sample", values_to = "counts", -ENSEMBL) %>%  
  left_join(as_tibble(sampleTableGMM, rownames = "sample")) %>% 
  filter(genotype %in% c('WT', 'OE')) %>%
  ggplot(aes(x = sample, y = counts, fill = genotype)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ ENSEMBL, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))
p
ggsave("plots/barPlot.GMM.top_9_WT_v_OE.pdf",p,width=16,height=12)
# 1. WT vs de
# 2. Wt vs comp
# 3. Wt vs OE
# Al in GMM (No Drug)

#====
# 4. WT GMM vs WT drug

sampleTableWT = subset(sampleTable,genotype == "WT")
exprnamesWT <-  (samplesFile %>% dplyr::filter(Genotype == "WT"))$SAMPLE

rownames(sampleTableWT) = exprnamesWT

# check this is right order!
#rownames(sampleTable) = colnames(countdata)

countdataWT = countdata[,exprnamesWT]
ddsWT <- DESeqDataSetFromMatrix(countData = countdataWT,
                                 colData   = sampleTableWT, 
                                 design    = ~ condition)
nrow(ddsWT)
ddsWT <- ddsWT[ rowSums(counts(ddsWT)) > 1, ]
nrow(ddsWT)
ddsWT <- estimateSizeFactors(ddsWT)
ddsWT <- estimateDispersions(ddsWT)
sizeFactors(ddsWT)

ddsWT <- DESeq(ddsWT)
res <- results(ddsWT)

# 4. WT GMM vs WT drug

res <- results(ddsWT,
               contrast = c('condition', 'GMM','Azole'))

res <- lfcShrink(ddsWT,
                 contrast = c('condition','GMM','Azole'), res=res, type = 'ashr')
res_tbl <- as_tibble(res, rownames = "ENSEMBL")
p<-EnhancedVolcano(res,
                   lab = rownames(res),
                   x = 'log2FoldChange',
                   title = 'WT GMM versus Azole',
                   y = 'pvalue')

ggsave("plots/Volcano.WT_GMM-vs-Azole.pdf",p,width=6,height=8)
res_tbl <- as_tibble(res, rownames = "ENSEMBL")
best_genes <- res_tbl %>%
  arrange(padj)  %>%
  head(9)
write_tsv(res_tbl,"results/WT.tsv")

p <- as_tibble(counts(ddsWT[best_genes$ENSEMBL, ], normalize = TRUE),
               rownames = 'ENSEMBL') %>% 
  pivot_longer(names_to = "sample", values_to = "counts", -ENSEMBL) %>%  
  left_join(as_tibble(sampleTableWT, rownames = "sample")) %>% 
  filter(condition %in% c('GMM', 'Azole')) %>%
  ggplot(aes(x = sample, y = counts, fill = genotype)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ ENSEMBL, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))
p
ggsave("plots/barPlot.WT.top_9_GMM_v_Azole.pdf",p,width=16,height=12)
#==
# 5. delta GMM vs delta drug

sampleTableDel = subset(sampleTable,genotype == "Del")
exprnamesDel <-  (samplesFile %>% dplyr::filter(Genotype == "Del"))$SAMPLE

rownames(sampleTableDel) = exprnamesDel

countdataDel = countdata[,exprnamesDel]
ddsDel <- DESeqDataSetFromMatrix(countData = countdataDel,
                                colData   = sampleTableDel, 
                                design    = ~ condition)
nrow(ddsDel)
ddsDel <- ddsDel[ rowSums(counts(ddsDel)) > 1, ]
nrow(ddsDel)
ddsDel <- estimateSizeFactors(ddsDel)
ddsDel <- estimateDispersions(ddsDel)
sizeFactors(ddsDel)

ddsDel <- DESeq(ddsDel)
res <- results(ddsDel)

# 5. delta GMM vs delta drug

res <- results(ddsDel,
               contrast = c('condition', 'GMM','Azole'))

res <- lfcShrink(ddsDel,
                 contrast = c('condition','GMM','Azole'), res=res, type = 'ashr')
res_tbl <- as_tibble(res, rownames = "ENSEMBL")
p<-EnhancedVolcano(res,
                   lab = rownames(res),
                   x = 'log2FoldChange',
                   title = 'Del GMM versus Azole',
                   y = 'pvalue')

ggsave("plots/Volcano.Del_GMM-vs-Azole.pdf",p,width=6,height=8)


best_genes <- res_tbl %>%
  arrange(padj)  %>%
  head(9)
write_tsv(res_tbl,"results/Delta.tsv")

p <- as_tibble(counts(ddsDel[best_genes$ENSEMBL, ], normalize = TRUE),
               rownames = 'ENSEMBL') %>% 
  pivot_longer(names_to = "sample", values_to = "counts", -ENSEMBL) %>%  
  left_join(as_tibble(sampleTableDel, rownames = "sample")) %>% 
  filter(condition %in% c('GMM', 'Azole')) %>%
  ggplot(aes(x = sample, y = counts, fill = genotype)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ ENSEMBL, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))
p
ggsave("plots/barPlot.Delta.top_9_GMM_v_Azole.pdf",p,width=16,height=12)
# 5. WT drug vs del drug
# 6. Wt drug vs comp drug
# 7 WT drug vs OE drug
exprnamesAzole <-  (samplesFile %>% dplyr::filter(Condition == "Azole"))$SAMPLE
sampleTableAzole = subset(sampleTable,condition == "Azole")
rownames(sampleTableAzole) = exprnamesAzole

# check this is right order!
#rownames(sampleTable) = colnames(countdata)

countdataAzole = countdata[,exprnamesAzole]

ddsAzole <- DESeqDataSetFromMatrix(countData = countdataAzole,
                                 colData   = sampleTableAzole, 
                                 design    = ~ genotype)
nrow(ddsAzole)
ddsAzole <- ddsAzole[ rowSums(counts(ddsAzole)) > 1, ]
nrow(ddsAzole)
ddsAzole <- estimateSizeFactors(ddsAzole)
ddsAzole <- estimateDispersions(ddsAzole)
sizeFactors(ddsAzole)

ddsAzole <- DESeq(ddsAzole)
res <- results(ddsAzole)

res <- results(ddsAzole,
               contrast = c('genotype', 'WT','Del'))
res <- lfcShrink(ddsAzole,
                 contrast = c('genotype','WT','Del'), res=res, type = 'ashr')
res_tbl <- as_tibble(res, rownames = "ENSEMBL")
#selectLab = c('AFUB_000050-T', 'AFUB_000060-T'),
p<-EnhancedVolcano(res,
                   lab = rownames(res),
                   x = 'log2FoldChange',
                   FCcutoff = 2, 
                   pCutoff = 1e-10,
                   labSize = 3.0,
                   selectLab = c('AFUB_000050-T', 'AFUB_000060-T'),
                   title = 'Azole WT versus Delta',
                   y = 'pvalue',
                   legendPosition = 'right',
                   labCol = 'goldenrod',
                   labFace = 'bold',
                   drawConnectors = TRUE,
                   widthConnectors = 0.5,
                   colAlpha = 4/5,
                  )
p

ggsave("plots/Volcano.Azole.WT-vs-Delta.callOut.pdf",p,width=16,height=12)

p<-EnhancedVolcano(res,
                   lab = rownames(res),
                   x = 'log2FoldChange',
                   FCcutoff = 2, 
                   pCutoff = 1e-10,
                   labSize = 3.0,
                   selectLab = c('AFUB_000050-T', 'AFUB_000060-T'),
                   title = 'Azole WT versus Delta',
                   y = 'pvalue',
                   legendPosition = 'right')
p

ggsave("plots/Volcano.Azole.WT-vs-Delta.pdf",p,width=16,height=12)

res_tbl <- as_tibble(res, rownames = "ENSEMBL")
best_genes <- res_tbl %>%
  arrange(padj)  %>%
  head(9)

write_tsv(res_tbl,"results/Azole.WT-vs-Delta.tsv")


p <- as_tibble(counts(ddsAzole[best_genes$ENSEMBL, ], normalize = TRUE),
               rownames = 'ENSEMBL') %>% 
  pivot_longer(names_to = "sample", values_to = "counts", -ENSEMBL) %>%  
  left_join(as_tibble(sampleTableAzole, rownames = "sample")) %>% 
  ggplot(aes(x = sample, y = counts, fill = genotype)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ ENSEMBL, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))
p
ggsave("plots/barPlot.Azole.top_9.pdf",p,width=16,height=12)

p <- as_tibble(counts(ddsAzole[best_genes$ENSEMBL, ], normalize = TRUE),
               rownames = 'ENSEMBL') %>% 
  pivot_longer(names_to = "sample", values_to = "counts", -ENSEMBL) %>%  
  left_join(as_tibble(sampleTableAzole, rownames = "sample")) %>% 
  filter(condition %in% c('WT', 'Del')) %>%
  ggplot(aes(x = sample, y = counts, fill = genotype)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ ENSEMBL, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))
p
ggsave("plots/barPlot.Azole.top_9_WT_v_Delta.pdf",p,width=16,height=12)

res <- results(ddsAzole,
               contrast = c('genotype', 'WT','Comp'))
res <- lfcShrink(ddsAzole,
                 contrast = c('genotype','WT','Comp'), res=res, type = 'ashr')
res_tbl <- as_tibble(res, rownames = "ENSEMBL")
p<-EnhancedVolcano(res,
                   lab = rownames(res),
                   x = 'log2FoldChange',
                   title = 'Azole WT versus Comp',
                   selectLab = c('AFUB_000050-T1', 'AFUB_000060-T1'),
                   y = 'pvalue')

ggsave("plots/Volcano.Azole.WT-vs-Comp.pdf",p,width=6,height=8)

p <- as_tibble(counts(ddsAzole[best_genes$ENSEMBL, ], normalize = TRUE),
               rownames = 'ENSEMBL') %>% 
  pivot_longer(names_to = "sample", values_to = "counts", -ENSEMBL) %>%  
  left_join(as_tibble(sampleTableAzole, rownames = "sample")) %>% 
  filter(condition %in% c('WT', 'Comp')) %>%
  ggplot(aes(x = sample, y = counts, fill = genotype)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ ENSEMBL, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))
p
ggsave("plots/barPlot.Azole.top_9_WT_v_Comp.pdf",p,width=16,height=12)

res <- results(ddsAzole,
               contrast = c('genotype', 'WT','OE'))
res <- lfcShrink(ddsAzole,
                 contrast = c('genotype','WT','OE'), res=res, type = 'ashr')
res_tbl <- as_tibble(res, rownames = "ENSEMBL")
p <- EnhancedVolcano(res,
                     lab = rownames(res),
                     x = 'log2FoldChange',
                     title = 'Azole WT versus OverExpression',
                     selectLab = c('AFUB_000050-T1', 'AFUB_000060-T1'),
                     y = 'pvalue')

ggsave("plots/Volcano.Azole.WT-vs-OE.pdf",p,width=6,height=8)

p <- as_tibble(counts(ddsAzole[best_genes$ENSEMBL, ], normalize = TRUE),
               rownames = 'ENSEMBL') %>% 
  pivot_longer(names_to = "sample", values_to = "counts", -ENSEMBL) %>%  
  left_join(as_tibble(sampleTableAzole, rownames = "sample")) %>% 
  filter(condition %in% c('WT', 'OE')) %>%
  ggplot(aes(x = sample, y = counts, fill = genotype)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ ENSEMBL, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))
p
ggsave("plots/barPlot.Azole.top_9_WT_v_OE.pdf",p,width=16,height=12)
