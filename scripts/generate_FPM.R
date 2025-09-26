#!/usr/bin/env Rscript

library(DESeq2)
library(tximport)
library(dplyr)
library(tidyverse)

countdata <- read.table("results/read_count_A1163.tsv", header=TRUE, row.names=1)

countdata <- countdata[ ,6:ncol(countdata)]
countdata <- as.matrix(countdata)
head(countdata)

samplesFile <- read_csv("samples.csv",col_names=TRUE) %>% arrange(SAMPLE)
exprnames <- samplesFile$SAMPLE

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

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData   = sampleTable,
                              design    = ~ condition + genotype)
nrow(dds)
write_tsv(as_tibble(fpm(dds),rownames='TRANSCRIPT'), 'results/FPM.tsv')
