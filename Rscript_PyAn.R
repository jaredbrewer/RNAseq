# Modification of Salmon_09202017.R for Integration into Python

library(DESeq2)
library(tximport)
library(readr)
library(tidyverse)
library(ensembldb)
library(org.Dr.eg.db)
library(AnnotationHub)
library(biomaRt)
library(GenomicFeatures)

# User inputs - need to figure out how to parse # 
setwd(fastq_dir + '/' + 'quant')
sampleTable <- read.table(fastq_dir + '/' + 'jared_samples.txt', header = T)



dir <- list.files()
samples <- as.data.frame(dir, header = T)
files <- file.path(dir, "quant.sf")
jared.files <- paste(".", dir, "quant.sf", sep="/")

rownames(sampleTable) <- sampleTable$ids
sampleTable[,c("ids","treat")]
names(files) <- sampleTable$ids

# This is all zebrafish specific code that I need to integrate with my previous Genus species input

zf_mart <- useMart("ensembl", "drerio_gene_ensembl")
zf_anno <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id_version", "chromosome_name", "start_position", "end_position", "strand", "description"), mart = zf_mart)
zf_anno$Chr <- paste("chr", zf_anno$chromosome_name, sep = "")
zf_anno$Start <- zf_anno$start_position
zf_anno$End <- zf_anno$end_position
zf_anno$strand <- gsub('-1', '-', zf_anno$strand)
zf_anno$strand <- gsub('1', '+', zf_anno$strand)
zf_anno$Strand <- zf_anno$strand
zf_ref <- zf_anno[, c(1,2,8,9,10,11,7)]
tx2gene <- zf_ref[, c(2,1)]
