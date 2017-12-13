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

quants <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(quants)

dds <- DESeqDataSetFromTximport(quants, sampleTable, ~treat)
dds$treat <- relevel(dds$treat, ref = "Mock")
dds <- DESeq(dds)

# results(dds, contrast = c("treat", "Mock", "pcaA"))
# Mock.pcaA[order(Mock.pcaA$pvalue),]
# Mock.WT$symbol <- mapIds(dr, keys = row.names(Mock.WT), column="SYMBOL", keytype = "ENSEMBL", multiVals = "first")
# Mock.WT$entrez <- mapIds(dr, keys = row.names(Mock.WT), column="ENTREZID", keytype = "ENSEMBL", multiVals = "first")
# Mock.WT$name <- mapIds(dr, keys = row.names(Mock.WT), column="GENENAME", keytype = "ENSEMBL", multiVals = "first")

jared.rld <- rlog(dds, blind = F)
plotPCA(jared.rld, intgroup = c("treat"))
salmon.volcano <- ggplot(pw.res, aes(log2FoldChange, -log10(pvalue), color = pvalue < 0.05)) + 
  geom_point() + scale_color_manual(values = setNames(c('red','black'),c(T, F))) + 
  geom_text_repel(data = head(pw.res, 20), aes(label = symbol)) + guides(color = F)

# I would bet that I can rewrite this entire section in Python and save myself the burden of crossing languages over and over.

# This is all zebrafish specific code that I need to integrate with my previous Genus species input

# zf_mart <- useMart("ensembl", "drerio_gene_ensembl")
# zf_anno <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id_version", "chromosome_name", "start_position", "end_position", "strand", "description"), mart = zf_mart)
# zf_anno$Chr <- paste("chr", zf_anno$chromosome_name, sep = "")
# zf_anno$Start <- zf_anno$start_position
# zf_anno$End <- zf_anno$end_position
# zf_anno$strand <- gsub('-1', '-', zf_anno$strand)
# zf_anno$strand <- gsub('1', '+', zf_anno$strand)
# zf_anno$Strand <- zf_anno$strand
# zf_ref <- zf_anno[, c(1,2,8,9,10,11,7)]
# tx2gene <- zf_ref[, c(2,1)]

# library(ensembldb)
# library(org.Dr.eg.db)
# library(AnnotationHub)
# ah <- AnnotationHub()
# query(ah, "EnsDb")
# ahDb <- query(ah, pattern = c("Danio Rerio", "EnsDb", 88))
# ahEdb <- ahDb[[1]]
# 
# ### Script section copied because new Ensembl versions have added .X to their transcripts - very annoying
# mart <- biomaRt::useMart(biomart = "ensembl", dataset = "drerio_gene_ensembl")
# t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name", "description", "transcript_biotype"), mart = mart)
# t2g$target_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep=".") # append version number to the transcript ID
# t2g[,c("ensembl_transcript_id","transcript_version")] <- list(NULL) # delete the ensembl transcript ID and transcript version columns
# t2g <- dplyr::rename( t2g, gene_symbol = external_gene_name, full_name = description, biotype = transcript_biotype )
# tx2gene <- t2g[, c(5,1)]
# write.csv(tx2gene, "./tx2gene.csv")
# 
# 
# dr.anno <- anno[, 8:13]
# dr.anno <- dr.anno[, c(2,3,5,4,6,1)]
# 
# anno.g <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id_version", "chromosome_name", "start_position", "end_position", "strand", "description"), mart = zf_mart)
# anno.g$GeneID <- anno.g$ensembl_gene_id
# anno.g$Chr <- paste("chr", anno.g$chromosome_name, sep = "")
# anno.g$Start <- anno.g$start_position
# anno.g$End <- anno.g$end_position
# anno.g$strand <- gsub('-1', '-', anno.g$strand)
# anno.g$strand <- gsub('1', '+', anno.g$strand)
# anno.g$Strand <- anno.g$strand
# g.anno <- anno.g[, 8:12]











