# Automatically Run Downstream Analysis from Salmon Alignment

# Dependencies: 
library(DESeq2)
library(tximport)
library(readr)
library(tidyverse)
library(ensembldb)
library(org.Dr.eg.db)
library(AnnotationHub)
library(biomaRt)
library(GenomicFeatures)
library(GOexpress)
library(ReactomePA)

# setwd() required for use, can I pass this from Python? Surely.

doEverything <- function(sampleTable, ref, treat, volcano = T, PCA = T, GO = T, React = T) 
{
  salmonAnalyzer()
  salmonCompare()
  salmonGraph()
  if (GO == T) {
    salmonGO()
  }
  if (PCA == T) {
    salmonPCA()
  }
  if (React == T) {
    salmonReact()
  }
} 

# Make this an option from the Python command line? I think that would be nifty. 

salmonAnalyzer <- function(dir, sampleTable, organism) {
  dir <- list.files()
  samples <- as.data.frame(dir, header = T)
  files <- file.path(dir, "quant.sf")
  quant.files <- paste(".", dir, "quant.sf", sep = "/")
  
  sampleTable <- read.table(sampleTable, header = T)
  rownames(sampleTable) <- sampleTable$ids
  sampleTable[,c("ids","treat")]
  names(files) <- sampleTable$ids
  
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
  
  import <- tximport(quant.files, type = "salmon", tx2gene = tx2gene)
  dds <- DESeqDataSetFromTximport(import, sampleTable, ~treat)
  dds$treat <- relevel(dds$treat, ref = "Mock")
  return(DESeq(dds))
}

# Pass organism from Python, pass dir, so the user just needs to provide a properly formatted sampleTable

salmonCompare <- function(dds = analysis, ref, treat, write = T) {
  results <- results(dds, contrast = c("treat", treat, ref))
  results <- results[order(results$pvalue),]
  results$symbol <- mapIds(dr, keys = row.names(results), column="SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  results$entrez <- mapIds(dr, keys = row.names(results), column="ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  results$name <- mapIds(dr, keys = row.names(results), column="GENENAME", keytype = "ENSEMBL", multiVals = "first")
  if (write == T) {write.csv(results, "comparison.csv")}
  return(results)
} 

# this will be simple, but include write.csv
salmonGraph <- function(dds = analysis, comparison = results, volcano = T, PCA = T) {
  if (volcano == T) {
    rld <- rlog(dds, blind = F)
    pca <- plotPCA(rld, intgroup = c("treat"))
    pca # This isn't quite working yet, but I'm not sure why.
  }
  if (PCA == T) {
    volcano <- ggplot(as.data.frame(comparison), aes(log2FoldChange, -log10(pvalue), color = pvalue < 0.05)) + 
      geom_point() + scale_color_manual(values = setNames(c('red','black'),c(T, F))) + 
      geom_text_repel(data = head(as.data.frame(comparison), 20), aes(label = symbol)) + guides(color = F)
    volcano
  }
} 
# Simple true/false. If the user wants better, they can do it themselves.
salmonGO <- function() {
  
} # Incorporate GOExpress here
salmonKEGG <- function () {
  
} # Maybe unnecessary, but I'll see if I can code it.
salmonReact <- function () {
  
} # Use ReactomePA to do better analysis for pathways and make pretty graphs.

ah <- AnnotationHub()
query(ah, "EnsDb")
ahDb <- query(ah, pattern = c("Danio Rerio", "EnsDb", 88))
ahEdb <- ahDb[[1]]
