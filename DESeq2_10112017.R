# DESeq2 Analysis of HISAT2 Data
# Created: 21 September 2017
# Last Edited: 21 September 2017
# Jared Brewer

setwd("~/RNAseq/HISAT2")

library(DESeq2)
library(Rsubread)

sorted.dir <- list.files(path = "./sorted", pattern = ".sorted.bam")
dedup.dir <- list.files(path = "./dedup", pattern = ".dedup.bam")

sorted.path <- paste("./sorted", sorted.dir, sep="/")
dedup.path <- paste("./dedup", dedup.dir, sep="/")

sorted.counts <- featureCounts(files = sorted.path, annot.ext = anno.g, isGTFAnnotationFile = F,
                              nthreads = 4, verbose = T)
dedup.counts <- featureCounts(files = dedup.path, annot.ext = anno.g, isGTFAnnotationFile = F,
                              nthreads = 4, verbose = T)

sampleTable <- read.table("~/RNAseq/HISAT2/jared_samples.txt", header = T)
rownames(sampleTable) <- sampleTable$ids
sampleTable[,c("ids","treat")]
rownames(sampleTable) <- colnames(sampleTable$counts)
names(sorted.dir) <- sampleTable$ids
names(dedup.dir) <- sampleTable$ids

sorted.countData <- sorted.counts$counts
dedup.countData <- dedup.counts$counts

sorted.dds <- DESeqDataSetFromMatrix(countData = sorted.countData, colData = sampleTable, design = ~ treat)
sorted.dds <- DESeq(sorted.dds)

dedup.dds <- DESeqDataSetFromMatrix(countData = dedup.countData, colData = sampleTable, design = ~ treat)
dedup.dds <- DESeq(dedup.dds)

sorted.MP <- results(sorted.dds, contrast = c("treat", "pcaA", "Mock"))
sorted.MP <- sorted.MP[order(sorted.MP$pvalue),]
sorted.MP$symbol <- mapIds(dr, keys = row.names(sorted.MP), column="SYMBOL", keytype = "ENSEMBL", multiVals = "first")
sorted.MP$entrez <- mapIds(dr, keys = row.names(sorted.MP), column="ENTREZID", keytype = "ENSEMBL", multiVals = "first")
sorted.MP$name <- mapIds(dr, keys = row.names(sorted.MP), column="GENENAME", keytype = "ENSEMBL", multiVals = "first")

sorted.MW <- results(sorted.dds, contrast = c("treat", "WT", "Mock"))
sorted.MW <- sorted.MW[order(sorted.MW$pvalue),]
sorted.MW$symbol <- mapIds(dr, keys = row.names(sorted.MW), column="SYMBOL", keytype = "ENSEMBL", multiVals = "first")
sorted.MW$entrez <- mapIds(dr, keys = row.names(sorted.MW), column="ENTREZID", keytype = "ENSEMBL", multiVals = "first")
sorted.MW$name <- mapIds(dr, keys = row.names(sorted.MW), column="GENENAME", keytype = "ENSEMBL", multiVals = "first")

sorted.PW <- results(sorted.dds, contrast = c("treat", "WT", "pcaA"))
sorted.PW <- sorted.PW[order(sorted.PW$pvalue),]
sorted.PW$symbol <- mapIds(dr, keys = row.names(sorted.PW), column="SYMBOL", keytype = "ENSEMBL", multiVals = "first")
sorted.PW$entrez <- mapIds(dr, keys = row.names(sorted.PW), column="ENTREZID", keytype = "ENSEMBL", multiVals = "first")
sorted.PW$name <- mapIds(dr, keys = row.names(sorted.PW), column="GENENAME", keytype = "ENSEMBL", multiVals = "first")

dedup.MP <- results(dedup.dds, contrast = c("treat", "pcaA", "Mock"))
dedup.MP <- dedup.MP[order(dedup.MP$pvalue),]
dedup.MP$symbol <- mapIds(dr, keys = row.names(dedup.MP), column="SYMBOL", keytype = "ENSEMBL", multiVals = "first")
dedup.MP$entrez <- mapIds(dr, keys = row.names(dedup.MP), column="ENTREZID", keytype = "ENSEMBL", multiVals = "first")
dedup.MP$name <- mapIds(dr, keys = row.names(dedup.MP), column="GENENAME", keytype = "ENSEMBL", multiVals = "first")

dedup.MW <- results(dedup.dds, contrast = c("treat", "WT", "Mock"))
dedup.MW <- dedup.MW[order(dedup.MW$pvalue),]
dedup.MW$symbol <- mapIds(dr, keys = row.names(dedup.MW), column="SYMBOL", keytype = "ENSEMBL", multiVals = "first")
dedup.MW$entrez <- mapIds(dr, keys = row.names(dedup.MW), column="ENTREZID", keytype = "ENSEMBL", multiVals = "first")
dedup.MW$name <- mapIds(dr, keys = row.names(dedup.MW), column="GENENAME", keytype = "ENSEMBL", multiVals = "first")

dedup.PW <- results(dedup.dds, contrast = c("treat", "WT", "pcaA"))
dedup.PW <- dedup.PW[order(dedup.PW$pvalue),]
dedup.PW$symbol <- mapIds(dr, keys = row.names(dedup.PW), column="SYMBOL", keytype = "ENSEMBL", multiVals = "first")
dedup.PW$entrez <- mapIds(dr, keys = row.names(dedup.PW), column="ENTREZID", keytype = "ENSEMBL", multiVals = "first")
dedup.PW$name <- mapIds(dr, keys = row.names(dedup.PW), column="GENENAME", keytype = "ENSEMBL", multiVals = "first")

sorted.rld <- rlog(sorted.dds, blind = F)
plotPCA(sorted.rld, intgroup = c("treat"))

dedup.rld <- rlog(dedup.dds, blind = F)
plotPCA(dedup.rld, intgroup = c("treat"))

library(pathview)
library(gage)

kegg.set <- kegg.gsets("dre")

sMW.fc <- sorted.MW$log2FoldChange
names(sMW.fc) <- sorted.MW$entrez
sMP.fc <- sorted.MP$log2FoldChange
names(sMP.fc) <- rownames(sorted.MP)
sPW.fc <- sorted.PW$log2FoldChange
names(sPW.fc) <- rownames(sorted.PW)
dMW.fc <- dedup.MW$log2FoldChange
names(dMW.fc) <- rownames(dedup.MW)
dMP.fc <- dedup.MP$log2FoldChange
names(dMP.fc) <- rownames(dedup.MP)
dPW.fc <- dedup.PW$log2FoldChange
names(sPW.fc) <- rownames(dedup.PW)

sMW.kegg <- gage(sMW.fc, gsets = kegg.set)
sMP.kegg <- gage(sMP.fc, gsets = kegg.set)
sPW.kegg <- gage(sPW.fc, gsets = kegg.set)
dMW.kegg <- gage(dMW.fc, gsets = kegg.set)
dMP.kegg <- gage(dMP.fc, gsets = kegg.set)
dPW.kegg <- gage(dPW.fc, gsets = kegg.set)



