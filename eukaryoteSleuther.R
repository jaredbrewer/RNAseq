# Sleuth Analysis of Kallisto Data
# Created: 29 September 2017
# Last Edited: 16 August 2022
# Jared Brewer

# You have to provide exactly two things:
  # The directory where you did everything and set the working directory to that location using setwd()
  # A sample table formatted with two columns: sample and condition and named sample_table.txt
# After running, you'll get an interactive window where you can play with your data and look for effects.
# There is merit in breaking out comparisons if you have a bunch of conditions...

setwd("~/Desktop/granuloma_rnaseq/")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install(c("rhdf5", "biomaRt", "ReactomePA"))

install.packages("devtools")
devtools::install_github("pachterlab/sleuth")

library(sleuth)
library(tidyverse)
library(biomaRt)

# drerio can be exchanged for any other [g]enus [species] and it should work. hsapiens, mmusculus, etc. 
# Non-animals may have to use the host "ensemblgenomes.org" - keep that in mind.
# This does not work at all for bacteria thanks to Ensembl discontinuing bacterial marts.

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'https://ensembl.org')
t2g <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- rename(t2g, target_id = ensembl_transcript_id_version, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

files.kal <- list.files(pattern = "./*_quant")
files.dirs <- sapply(files.kal, function(id) file.path(".", id))
sample.table <- read.table("./sample_table.txt", header = T)
sample.table <- dplyr::mutate(sample.table, path = files.dirs)

so <- sleuth_prep(sample.table, target_mapping = t2g, extra_bootstrap_summary = T, read_bootstrap_tpm = T)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
# so <- sleuth_wt(so, which_beta = 'condition', which_model = 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, pval <= 0.05)

sleuth_live(so)

de <- subset(sleuth_table, pval < 0.05)
paths <- enrichPathway(gene = trimws(de$entrez), pvalueCutoff = 0.2, readable = T, organism = "human")





