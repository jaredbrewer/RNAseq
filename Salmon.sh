#! /bin/sh

cd ~/RNAseq/

export PATH=$PATH:~/RNAseq/Salmon

mkdir ~/RNAseq/Salmon/index
cd ~/RNAseq/Salmon/index

sudo DYLD_FALLBACK_LIBRARY_PATH=~/RNAseq/Salmon ./salmon -h

# Index the reference transcriptome 

salmon index -t /Users/jared/RNAseq/Salmon/GCF_GRCz11_RNA.fna -i ~/RNAseq/Salmon/refseq -p 4

# Map reads to the reference transcriptome
for hy in $(ls *.fastq.gz | sed 's/_R1_001.fastq.gz//')
do
salmon quant -i ~/RNAseq/Salmon/refseq -l A \
	-r ~/RNAseq/reads/${hy}_R1_001.fastq.gz \
	-p 4 -o /Volumes/JB_HD2/RNAseq/Salmon/refseq/${hy}_quant \
	--seqBias --gcBias -l SF 
done
# 