# RNAseq Analysis Pipeline with Kallisto/Sleuth

This pipeline is meant to simplify the process of analyzing (bulk) RNA sequencing data using a user-friendly interface. The Python script (__eukaryoteKallisto.py__) guides the user through the process of providing all the needed information while the R script (__eukaryoteSleuther.R__) requires minimal changes in order to function.

The user needs to provide:
* FASTQ files
* Sample table in the following format:

	| sample | condition |
	|--------|-----------|
	| [filename] | [simple name] |

The script will download the needed reference cDNA, make the index file, and analyze the provided FASTQ files. Then, the directory has to be provided in R to find the files and the script can be run for visual analysis in Sleuth.
