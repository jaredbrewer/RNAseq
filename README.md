# RNAseq Analysis Pipeline with Kallisto/Sleuth

This pipeline is meant to simplify the process of analyzing (bulk) RNA sequencing data using a user-friendly interface. The Python script (__eukaryoteKallisto.py__) guides the user through the process of providing all the needed information while the R script (__eukaryoteSleuther.R__) requires minimal changes in order to function.

The user needs to provide:
* FASTQ files
* Sample table in the following format:

	| sample | condition |
	|--------|-----------|
	| [filename] | [simple name] |

The script will download the needed reference cDNA, make the index file, and analyze the provided FASTQ files. Then, the directory has to be provided in R to find the files and the script can be run for visual analysis in Sleuth.

This repository includes an macOS UNIX executable binary of Kallisto based on version 0.48 that includes the required support for HDF5 in order to analyze the abundance data in Sleuth. Users on Windows are encouraged to download versions of Kallisto prior to version 0.46.2 from the [Pachter Lab](https://github.com/pachterlab/kallisto/releases). Linux users may also use their earlier binaries or compile their own by changing the first option in the CMakeLists.txt file from ```option(USE_HDF5 "Compile with HDF5 support" OFF)``` to ```option(USE_HDF5 "Compile with HDF5 support" ON)```. (License included for the purposes of redistribution of these binaries).
