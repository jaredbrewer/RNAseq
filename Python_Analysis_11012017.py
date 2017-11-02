#!/usr/bin/env python3

import system, re, os, Bio
from biomart import BiomartServer

fastq_dir = input("Enter the directory of your FASTQ files:")
subprocess.call(['cd', fastq_dir])

if not os.path.exists(fastq_dir + 'index/'):
        os.makedirs(fastq_dir + 'index/')

organism_name = input("Input the 'Genus species' for your reference organism:")

# This install HomeBrew, which is then used to install Salmon in the next step.
subprocess.call(['./brew_installer.py'])

# Installs the latest version of Salmon from HomeBrew
subprocess.call(['brew', 'install', 'salmon'])

# Installs the latest version of BioMart, for transcriptome retrieval. 
subprocess.call(['pip', 'install', 'biomart'])

# This does give a warning if Salmon is already installed and up to date, but that is not an issue.

salmon = salmon
index = index
server = BiomartServer("http://ensembl.org/biomart")
ensembl = server.datasets['ENSEMBL_MART_ENSEMBL']

# May need to try MySQL retrieval of the transcriptomes if BioMart keeps causing a fuss.