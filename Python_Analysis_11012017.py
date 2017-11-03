#!/usr/bin/env python3

import re, os, ftplib, subprocess, glob, sys, shutil

# This gives the script some self awareness. It finds itself and changes the working directory to that path (temporarily).
# This is important for executing the brew_installer.sh script.
script_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(script_path)
# print(os.getcwd())

# This will check to see if several important system programs are installed in hierarchical order
# If they are not, then it executes a script to install them, may require user password. 
if not shutil.which('xcode-select'):
	if not shutil.which('brew'):
		if not shutil.which('salmon'):
			subprocess.call(['./brew_installer.sh'])

from termcolor import colored, cprint
import Bio

fastq_dir = input("Enter the directory of your FASTQ files:")
fastq_dir = fastq_dir.strip()
try:
	os.chdir(fastq_dir)
except FileNotFoundError: 
	text = colored("Looks like that directory does not exist - restart the script and try dragging and dropping the folder directly into Terminal.", "red")
	print(text)
	sys.exit(1)
except PermissionError:
	print(text)
	sys.exit(1)

if not os.path.exists(fastq_dir + '/index/'):
        os.makedirs(fastq_dir + '/index/')
if not os.path.exists(fastq_dir + '/quant/'):
        os.makedirs(fastq_dir + '/quant/')

organism_name = input("Input the 'Genus species' for your reference organism:")
org_split = organism_name.lower().split()
org_dir = "_".join(org_split)

pattern = '*.cdna.all.fa.gz'
ref_cdna = 'ref_cdna.fa.gz'

# This is where the user-friendliness comes in: use Ensembl's highly regular patterning - 
# to fetch the needed file for supported organisms. 

try: 
	with ftplib.FTP('ftp.ensembl.org') as ftp:
		ftp.login('anonymous')
		ftp.cwd('/pub/release-90/fasta/{}/cdna/'.format(org_dir))
		for filename in ftp.nlst(pattern): 
			fhandle = open(filename, 'wb') 
			ftp.retrbinary('RETR ' + filename, fhandle.write)
			fhandle.close()
except ftplib.error_perm:
	text = colored("It appears that your organism is not supported. Try running again with a closely related species or check spelling.", 'red')
	print(text)
	sys.exit(1)

# Define some potentially useful variables that I can plug into subprocess.
salmon = 'salmon'
index = 'index'
quant = 'quant'

# Fetch the index file (name independent, it just needs to end in the pattern defined above)
# All Ensembl cDNA files should have that precise formatting.

ref_cdna = glob.glob(pattern)
ref_cdna = ref_cdna[0]

# Find all the FASTQ files in the directory.
fastqs = glob.glob('*fastq.gz')

subprocess.call([salmon, index, '-t', fastq_dir + '/' + ref_cdna, '-i', fastq_dir + '/index', '-p', '4'])
for fastq in fastqs:
	subprocess.call([salmon, quant, '-i', fastq_dir + '/index/', 
	'-r', fastq_dir + '/' + fastq, 
	'-p', '4', '-o', fastq_dir + '/' + 'quant' + '/' + fastq + '_quant', 
	'--seqBias', '--gcBias', '-l', 'SF'])
	
#### Need to add in analysis portion next, but that will be significantly more difficult. ####

