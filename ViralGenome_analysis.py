#!/usr/bin/env python
# ViralGenome_Analysis.py

import os, cmath, math, sys, glob, subprocess, re, argparse, shutil, datetime, csv
import numpy as np
from matplotlib_venn import venn2
import pandas as pd
from collections import defaultdict
from operator import itemgetter
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import defaultdict
from optparse import OptionParser
csv.register_dialect("textdialect",delimiter='\t')

global sampleName
global outfilepath
global logFile
global logOpen

### Parsing arguments ###

parser = argparse.ArgumentParser(description="ViralGenome_Analysis.py", epilog="Example: ViralGenome_Analysis.py -i rawdata/example_MMhur_R1.fastq rawdata/example_MMhur_R2.fastq -v jfh1 -n MMhur_JFH1 -o results")
parser.add_argument('-i', metavar=("INPUT1", "INPUT2"), nargs=2, help="2 input fastq or fastq.gz files separated by a space", required=True)
parser.add_argument('--trimmed', action='store_true', help="flag if files are already trimmed")
parser.add_argument('-v', metavar='VIRUS', help="Name of viral genome to be mapped against", required=True)
parser.add_argument('-n', metavar='NAME', help="Name of output directory", required=True)
parser.add_argument('-o', metavar='OUTPUT', help="Name of directory where output directory will be made", required=True)
parser.add_argument('-f', metavar='N', type=int, help="First base to keep on 5' end of each read. Default is 14.", default=14)
parser.add_argument('-a', metavar='ADAPTER', help="3' adapter to trim from the end of each read. Default is AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG.", default='AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG')
parser.add_argument('-t', metavar='THRESHOLD', type=int, help="Stringency of RT stop filtering, where higher numbers indicate greater stringency. Default is 3.", default=3)
parser.add_argument('-m', metavar='MAPQ', type=int, help="Minimum MAPQ score allowed. Default is 42.", default=42)
parser.add_argument('-q', metavar='Q', type=int, help="Minimum quality score to keep during filtering. Default is 25.", default=25)
parser.add_argument('-p', metavar='P', type=int, help="Percentage of bases that must have quality > q during filtering. Default is 80.", default=80)
parser.add_argument('-l', metavar='L', type=int, help="Minimum length of read. Default is 15.", default=15)

args = parser.parse_args()

# input files
reads = args.i
for fn in reads:
	if not glob.glob(fn):
		print "Error: input file " + fn + " not accessible. Exiting."
		exit()
		
# sample name and output directory
sampleName = args.n
outfilepath = args.o
if not glob.glob(outfilepath):
	print "Error: output directory " + outfilepath + " not accessible. Exiting."
	exit()
outfilepath = outfilepath + '/%s/'%sampleName
if not glob.glob(outfilepath): os.system("mkdir " + outfilepath)

# virus files
virus = args.v
viral_index = os.getcwd() + '/docs/%s/%s' % (virus, virus)
viral_genome = os.getcwd()+'/docs/%s/%s.fa' % (virus, virus)
if not glob.glob(viral_genome):
	print "Error: viral genome file %s not accessible. Exiting." % viral_genome
	exit()

if not glob.glob(viral_index + "*.bt2"):
	print "Error: viral genome index files starting with %s not accessible. Exiting." % viral_index
	exit()

# Create log and start pipeline
logFile = outfilepath + "runLog"
logOpen = open(logFile, 'w')

### Parameters ###
iCLIP3pBarcode = args.a # Barcode sequence to trim from reads.
mapq = args.m # Minimum MAPQ score allowed
q = args.q # Minimum quality score to keep during filtering.
p = args.p # Percentage of bases that must have quality > q during filtering.
l = args.l +  args.f - 1 # Minimum length of read + 5' adapter
iCLIP5pBasesToTrim = args.f # Number of reads to trim from 5' end of clip reads + 1 (this number represents the first base kept)
k = '1' # k = N distinct, valid alignments for each read in bt2 mapping.
threshold = args.t # Sum of RT stops (for both replicates) required to keep it. 
expand = 15 # Bases to expand around RT position after RT stops are merged.
CLIPPERoutNameDelim = '_' # Delimiter that for splitting gene name in the CLIPper windows file.

### start running pipeline ###
now=datetime.datetime.now()
logOpen.write("Timestamp:%s\n"%str(now))
logOpen.write("\n###Parameters used###\n")
logOpen.write("3' barcode:%s\n'"%iCLIP3pBarcode)
logOpen.write("Minimum quality score (q):%s\n"%q)
logOpen.write("Percentage of bases with > q:%s\n"%p)
logOpen.write("5' bases to trim:%s\n'"%iCLIP5pBasesToTrim)
logOpen.write("k distinct, valid alignments for each read in bt2 mapping:%s\n"%k)
logOpen.write("Threshold for minimum number of RT stops:%s\n"%threshold)
logOpen.write("Bases for expansion around conserved RT stops:%s\n"%expand)
logOpen.write("\n\n\n")

def remove_dup(reads, q, p):
	uniq_reads = []
	for inread in reads:
		outread = outfilepath + os.path.basename(inread)
		is_compressed = True if inread[-3:] == '.gz' else False
		if is_compressed:
			outread=outread.replace(".fastq.gz", "_nodup.fasta")		
		else:
			outread=outread.replace(".fastq", "_nodup.fasta")
			
		uniq_reads.append(outread)
		if glob.glob(outread):
			uniq_reads=uniq_reads+[outread]
			print "Filtering and duplicate removal already done."
			logOpen.write("Filtering and duplicate removal already done.\n")
			continue
			
		if is_compressed:
			cmd_1 = "gunzip -c {}".format(inread)
		else:
			cmd_1 = "cat {}".format(inread)
			
		cmd_2 = "fastq_quality_filter -Q33 -q{} -p{}".format(q, p)
		cmd_3 = "fastx_collapser -Q33 > {}".format(outread)
		full_cmd = ' | '.join([cmd_1, cmd_2, cmd_3])
		print full_cmd
		os.system(full_cmd)
	return uniq_reads

if not args.trimmed: 
	print "Removing duplicates"
	dup_removed_reads = remove_dup(reads, q, p)

def trim(reads, adapter3p, l, n):
	# Usage: Trims a specified number of bases from the 5' end of each read.
	# Input: List of fastq files.
	# Output: List of 5p trimmed files.
	program = os.getcwd() + "/bin/fasta_to_fastq.pl"
	trimmedReads = []
	for inread in reads:
		outread = inread.replace("_nodup.fasta", "_trimmed.fastq")
		trimmedReads.append(outread)
		if glob.glob(outread):
			print "5' and 3' barcode trimming already done."
			logOpen.write("5' and 3' barcode trimming already done.\n")
			continue
		
		cmd_1 = "perl {} {}".format(program, inread)
		cmd_2 = "fastx_clipper -n -l{} -Q33 -a {}".format(l, adapter3p)
		cmd_3 = "fastx_trimmer -f{} -Q33 > {}".format(n, outread)
		full_cmd = ' | '.join([cmd_1, cmd_2, cmd_3])
		print full_cmd
		os.system(full_cmd)
		logOpen.write("Perform 5' and 3' barcode trimming.\n")
	return trimmedReads

if not args.trimmed: 
	print "Trimming 5' and 3'"
	processed_reads = trim(dup_removed_reads, iCLIP3pBarcode, l, iCLIP5pBasesToTrim)
else: processed_reads = reads

def runBowtie(processed_reads, index):
	# Usage: Read mapping.
	# Input: Fastq files of replicate trimmed read files.
	# Output: Path to samfile for each read.

	viral_sam = []
	for infastq in processed_reads:
		mapped = infastq.replace(".fastq", "_mappedTo_%s.sam" % virus)
		mapped = outfilepath + os.path.basename(mapped)
		viral_sam.append(mapped)
			
		cmd = "bowtie2 -p 8 -x {} {} -S {} > {} 2>&1".format(index, infastq, mapped, mapped + '_stats.txt')
		print cmd
		os.system(cmd)
		
		logOpen.write("Perform mapping.")
	return viral_sam

print "\nRun mapping to viral index."  
viral_sam = runBowtie(processed_reads, viral_index)

def run_samtools(samfiles, mapq):
	# Usage: Samfile processing (also takes unique mappers only)
	# Input: Sam files from Bowtie mapping.
	# Output: Sorted bedFiles.
	program = 'samtools'
	program2 = 'bamToBed'
	out_bedfiles = []
	for samfile in samfiles:
		bamfile_sort = samfile.replace('.sam', '_sorted') 
		bedfile = bamfile_sort.replace('_sorted', '.bed') 
		
		out_bedfiles.append(bedfile)
		if glob.glob(bedfile):
			print "Samtools already done."
			logOpen.write("Samtools already done.\n")
			continue
			
		cmd_1 = "cat {} | samtools view -q {} -Suo - - | samtools sort - {}".format(samfile, mapq, bamfile_sort)
		cmd_2 = "bamToBed -i {} > {}".format(bamfile_sort + '.bam', bedfile)
		# cmd_3 = "samtools flagstat {} > {}".format(bamfile_sort + '.bam', bamfile_sort + '_stats.txt')
		# cmd_4 = "rm -f {}".format(samfile)
		print cmd_1
		os.system(cmd_1)
		print cmd_2
		os.system(cmd_2)
		
	return out_bedfiles

print "\nRun samtools."
logOpen.write("Run samtools.\n")
viral_bed = run_samtools(viral_sam, mapq)

# Index
viral_genome_build = np.genfromtxt(viral_genome, dtype='string')
viral_genome_size = len(viral_genome_build[1])

# Data frames to store HCV data
storage_jfh1 = pd.DataFrame()
indexPositions = np.arange(0, viral_genome_size, 1)
storeGenomeHits = pd.DataFrame(index=indexPositions)

coverageHistograms=[]
for bedFile in viral_bed:
	posToCount_plus = defaultdict(lambda: 0)
	posToCount_minus = defaultdict(lambda: 0)
	handle = os.path.basename(bedFile)
	with open(bedFile, 'r') as ifile:
		reader = csv.reader(ifile, 'textdialect')
		for row in reader:
			if row[-1] == '+': posToCount_plus[int(row[1])] += 1 # Plus strand --> Start of read
			else: posToCount_minus[int(row[2]) - 2] += 1 # Minus strand --> End of read - 1 (since 1-based coord)
		
		genomeHits_plus = []
		genomeHits_minus = []
		for i in range(viral_genome_size):
			genomeHits_plus.append(posToCount_plus[i])
			genomeHits_minus.append(posToCount_minus[i])
		storeGenomeHits['Plus strand ' + handle] = genomeHits_plus
		storeGenomeHits['Minus strand ' + handle] = genomeHits_minus

bases = [viral_genome_build[1][pos] for pos in range(viral_genome_size)]
storeGenomeHits['Base']=bases
storeGenomeHits.to_csv(outfilepath + '/numReads.txt')







# Removing things
os.chdir(outfilepath)
os.system("rm -f *.sam")

logOpen.close()



