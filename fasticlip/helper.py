# helper.py

import os, cmath, math, sys, glob, subprocess, re, argparse, shutil, datetime, csv, commands
import numpy as np
import cfg
import pandas as pd
from collections import defaultdict
from operator import itemgetter
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from optparse import OptionParser
mpl.rcParams['savefig.dpi'] = 200
mpl.rcParams['path.simplify'] = True
csv.register_dialect("textdialect",delimiter='\t')
a=[]
b=[]

################
### TRIMMING ###
################

def log(str):
	cfg.logOpen.write(str + '\n')
	
def remove_dup(reads, q, p):
	uniq_reads = []
	for inread in reads:
		outread = cfg.outfilepath + os.path.basename(inread)
		is_compressed = True if inread[-3:] == '.gz' else False
		if is_compressed:
			outread=outread.replace(".fastq.gz", "_nodup.fasta")        
		else:
			outread=outread.replace(".fastq", "_nodup.fasta")
			
		uniq_reads.append(outread)
		if glob.glob(outread):
			log("Filtering and duplicate removal already done.")
			continue
			
		if is_compressed:
			cmd_1 = "gunzip -c {}".format(inread)
		else:
			cmd_1 = "cat {}".format(inread)
			
		cmd_2 = "fastq_quality_filter -Q33 -q{} -p{}".format(q, p)
		cmd_3 = "fastx_collapser -Q33 > {}".format(outread)
		full_cmd = ' | '.join([cmd_1, cmd_2, cmd_3])
		if cfg.verbose: log(full_cmd)
		os.system(full_cmd)
	return uniq_reads

def trim(reads, adapter3p, l, n):
	# Usage: Trims a specified number of bases from the 5' end of each read.
	# Input: List of fastq files.
	# Output: List of 5p trimmed files.
	program = cfg.home + "/bin/fasta_to_fastq.pl"
	trimmedReads = []
	for inread in reads:
		outread = inread.replace("_nodup.fasta", "_trimmed.fastq")
		trimmedReads.append(outread)
		if glob.glob(outread):
			log("5' and 3' barcode trimming already done.")
			continue
		
		cmd_1 = "perl {} {}".format(program, inread)
		cmd_2 = "fastx_clipper -n -l{} -Q33 -a {}".format(l, adapter3p)
		cmd_3 = "fastx_trimmer -f{} -Q33 > {}".format(n, outread)
		full_cmd = ' | '.join([cmd_1, cmd_2, cmd_3])
		if cfg.verbose: log(full_cmd)
		os.system(full_cmd)
	return trimmedReads

################
### MAPPING  ###
################

def run_mapping(processed_reads, exoViruses, repeat_index, endoVirus_index, trna_index, bowtie_index, num_cores):
	# Usage: Read mapping.
	# Input: Fastq files of replicate trimmed read files.
	# Output: Path to samfile for each read.

	viral_sam = []
	rep_sam = []
	endoVirus_sam = []
	trna_sam = []
	genome_sam = []
	for infastq in processed_reads:
		# NAMES 
		for v in exoViruses:
			v_mapped = infastq.replace(".fastq", "_mappedTo{}.sam".format(v))
			viral_sam.append(v_mapped)
			
		rep_mapped = infastq.replace(".fastq", "_mappedToRepeat.sam")
		rep_unmapped = infastq.replace(".fastq", "_notMappedToRepeat.fastq")

		rep_sam.append(rep_mapped)
		endoVirus_mapped = rep_unmapped.replace("_notMappedToRepeat.fastq", "_mappedToendoVirus.sam")
		endoVirus_sam.append(endoVirus_mapped)
		
		endoVirus_unmapped = rep_unmapped.replace("_notMappedToRepeat.fastq", "_notMappedToendoVirus.fastq")
		trna_mapped = endoVirus_unmapped.replace("_notMappedToendoVirus.fastq", "_mappedToTrna.sam")
		trna_sam.append(trna_mapped)

		trna_unmapped = endoVirus_unmapped.replace("_notMappedToendoVirus.fastq", "_notMappedToTrna.fastq")
		genome_star_prefix = trna_unmapped.replace("_notMappedToTrna.fastq", "")
		genome_star_output = trna_unmapped.replace("_notMappedToTrna.fastq", "Aligned.out.sam")
		genome_mapped = trna_unmapped.replace("_notMappedToTrna.fastq", "_mappedToGenome.sam")
		genome_sam.append(genome_mapped)
		
		# MAPPING
		
		if glob.glob(genome_mapped) and os.stat(genome_mapped).st_size > 0:
			log("Bowtie/STAR already done.")
			continue
					
		if exoViruses:
			log("Mapping {} to exoViruses".format(infastq))

			for v in exoViruses:    
				viral_index = cfg.home + "/docs/viral/{}".format(v)
				if 'notMappedToViral' not in infastq:  # first exoVirus mapping
					v_mapped = infastq.replace(".fastq", "_mappedTo{}.sam".format(v))
					v_unmapped = infastq.replace(".fastq", "_notMappedToViral_new.fastq".format(v))
				else:  # all subsequent exoVirus mappings
					v_mapped = infastq.replace("_notMappedToViral.fastq", "_mappedTo{}.sam".format(v))
					v_unmapped = infastq.replace("_notMappedToViral.fastq", "_notMappedToViral_new.fastq")              

				cmd = "bowtie2 -p {} -x {} {} --un {} -S {} > {} 2>&1".format(num_cores, viral_index, infastq, v_unmapped, v_mapped, v_mapped + '_stats.txt')
				if cfg.verbose: log(cmd)
				os.system(cmd)
		
				infastq = v_unmapped.replace("_notMappedToViral_new.fastq", "_notMappedToViral.fastq")
				os.system("mv {} {}".format(v_unmapped, infastq))
				
			rep_mapped = infastq.replace("_notMappedToViral.fastq", "_mappedToRepeat.sam")
			rep_unmapped = infastq.replace("_notMappedToViral.fastq", "_notMappedToRepeat.fastq")

		log("Mapping {} to repeat".format(infastq))
		cmd = "bowtie2 -p {} -x {} {} --un {} -S {} > {} 2>&1".format(num_cores, repeat_index, infastq, rep_unmapped, rep_mapped, rep_mapped + '_stats.txt')
		if cfg.verbose: log(cmd)
		os.system(cmd)
		
		log("Mapping {} to endovirus".format(infastq))
		cmd = "bowtie2 -p {} -x {} {} --un {} -S {} > {} 2>&1".format(num_cores, endoVirus_index, rep_unmapped, endoVirus_unmapped, endoVirus_mapped, endoVirus_mapped + '_stats.txt')

		if cfg.verbose: log(cmd)
		os.system(cmd)
		
		log("Mapping {} to tRNA".format(infastq))
		cmd = "bowtie2 -p {} -x {} {} --un {} -S {} > {} 2>&1".format(num_cores, trna_index, endoVirus_unmapped, trna_unmapped, trna_mapped, trna_mapped + '_stats.txt')
		if cfg.verbose: log(cmd)
		os.system(cmd)
		
		log("Mapping {} to genome".format(infastq))
		cmd = "bowtie2 -p {} -x {} {} -S {} > {} 2>&1".format(num_cores, bowtie_index, trna_unmapped, genome_mapped, genome_mapped + '_stats.txt')
		if cfg.verbose: log(cmd)
		os.system(cmd)
		
	return viral_sam, rep_sam, endoVirus_sam, trna_sam, genome_sam

def run_samtools(samfiles, flags=""):
	# Usage: Samfile processing (also takes unique mappers only)
	# Input: Sam files from Bowtie mapping.
	# Output: Sorted bedFiles.
	program = 'samtools'
	program2 = 'bamToBed'
	out_bedfiles = []
	for samfile in samfiles:
		bamfile_sort = samfile.replace('.sam','_sorted') 
		bedfile = bamfile_sort.replace('_sorted', '.bed') 
		out_bedfiles.append(bedfile)
				
		cmd_1 = "cat {} | samtools view {} -Su -F 0x4 - -o - | samtools sort - {}".format(samfile, flags, bamfile_sort)
		cmd_2 = "bamToBed -i {} > {}".format(bamfile_sort + '.bam', bedfile)
		if cfg.verbose: log(cmd_1)
		os.system(cmd_1)
		if cfg.verbose: log(cmd_2)
		os.system(cmd_2)
		
	return out_bedfiles

################
### ANALYSIS ###
################

def run_iclipro(samfiles):
	# Usage: Run iCLIPro for read overlap testing; http://www.biolab.si/iCLIPro/doc/
	# Input: mapq filtered and sorted bam files
	# Output: read overlap heatmaps and some summary stats
	
	for samfile in samfiles:
		bamfile_sort = samfile.replace('.sam','_sorted.bam') 
	
		cmd_1 = "iCLIPro -r 50 -b 300 -g \"L18:18,L19:19,L20:20,L21:21,L22:22,L23:23,L24:24,L25:25,L26:26,L27:27,L28:28,L29:29,L30:30,L31:31,L32:32,L33:33,L34:34,L35:35,L36:36,L37:37,L38:38,L39:39,R:42\" -p \"L18-R,L19-R,L20-R,L21-R,L22-R,L23-R,L24-R,L25-R,L26-R,L27-R,L28-R,L29-R,L30-R,L31-R,L32-R,L33-R,L34-R,L35-R,L36-R,L37-R,L38-R,L39-R\" -o {} {} > {}".format(bamfile_sort + "_iclipro", bamfile_sort, bamfile_sort + "_iclipro_stats.txt")
	
		if cfg.verbose: log(cmd_1)
		rt = os.system(cmd_1)
		if rt != 0: log("iCLIPro error; maybe too few reads?")
	
def trna_isotype_count(trna_samfiles, minpass, threshold):
	# Usage: get reads that are unique to only one AA of tRNA
	# retain reads that are uniquely mappable to one isotype of tRNA and make count table per isotype
	# Input: Sam files
	# Output: 1) list of de-dupped tRNA read names (to plug into file line counting later on in fig 1), 2) count table for isotypes
	
	out_trna_reads = []
	allreps_isotype_to_histo = defaultdict(lambda: defaultdict(lambda: []))
	for samfile in trna_samfiles:
		trna_isotype_reads = samfile.replace('.sam', '_trna_isotypes_reads.txt')
		trna_isotype_counts = samfile.replace('.sam', '_trna_isotypes.counts')
		trna_reads = samfile.replace('.sam', '_trnaReads.txt')
		trna_sam_txt = samfile.replace('.sam', '.samtxt')
		trna_histo = samfile.replace('.sam', '_histo.txt')
		out_trna_reads.append(trna_reads)
		
		# process the reads
		cmd_1 = "samtools view -S -F4 {} | awk -F'[\\t_]' ' {{ split($5,arr,\"-\"); print $1, arr[2] }} ' | sort -k1 -u | cut -f2 > {} ".format(samfile, trna_isotype_reads) #File with de-dupped reads per isotype; only taking mapped reads
		cmd_2 = "cat {} | awk ' $1 != prev {{ if (count==1) print line; count=0 }} {{prev=$1; line=$0; ++count }} END {{if (count==1) print }} ' | cut -d\" \" -f2 | sort | uniq -c > {} ".format(trna_isotype_reads, trna_isotype_counts) #Output the number of reads unique to each isotype
		cmd_3 = "cat {} | cut -f1 -d\" \" | sort -u > {} ".format(trna_isotype_reads, trna_reads) #Outputs dedupped list of read names mapping to tRNA
		
		os.system(cmd_1)
		os.system(cmd_2)
		os.system(cmd_3)
		
		# make histogram
		cmd_4 = "samtools view -S -F4 {} > {} ".format(samfile, trna_sam_txt) #File with de-dupped reads per isotype; only taking mapped reads
		os.system(cmd_4)
		
		isotype_to_histo = defaultdict(lambda: defaultdict(lambda: 0))
		ifile = open(trna_sam_txt, 'r')
		reader = csv.reader(ifile, 'textdialect')
		for row in reader:
			isotype = row[2].split('-')[-1]
			pos = int(row[3])
			isotype_to_histo[isotype][pos] += 1
		ifile.close()
		
		ofile = open(trna_histo, 'w')
		writer = csv.writer(ofile, 'textdialect')
		for trna in sorted(isotype_to_histo.keys()):
			row = [trna]
			for pos in range(80):
				row.append(isotype_to_histo[trna][pos])  # histogram for that particular one
				allreps_isotype_to_histo[trna][pos].append(isotype_to_histo[trna][pos])
			writer.writerow(row)
		ofile.close()       
		
		os.system("rm -f {}".format(trna_sam_txt))
	
	allreps_trna_histo = cfg.outfilepath+cfg.sampleName+'_threshold=%s'%threshold+'_repeat_allreads_trna_histo.txt'
	ofile = open(allreps_trna_histo, 'w')
	writer = csv.writer(ofile, 'textdialect')
	for trna in sorted(allreps_isotype_to_histo.keys()):
		row = [trna]
		for pos in range(80):
			nums = np.array(allreps_isotype_to_histo[trna][pos])
			if len(nums[nums >= threshold]) >= minpass: row.append(sum(nums))
			else: row.append(0)
		writer.writerow(row)
	ofile.close()       
		
	return out_trna_reads
		
##################
### PROCESSING ###
##################

#def remove_RepeatMaskerRegions(gen_bed, blacklistregions, repeatregions):
def remove_RepeatMaskerRegions(gen_bed, repeatregions):
	# Usage: Remove repeat regions from bedfile following mapping.
	# Input: .bed file after mapping (duplicates removed by samtools).
	# Output: Bedfile with repeat regions removed.

	masked=[]
	in_no_repeat_RTstop_collapsed=[]
	for bedIn in gen_bed:
		no_repeat = bedIn.replace('.bed', '_noRepeat.bed')
		repeat = bedIn.replace('.bed', '_repeat.bed')
			
		cmd_2 = "bedtools intersect -a {} -b {} -wa -v -sorted -s > {}".format(bedIn, repeatregions, no_repeat)
		cmd_3 = "bedtools intersect -a {} -b {} -wa -wb -sorted -s > {}".format(bedIn, repeatregions, repeat)

		if cfg.verbose: log(cmd_2)
		os.system(cmd_2)
		if cfg.verbose: log(cmd_3)
		os.system(cmd_3)
		masked.append(no_repeat)
		
	return masked
	

def separateStrands(mappedReads):
	# Usage: separate positive and negative strands.
	# Input: Paths to bed file from Samtools.
	# Output: Paths to bed files isolated by strand.
	negativeStrand=[]
	positiveStrand=[]
	for mapFile in mappedReads:
		with open(mapFile, 'r') as infile:
			neg_strand=mapFile.replace('.bed','_neg.bed')
			pos_strand=mapFile.replace('.bed','_pos.bed')
			negativeStrand=negativeStrand+[neg_strand]
			positiveStrand=positiveStrand+[pos_strand]
				
			neg = open(neg_strand, 'w')
			pos = open(pos_strand, 'w')
			for line in infile: 
				if str(line.strip().split('\t')[5]) == '-':
					neg.write(line)
				elif str(line.strip().split('\t')[5]) == '+':
					pos.write(line)
	return (negativeStrand,positiveStrand)

def modifyNegativeStrand(negativeStrandReads):
	# Usage: For negative stranded reads, ensure 5' position (RT stop) is listed first.
	# Input: Bed file paths to all negative stranded.
	# Output: Paths to modified bed files.
	negativeStrandEdit=[]
	for negativeRead in negativeStrandReads:
		neg_strand_edited=negativeRead.replace('_neg.bed','_negEdit.bed')
		negativeStrandEdit=negativeStrandEdit+[neg_strand_edited]
		
		neg_edit = open(neg_strand_edited, 'w')
		with open(negativeRead, 'r') as infile:
			for line in infile: 
				chrom,start,end,name,quality,strand=line.strip().split('\t')
				neg_edit.write('\t'.join((chrom, str(int(end)), str(int(end) + 29), name, quality, strand)) + '\n') # end-1 since it's [start, end)
	return negativeStrandEdit


def isolate5prime(strandedReads):
	# Usage: Isolate only the Chr, 5' position (RT stop), and strand. Shifts RT stops to correct positon based on strand. -1 position for (+) strand RT stops and +1 position for (-) strands RT stops.
	# Input: Bed file paths to strand separated reads.
	# Output: Paths RT stop files.
	RTstops=[]
	for reads in strandedReads:
		RTstop=reads.replace('.bed','_RTstop.bed')
		with open(reads, 'r') as infile:
			RTstops=RTstops+[RTstop]
			f = open(RTstop, 'w')
			for line in infile: 
				chrom,start,end,name,quality,strand=line.strip().split('\t')
				#f.write('\t'.join((chrom,start,strand))+'\n')
				#f.write('\t'.join((chrom,str(int(start)-1),strand))+'\n')
				if strand=="+":
					f.write('\t'.join((chrom,str(int(start)-1),strand))+'\n')
				else:
					f.write('\t'.join((chrom,str(int(start)+1),strand))+'\n')
	return RTstops

def fileCat(destinationFile,fileList, sort=False):
	f = open(destinationFile + "_temp", "w")
	for tempfile in fileList:
		readfile = open(tempfile, "r")
		f.write(readfile.read())
		readfile.close()
	f.close()
	
	if not sort:
		os.system("sort -k1,1 -k2,2n {} > {}".format(destinationFile + "_temp", destinationFile))
		os.system("rm -f {}".format(destinationFile + "_temp"))

def RTcounts(RTfile):
	posRT_R1=pd.DataFrame(pd.read_table(RTfile,dtype={'0': np.str, '1': np.int64, '2': np.str}, index_col=None,header=None,sep='\t',low_memory=False))
	posRT_R1.columns=['Chr','Start','Strand']
	cts=posRT_R1.groupby(['Chr','Start']).size()
	return cts

def countPassed(x, n):
	ct = 0
	for i in x:
		if i >= n: ct += 1
	return ct
	
def mergeRT(RTstopFiles, outfilename, statsfilename, minpass, threshold, strand):
	# Usage: Merge RT stops between replicates and keep only those positions that exceed threshold.
	# Input: Files with RT stops for each replicate, outfile, threshold, strand, and bases to expand around RT stop.
	# Output: None. Writes merged RT stop file.
	
	cts = [RTcounts(fn) for fn in RTstopFiles if os.stat(fn).st_size > 0]
	if not cts: 
		f = open(outfilename, 'w')
		fs = open(statsfilename, 'w')
		f.close()
		fs.close()
		return
	m = pd.concat(cts, axis=1, join='inner')
	numpass = m.apply(lambda x: countPassed(x, threshold), axis=1)
	m = m[numpass >= minpass]

	m_filter = m.copy(deep=True)
	m_filter['sum'] = m.apply(sum, axis=1)
	m_filter['mean'] = m.apply(np.mean, axis=1)
	m_filter['stdev'] = m.apply(np.std, axis=1)
	
	f = open(outfilename, 'w')
	fw = csv.writer(f, 'textdialect')
	fs = open(statsfilename, 'w')
	fsw = csv.writer(fs, 'textdialect')
	for i in m_filter.index:
		if str(i[0]).find('chr') == -1:
			chrom = 'chr' + str(i[0]) 
		else:
			chrom = i[0]
		RT = i[1]
		count = m_filter.loc[i,'sum']
		mean = str(m_filter.loc[i, 'mean'])
		stdev = str(m_filter.loc[i, 'stdev'])
		read = [chrom, str(int(RT)), str(int(RT)+1), 'CLIPread', '255', strand]
		sread = [chrom, str(int(RT)), str(int(RT)+1), 'CLIPread', '255', strand]
		sread.extend(m.loc[i])
		sread.extend([mean, stdev])
		fsw.writerow(sread)
		for _ in range(int(count)):
			fw.writerow(read)
	f.close()
	fs.close()
	
def filter_snoRNAs(negAndPosMerged, snoRNAmasker, miRNAmasker):
	# Usage: Filter snoRNA and miRNAs from protein coding reads.
	# Input: .bed file with protein coding reads.
	# Output: snoRNA and miR filtered .bed file.
	program='intersectBed'
	proteinWithoutmiRNAs = negAndPosMerged.replace('.bed','_snoRNAremoved_miRNAremoved.bed')

	cmd1 = "bedtools intersect -a {} -b {} -wa -v -s | sort -k1,1 -k2,2n".format(negAndPosMerged, snoRNAmasker)
	cmd2_1 = "bedtools intersect -a - -b {} -wa -v -s -sorted".format(miRNAmasker)
	cmd2_2 = "awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {print $1,$2,$3,$4 \"_\" NR,$5,$6}'"
	cmd2_3 = proteinWithoutmiRNAs
	cmd = cmd1 + ' | ' + cmd2_1 + ' | ' + cmd2_2 + ' > ' + cmd2_3
	if cfg.verbose: log(cmd)
	os.system(cmd)
	
	return proteinWithoutmiRNAs


def annotate_genes(sno_mirna_filtered_reads, geneStartStopRepoBed):
	# Usage: Annotate all reads that match ENSG genes; delete the rest.
	# Input: .bed file with reads that are snoRNA and miRNA filtered.
	# Output: .bed file with reads that are within ENSG genes. also annotated (BED6 + BED6)
	
	out_reads = sno_mirna_filtered_reads.replace('.bed', '_ens_annotated.bed')
	cmd = "bedtools intersect -a {} -b {} -wa -wb -s -sorted > {}".format(sno_mirna_filtered_reads, geneStartStopRepoBed, out_reads)
	os.system(cmd)
	
	with open(out_reads, 'r') as ifile, open(out_reads + '_uniq', 'w') as ofile:
		reader = csv.reader(ifile, 'textdialect')
		writer = csv.writer(ofile, 'textdialect')
		seen_reads = set()
		for row in reader:
			clip_read_name = row[3]
			if clip_read_name not in seen_reads:
				writer.writerow(row)
				seen_reads.add(clip_read_name)
		
	os.system("mv {} {}".format(out_reads + '_uniq', out_reads))
	return out_reads
	
def make_gene_list_from_annotation(gene_reads):
	outfile = gene_reads + '.geneNames'
	cmd = "cut -f10 {} > {}".format(gene_reads, outfile)
	os.system(cmd)
	return outfile
	
def getBedCenterPoints(inBed, namecol):
	# Usage: Obtain center coordinates of bedFile.
	# Input: BedFile.   
	# Output: Center coordinates returned.
	outBed = inBed.replace('.bed','_centerCoord.bed')
	with open(inBed, 'r') as ifile, open(outBed, 'w') as ofile:
		reader = csv.reader(ifile, 'textdialect')
		writer = csv.writer(ofile, 'textdialect')
		for row in reader:
			writer.writerow([row[0], int(row[1]), int(row[1]) + 1, row[namecol], row[4], row[5]])
	return outBed

def cleanBedFile(inBed, sizesFile):
	# Usage: Recover only first 6 fields from a bed file. Sorting was done in fileCat (BD 01/28/17)
	# Input: BedFile.
	# Output: bedFile with correct number of fields. Added in removal of chromosomes not in genome file.

	chrom = set()
	tempfile = inBed.replace('.bed', '_temp.bed')
	with open(sizesFile, 'r') as genome_file:
		lines = list(csv.reader(genome_file, delimiter = "\t"))
		for line in lines:
			chrom.add(line[0])	

	with open(inBed, 'r') as cur_file, open(tempfile, 'w') as new_file:
		lines = list(csv.reader(cur_file, delimiter = "\t"))
		for line in lines:
			if line[0] in chrom:
				for i in range(5):
					new_file.write(line[i] + "\t")
				new_file.write(line[5] + "\n")

	os.system("mv {} {}".format(tempfile, inBed))

	cleanedBed = inBed.replace('.bed','_cleaned_sorted.bed')    
	cmd = "cut -f1-6 {} > {}".format(inBed, cleanedBed)
	os.system(cmd)
	return cleanedBed

def makeBedGraph(negAndPosMerged, sizesFile):
	# Usage: From a bedFile, generate a plus and minus and total bedGraph and bigWig.
	# Input: BedFile.
	# Output: BedGraph file.
	cleanBed = cleanBedFile(negAndPosMerged, sizesFile)
	for strand in ["", "+", "-"]:
		word = ""
		if strand == "+": word = "_plus"
		elif strand == "-": word = "_minus"
		outname = cleanBed.replace('.bed', '{}.bedgraph'.format(word))
		outname2 = cleanBed.replace('.bed', '{}.bw'.format(word))
		if strand == "":
			cmd1 = "bedtools genomecov -bg -split -i {} -g {} > {}".format(cleanBed, sizesFile, outname)
		else:
			cmd1 = "bedtools genomecov -bg -split -i {} -g {} -strand {} > {}".format(cleanBed, sizesFile, strand, outname)
		cmd2 = cfg.home + "/bin/bedGraphToBigWig {} {} {}".format(outname, sizesFile, outname2)
		if cfg.verbose: log(cmd1)
		os.system(cmd1)
		if cfg.verbose: log(cmd2)
		os.system(cmd2)
	return cleanBed.replace('.bed', '.bedgraph')

#####################
### PARTITIONING  ###
#####################

def getReadTypes(gene_reads,pathToGeneLists):
	# Usage: Given a list of genes, return all reads for the associated genes.
	# Input: Gene list and the path to lowFDR read file.
	# Output: List of reads assocaited with the given genes.
	genelist = []
	for path in pathToGeneLists:
		outfile = path + '_reads.bed'
		cmd = "grep -F -f {} {} > {}".format(path, gene_reads, outfile)
		os.system(cmd)
		genelist=genelist+[outfile]
	return genelist

def compareLists(list1,list2,outname):
	# Usage: Compare gene lists and output matches to the file. 
	# Input: Two gene lists.
	# Output: Path file containing the matching genes.
	f=open(list1,'r')
	g=open(list2,'r')
	commonGenes=set(f.readlines()) & set(g.readlines())
	geneCategory=outname.split('.')[1]
	outputName=cfg.outfilepath+'clipGenes_'+geneCategory
	outfh=open(outputName,'w')
	for gene in commonGenes:
		outfh.write(gene)
	outfh.close()
	return outputName

def getGeneTypes(GeneList,geneAnnot):
	# Usage: Get all genes listed under each type, compare to CLIPper targets. (--clipper option)
	# Input: .bed file passed into CLIPper and the CLIPper windows file.
	# Output: Path to file containing genes of each type.
	geneTypes=[]
	for genepath in geneAnnot:
		if 'snoRNA' not in genepath:
			genes = compareLists(GeneList, genepath, os.path.split(genepath)[1])
			geneTypes.append(genes)
	return geneTypes
	
def get_gene_counts(bedFile):
	bf=pd.DataFrame(pd.read_table(bedFile,header=None))
	bf.columns=['Chr','Start','Stop','name','Q','Strand']
	bf['geneName']=bf['name'].apply(lambda x: x.split('_')[0])
	geneCounts=bf.groupby('geneName').size()
	geneCounts.sort_values(ascending=False)
	return geneCounts

def getSnoRNAreads(negAndPosMerged,snoRNAindex):
	program='intersectBed'      
	bedFile=cfg.outfilepath+'clipGenes_snoRNA_reads.bed'
	cmd = "bedtools intersect -a {} -b {} -s -wa -wb -sorted > {}".format(negAndPosMerged, snoRNAindex, bedFile)
	os.system(cmd)
	return bedFile

def countSnoRNAs(bedFile_sno):
	bf=pd.DataFrame(pd.read_table(bedFile_sno,header=None))
	bf.columns=['Chr','Start','End','name','Q','Strand','Chr_snoRNA','Start_snoRNA','Stop_snoRNA','name_snoRNA','Type','strand_snoRNA']
	geneCounts=bf.groupby('name_snoRNA').size()
	geneCounts.sort_values(ascending=False)
	return geneCounts

def countRemainingGeneTypes(remaining):
	log("countRemainingGeneTypes") #remove later
	for bedFile in remaining:
		try:
			bf=pd.DataFrame(pd.read_table(bedFile,header=None))
			bf.columns=['Chr','Start','End','ReadName','Q','Strand','winChr','winStart','winEmd','winaName','winP','winStrand']
			# *** THIS MAY DEPEND UPON THE VERSION OF CLIPPER USED ***
			bf['geneName']=bf['winaName'].apply(lambda x: x.split('_')[0])
			geneCounts=bf.groupby('geneName').size()
			geneCounts.sort_values(ascending=False) 
						
			head,fname=os.path.split(bedFile)
			geneType=fname.split("_")[1]
			cfg.outfilepathToSave=cfg.outfilepath+'/PlotData_ReadsPerGene_%s.txt'%geneType
			geneCounts.to_csv(cfg.outfilepathToSave)
			
		except ValueError:
			log("No reads in {}".format(bedFile))

#####################
### INTRON/UTRS   ###
#####################           

def extract_regions(bedIn, fivePUTRBed, threePUTRBed, exonBed):
	# Usage: Extract all region specific reads from the input file.
	# Input: .bed files
	# Output: Mutually exclusive partitions of the input file.
	
	exonreads = bedIn.replace('.bed', '_exons.bed')
	intronreads = bedIn.replace('.bed', '_introns.bed')
	cmd1 = "bedtools intersect -a {} -b {} -wa -s -u > {}".format(bedIn, exonBed, exonreads)
	cmd2 = "bedtools intersect -a {} -b {} -wa -s -v > {}".format(bedIn, exonBed, intronreads)
	os.system(cmd1)  # Exons
	os.system(cmd2)  # Introns
	
	fivePreads = bedIn.replace('.bed', '_5p.bed')
	threePreads = bedIn.replace('.bed', '_3p.bed')
	cmd1 = "bedtools intersect -a {} -b {} -wa -s -u > {}".format(exonreads, fivePUTRBed, fivePreads)
	cmd2 = "bedtools intersect -a {} -b {} -wa -s -u > {}".format(exonreads, threePUTRBed, threePreads)
	os.system(cmd1)  # 5'
	os.system(cmd2)  # 3'
	
	cdsreads = bedIn.replace('.bed', '_cds.bed')
	cmd1 = "cat {} {} | bedtools intersect -a {} -b stdin -wa -s -v > {}".format(fivePUTRBed, threePUTRBed, exonreads, cdsreads)
	os.system(cmd1)  # CDS
	
	return (exonreads, intronreads, fivePreads, threePreads, cdsreads)

def read_readspergene_file(data_dict, pos, fn):
	with open(fn, 'r') as ifile:
		for line in ifile:
			[gene, cts] = line.rstrip().split(',')
			data_dict[gene][pos] = int(cts)
	return data_dict
	
def gene_binding_by_region():
	fivePfile = cfg.outfilepath + '/PlotData_ReadsPerGene_5pUTR.txt'
	threePfile = cfg.outfilepath + '/PlotData_ReadsPerGene_3pUTR.txt'
	CDSfile = cfg.outfilepath + '/PlotData_ReadsPerGene_CDS.txt'
	exonfile = cfg.outfilepath + '/PlotData_ReadsPerGene_Exons.txt'
	intronfile = cfg.outfilepath + '/PlotData_ReadsPerGene_Introns.txt'
	allfile = cfg.outfilepath + '/PlotData_ReadsPerGene_proteinCoding.txt'
	
	gene_to_counts = defaultdict(lambda: [0,0,0,0,0,0])  # all, exon, intron, 5P, CDS, 3P
	if glob.glob(allfile):
		gene_to_counts = read_readspergene_file(gene_to_counts, 0, allfile)
	if glob.glob(exonfile):
		gene_to_counts = read_readspergene_file(gene_to_counts, 1, exonfile)
	if glob.glob(intronfile):
		gene_to_counts = read_readspergene_file(gene_to_counts, 2, intronfile)
	if glob.glob(fivePfile):
		gene_to_counts = read_readspergene_file(gene_to_counts, 3, fivePfile)
	if glob.glob(CDSfile):
		gene_to_counts = read_readspergene_file(gene_to_counts, 4, CDSfile)
	if glob.glob(threePfile):
		gene_to_counts = read_readspergene_file(gene_to_counts, 5, threePfile)
	
	ofn = cfg.outfilepath + '/PlotData_proteinCoding_byRegion.txt'
	with open(ofn, 'w') as ofile:
		writer = csv.writer(ofile)
		writer.writerow(['Gene', 'All', 'Exonic', 'Intronic', '5P', 'CDS', '3P'])
		for gene in gene_to_counts:
			row = [gene]
			row.extend(gene_to_counts[gene])
			writer.writerow(row)
	
#####################
###  METAGENES    ###
#####################   

def getGeneStartStop(bedFile,geneRef):
	try:
		bf=pd.DataFrame(pd.read_table(bedFile,header=None))
		bf.columns=['Chr','Start','End','ReadName','Q','Strand','winChr','winStart','winEmd','winaName','winP','winStrand']
		bf['geneName']=bf['winaName'].apply(lambda x: x.split('_')[0].split('.')[0])
		merge=pd.merge(geneRef,bf,left_on='Ensembl Gene ID',right_on='geneName')
		ncRNA_startStop=merge[['Ensembl Gene ID','Gene Start (bp)','Gene End (bp)','Start','End','Strand']]
		outfilepathToSave=bedFile.replace(".bed",".geneStartStop")
		ncRNA_startStop.to_csv(outfilepathToSave)
	except ValueError:
		log("No reads in {}".format(bedFile))
		
def makeRepeatAnnotation(repeatGenomeBuild,repeatAnnotation):
	repeat_genome=np.genfromtxt(repeatGenomeBuild,dtype='string')
	repeat_genome_bases=repeat_genome[1]
	repeat_genome_size=len(repeat_genome[1])
	repeatAnnotDF=pd.DataFrame(pd.read_table(repeatAnnotation,header=None))
	repeatAnnotDF.columns=['Name','Length','IndexStart','IndexEnd']
	repeatAnnotDF['End_for_extraction']=repeatAnnotDF['IndexEnd']+1 # Python list extraction is not end index inclusive; to extract sequence, use end + 1. 
	return (repeat_genome_bases,repeatAnnotDF)

############
## PLOT 1 ##
############

def lineCount(filename):
	if cfg.verbose: log(filename)
	if filename[-3:] == '.gz':
		cmd_1 = "gunzip -c {}".format(filename)
	else:
		cmd_1 = "cat {}".format(filename)
	cmd = cmd_1 + ' | wc -l'
	return int(commands.getstatusoutput(cmd)[1])

def plot_ReadAccounting(nsamp, reads, threshold_nr, index_tag):

	bases = [os.path.splitext(os.path.basename(reads[i]))[0] for i in range(nsamp)]
	bases = [base.replace('.fastq.gz', '').replace('.fastq', '') for base in bases]
	bases = [cfg.outfilepath + base for base in bases]
	
	filesToCount = reads
	filesToCount.extend([base + '_trimmed.fastq' for base in bases])
	filesToCount.extend([base + '_trimmed_mappedToRepeat.bed' for base in bases])
	filesToCount.extend([base + '_trimmed_mappedToTrna.bed' for base in bases])
	filesToCount.extend([base + '_trimmed_mappedToGenome.bed' for base in bases])
	filesToCount.extend([base + '_trimmed_mappedToGenome_noRepeat.bed' for base in bases])
	filesToCount.extend([base + '_trimmed_mappedToendoVirus.bed' for base in bases])
	filesToCount.extend([base + '_trimmed_mappedToDV.bed' for base in bases])
	filesToCount.extend([base + '_trimmed_mappedToZV.bed' for base in bases])
	filesToCount.extend([base + '_trimmed_mappedToHCV_JFH1.bed' for base in bases])
	
	fileNameBases = ['Raw', 'No dupes', 'Repeat Mapped', 'tRNA Mapped', 'Genome Mapped', 'Blacklist Masked', 'EndoVirus', 'DV', 'ZV', 'HCV_JFH1']
	fileNames = [b + " (R{})".format(i + 1) for b in fileNameBases for i in range(nsamp)]
	
	counts = []
	counter = 0
	for fileString in filesToCount:
		temp = lineCount(fileString)
		if '.fastq' in fileString:
			temp = temp/4 # Fastq files
		counts = counts+[temp]
		counter += 1

	ind = np.arange(len(counts)) + 0.5
	plt.barh(ind,list(reversed(np.log10(np.array(counts)))),align='center',color='blue')
	plt.xlabel('log10(Counts per file)',fontsize=5)
	locs,pltlabels = plt.xticks(fontsize=5)
	plt.setp(pltlabels, rotation=90, fontsize=5)
	plt.yticks(ind,list(reversed(fileNames)),fontsize=5)
	plt.tick_params(axis='yticks',labelsize=5) 
	ax=plt.gca()
	for line in ax.get_yticklines():
		line.set_markersize(0)
	plt.title('Read counts',fontsize=5)
	
	if cfg.verbose: log(str(fileNames))
	if cfg.verbose: log(str(counts))
	readDF=pd.DataFrame()
	readDF['File_name'] = fileNames
	readDF['Reads_per_file'] = counts
	outfilepathToSave=cfg.outfilepath + '/PlotData_ReadsPerPipeFile.txt'
	readDF.to_csv(outfilepathToSave)

def plot_BoundGeneTypes():
	record=pd.DataFrame()   
	
	# Exclude specific files
	geneListToPlot = [] 
	for f in glob.glob(cfg.outfilepath+'PlotData_ReadsPerGene_*'):
		want = True
		for g in ['5pUTR', '3pUTR', 'Introns', 'Exons', 'CDS']:
			if g in f: want = False
		if want: geneListToPlot.append(f)
		
	for boundGenes in geneListToPlot:
		glist=pd.read_csv(boundGenes,header=None)
		glist.columns=['GeneName','Count']
		gName=boundGenes.split('_')[-1]
		record.loc[gName,'genesBound']=glist.shape[0]
		record.loc[gName,'totalReads']=glist['Count'].sum()
	record.sort_values('genesBound',inplace=True)
	outfilepathToSave=cfg.outfilepath + '/PlotData_ReadAndGeneCountsPerGenetype.txt'
	record.to_csv(outfilepathToSave)
	ind = np.arange(record.shape[0]) + 0.5
	plt.bar(ind,record['genesBound'],align='center',color='blue')
	locs,pltlabels = plt.yticks(fontsize=5)
	locs,pltlabels = plt.xticks(ind,record.index,fontsize=5)
	plt.setp(pltlabels, rotation=90, fontsize=5)
	plt.tick_params(axis='xticks',labelsize=5) 
	ax=plt.gca()
	for line in ax.get_xticklines():
		line.set_markersize(0)
	plt.ylabel('Number of genes bound',fontsize=5)
	plt.tick_params(axis='yticks',labelsize=5)
	plt.title('Bound genes by class',fontsize=5)
	
def readUTRfile(path):
	geneCounts=pd.read_csv(path,header=None)
	geneCounts.columns=['Gene_name','Count']
	return geneCounts

def plot_readsBymRNAregion(ax): 
	fivePfile = cfg.outfilepath+'/PlotData_ReadsPerGene_5pUTR.txt'
	threePfile = cfg.outfilepath+'/PlotData_ReadsPerGene_3pUTR.txt'
	CDSfile = cfg.outfilepath+'/PlotData_ReadsPerGene_CDS.txt'
	pc_5pReads = readUTRfile(fivePfile)['Count'].sum() if glob.glob(fivePfile) else 0
	pc_3pReads = readUTRfile(threePfile)['Count'].sum() if glob.glob(threePfile) else 0
	pc_CDSReads = readUTRfile(CDSfile)['Count'].sum() if glob.glob(CDSfile) else 0
	non_intronic=pc_5pReads+pc_3pReads+pc_CDSReads
	allProteinCoding=cfg.outfilepath +'clipGenes_proteinCoding_reads_centerCoord.bed'
	all_pc=pd.DataFrame(pd.read_table(allProteinCoding,header=None))
	pc_allReads=all_pc.shape[0]
	v=[float(pc_allReads-non_intronic)/pc_allReads,float(pc_5pReads)/pc_allReads,float(pc_CDSReads)/pc_allReads,float(pc_3pReads)/pc_allReads]
	pie_wedges=ax.pie(v,labels=["Intronic","5p UTR","CDS","3pUTR"],labeldistance=1.1,autopct='%1.1f%%')
	plt.rcParams['font.size']=2
	for wedge in pie_wedges[0]:
		wedge.set_edgecolor('black')
		wedge.set_lw(1)
		
def plot_figure1(nsamp, reads, threshold_nr, index_tag, exoViruses, reads_files):
	plt.subplot(2,2,1) 
	plot_ReadAccounting(nsamp, reads, threshold_nr, index_tag)
	ax = plt.subplot(2,2,2)
	plot_readsBymRNAregion(ax)
	plt.subplot(2,2,3)
	plot_BoundGeneTypes()
	plt.subplot(2,2,4)
	make_reads_by_type_pie2(nsamp, reads, threshold_nr, index_tag)
	plt.tight_layout()
	
############
## PLOT 2 ##
############


def plot_mRNAgeneBodyDist():
	averageGraph=cfg.outfilepath+'clipGenes_proteinCoding_reads_centerCoord_UTRs_scaled_cds200_abt0_averageGraph.txt'
	hmap=pd.DataFrame(pd.read_table(averageGraph,header=None,skiprows=1))
	hmap=hmap.set_index(0)
	avgTrace=hmap.loc['treat',:]
	plt.plot(avgTrace,color='blue',linewidth='2')
	plt.vlines(200,0,np.max(avgTrace),linestyles='dashed')
	plt.vlines(400,0,np.max(avgTrace),linestyles='dashed')
	plt.ylim(0,np.max(avgTrace))
	plt.tick_params(axis='x',labelbottom='off') 
	plt.xlabel('mRNA gene body (5pUTR, CDS, 3pUTR)')
	plt.ylabel('Read density')
	plt.tick_params(axis='y',labelsize=5) 
	plt.title('CLIP signal across average mRNA transcript.',fontsize=5)

def convertENBLids(enst_name, ensemblGeneAnnot):
	ensg_name=ensemblGeneAnnot.loc[enst_name,'name2']
	return ensg_name

def getUTRbindingProfile(utr,hmap_m):
	if utr=='5p':
		ix=(hmap_m[range(201,601)].sum(axis=1)==0)&(hmap_m[range(1,201)].sum(axis=1)>0)
		screen=readUTRfile(cfg.outfilepath+'/PlotData_ReadsPerGene_5pUTR.txt')
	elif utr=='3p':
		ix=(hmap_m[range(1,401)].sum(axis=1)==0)&(hmap_m[range(401,601)].sum(axis=1)>0)
		screen=readUTRfile(cfg.outfilepath+'/PlotData_ReadsPerGene_3pUTR.txt')
	elif utr=='5p3p':
		ix=(hmap_m[range(201,401)].sum(axis=1)==0)&(hmap_m[range(401,601)].sum(axis=1)>0)&(hmap_m[range(1,201)].sum(axis=1)>0)
		screen=readUTRfile(cfg.outfilepath+'/PlotData_ReadsPerGene_5p3pUTR.txt')
	else:  # utr=='CDS'
		ix=(hmap_m[range(1,201)].sum(axis=1)==0)&(hmap_m[range(401,601)].sum(axis=1)==0)&(hmap_m[range(201,401)].sum(axis=1)>0)
		screen=readUTRfile(cfg.outfilepath+'/PlotData_ReadsPerGene_CDS.txt')
		
	# Ensure all genes are also identified in pre-allocated gene lists.
	screen['Gene_name'] = screen['Gene_name'].map(lambda x: x.split('.')[0]) # remove transcript number from gene name
	hmap_m_utrSpec=hmap_m.ix[ix,:]
	hmap_m_utrSpec_filter=pd.merge(hmap_m_utrSpec,screen,left_on='ENSG_ID',right_on='Gene_name',how='inner')
	sums=hmap_m_utrSpec_filter[range(1,601)].sum(axis=1)
	hmap_m_utrSpec_filter=hmap_m_utrSpec_filter.loc[np.argsort(sums),:]
	return hmap_m_utrSpec_filter

def plot_geneBodyPartition(ensemblGeneAnnot):
	treatMatrix=cfg.outfilepath+'clipGenes_proteinCoding_reads_centerCoord_UTRs_scaled_cds200_abt0_treatmatrix.txt'
	hmap=pd.DataFrame(pd.read_table(treatMatrix,header=None,skiprows=1))
	
	# Ensure genes recovered from this analysis are independently identified 
	hmap['ENSG_ID']=hmap.ix[:,0].apply(convertENBLids, args=(ensemblGeneAnnot,))
	bound_pc = cfg.outfilepath+'clipGenes_proteinCoding'
	pc_genes=pd.DataFrame(pd.read_table(bound_pc,header=None,))
	pc_genes.columns=['ENSG_ID']
	hmap_m=pd.merge(hmap,pc_genes,left_on='ENSG_ID',right_on='ENSG_ID',how='inner') 
	
	# Isolate intronic bound genes.
	tosave=cfg.outfilepath+'PlotData_ExclusiveBound_Intronic.tdt' 
	intronicBoundGenes=list(set(pc_genes['ENSG_ID'])-set(hmap_m['ENSG_ID']))
	np.savetxt(tosave,np.array(intronicBoundGenes),fmt="%s")
	
	# UTR specific genes.
	geneTypes=['5p','cds','3p','5p3p'] 
	depth=50
	for i in range(0,4):    
		if i == 0: name = "5pUTR"
		elif i == 1: name = "CDS"
		elif i == 2: name = "3pUTR"
		else: name = '5p3pUTR'
		fn = cfg.outfilepath+'/PlotData_ReadsPerGene_'+name+".txt"
		if not glob.glob(fn): continue
		
		plt.subplot2grid((2,4),(1,i),colspan=1)
		utrMatrix=getUTRbindingProfile(geneTypes[i],hmap_m)
		tosave=cfg.outfilepath+'PlotData_ExclusiveBound_%s.txt'%geneTypes[i] 
		np.savetxt(tosave,np.array(list(set(utrMatrix['ENSG_ID']))),fmt="%s")
		dataToPlot=utrMatrix[range(1,601)]
		
		if not dataToPlot.empty:
			p=plt.pcolormesh(np.array(dataToPlot)[-depth:-1,:],cmap='Blues')
			plt.title(geneTypes[i],fontsize=5)
			plt.vlines(200,0,depth,linestyles='dashed')
			plt.vlines(400,0,depth,linestyles='dashed')
			plt.tick_params(axis='x',labelbottom='off') 
			plt.tick_params(axis='y',labelleft='off') 
			plt.ylim(0,depth)
			plt.ylabel('Ranked genes (highest on bottom)',fontsize=5)
			plt.xticks(visible=False)
			plt.yticks(visible=False)
			plt.title('%s specific genes: %s'%(geneTypes[i],np.unique(utrMatrix['ENSG_ID']).shape[0]),fontsize=5)
	 
############
## PLOT 3 ##
############

def read_csv_by_line(fn):
	RTpositions_plus = []
	RTpositions_minus = []
	start = 0
	end = 0
	with open(fn, 'r') as ifile:
		reader = csv.reader(ifile)
		reader.next()
		for row in reader:
			if row[6] == '+': RTpositions_plus.append(int(row[7]))
			else: RTpositions_minus.append(int(row[7]))
			
		ifile.seek(0)
		reader.next()  # skip header
		row = reader.next()
		# added "-1"
		start = int(row[9])-1
		end = int(row[8])
	return np.array(RTpositions_plus), np.array(RTpositions_minus), start, end
	
def get_histogram(RTpositions, bins, repeat_genome_bases, start, end):
	hist, bins = np.histogram(RTpositions, bins=bins)
	center = (bins[:-1] + bins[1:])/2
	
	# Normalize
	histPlot = np.array(hist, dtype=float)
	if len(RTpositions) > 0:
		histPlot = np.array(histPlot/float(len(RTpositions)), dtype=float)
		
	# Make sure same length
	sequence = repeat_genome_bases[start:end+1]
	hist = hist[:len(sequence)]
	histPlot = hist[:len(sequence)]
	center = center[:len(sequence)]
	
	return hist, histPlot, sequence, center 


def plot_RT_stops(pos_reads, neg_reads, name):
	pos_neg_map = {}
	for read in pos_reads:
		if read in pos_neg_map:
			pos_neg_map[read][0] += 1
		else:
			pos_neg_map[read] = [0 for x in range(2)]
			pos_neg_map[read][0] = 1
			pos_neg_map[read][1] = 0

	for read in neg_reads:
		if read in pos_neg_map:
			pos_neg_map[read][1] += 1
		else:
			pos_neg_map[read] = [0 for x in range(2)]
			pos_neg_map[read][0] = 0
			pos_neg_map[read][1] = 1

	if pos_reads and neg_reads:
		num_bins = max(max(pos_reads), max(neg_reads)) 
	elif pos_reads:
		num_bins = max(pos_reads)
	elif neg_reads:
		num_bins = max(neg_reads)
	else:
		num_bins = 0
	num_bins += 200
	bins_initial = range(num_bins)
	
	hist_pos, bins_pos = np.histogram(pos_reads, bins = bins_initial)
	center_pos = (bins_pos[:-1] + bins_pos[1:])/2
	histPlot_pos = np.array(hist_pos, dtype=float)

	hist_neg, bins_neg = np.histogram(neg_reads, bins = bins_initial)
	center_neg = (bins_neg[:-1] + bins_neg[1:])/2

	histPlot_neg = np.array(hist_neg * -1, dtype=float)
	neg_min = min(histPlot_neg)
	pos_max = max(histPlot_pos)
	if abs(neg_min) > 5 * pos_max:
		pos_max = abs(neg_min) / 5
	elif pos_max > 5 * abs(neg_min):
		neg_min = pos_max / -5
	if pos_max < 1:
		pos_max = 1
	if neg_min > -1:
		neg_min = -1

	plt.bar(center_pos.tolist(), histPlot_pos.tolist(), align = 'center', width = 1, color = 'blue',
		edgecolor = 'blue', alpha = 0.5)
	plt.bar(center_neg.tolist(), histPlot_neg.tolist(), align = 'center', width = 1, color = 'red', 
		edgecolor = 'red', alpha = 0.5)

	plt.xticks(rotation='vertical')#,verticalalignment='bottom')
	plt.xlim(0, num_bins)
	plt.ylim(1.2 * neg_min, 1.2 * pos_max)
	name_=name.replace(".txt", "")
	plt.title("%s" % name_)


def plot_repeatRNA(repeatGenomeBuild):
	repeat_genome = np.genfromtxt(repeatGenomeBuild,dtype='string')
	repeat_genome_bases = repeat_genome[1]
	repFiles=glob.glob(cfg.outfilepath + '/PlotData_RepeatRNAreads_*')
	
	plotDim=math.ceil(math.sqrt(len(repFiles)))
	i=0
	import pickle
	for path in repFiles:
		name=path.split('RepeatRNAreads_')[-1]
		try:
			(RTpositions_plus, RTpositions_minus, start, end) = read_csv_by_line(path)
		except:
			log("No reads for repeatRNA {}".format(name))
			continue
			
		# Histogram of RT stops across gene body
		bins = range(start, end + 2, 1)
		width = 0.7*(bins[1]-bins[0])
		(hist_plus, histPlot_plus, sequence_plus, center_plus) = get_histogram(RTpositions_plus, bins, repeat_genome_bases, start, end)
		(hist_minus, histPlot_minus, sequence_minus, center_minus) = get_histogram(RTpositions_minus, bins, repeat_genome_bases, start, end)
	
		# Subplot
		plt.subplot(plotDim, plotDim, i+1)
		plot_RT_stops(RTpositions_plus.tolist(), RTpositions_minus.tolist(), name)

		# Record data
		storageDF = pd.DataFrame()
		storageDF['Sequence'] = pd.Series(list(sequence_plus))
		storageDF['RT_stops'] = np.array(list(hist_plus))
		storageDF['RT_stops_minus'] = np.array(list(hist_minus))
		outfilepathToSave = cfg.outfilepath + '/PlotData_RepeatRNAHist_%s'%name
		storageDF.to_csv(outfilepathToSave)
		i += 1

	plt.figtext(0.42, 1.02, 'RT Stops (+ blue, - red)', fontdict=None, fontsize='9')
	plt.figtext(0.46, 0.55, 'Genomic Coordinate', fontdict=None, fontsize='7')
	plt.figtext(-0.01, 0.815, 'No. of Hits', fontdict=None, fontsize='7', rotation='vertical')
	
############
## PLOT 4 ##
############

def plot_rDNA(start18s, end18s, start5s, end5s, start28s, end28s, rRNAend):
	plt.subplot2grid((3,3),(0,0),colspan=3)
	name='rDNA'
	rDNA=glob.glob(cfg.outfilepath + 'PlotData_RepeatRNAreads_rDNA.txt')
	
	(RTpositions, RTpositions_minus, start, end) = read_csv_by_line(rDNA[0])

	bins=range(start,end+2,1)
	hist,bins=np.histogram(RTpositions,bins=bins)
	width=0.7*(bins[1]-bins[0])
	center=(bins[:-1]+bins[1:])/2
	histPlot=np.array(hist,dtype=float)
	histPlot=np.array(histPlot/float(len(RTpositions)),dtype=float)
	plt.bar(center,histPlot,align='center',width=width,color='blue',alpha=0.45)
	plt.tick_params(axis='x',labelsize=2.5) 
	plt.tick_params(axis='y',labelsize=2.5)  
	plt.title('RT stops for %s: %s'%(name,len(RTpositions)),fontsize=5)
	plt.xlim(start,end)    
	
	# Features of rDNA with respect to start of the bowtie index (index=0)
	rRNAstart=start
	plt.axvspan(start18s+rRNAstart,end18s+rRNAstart,facecolor='g',alpha=0.5)
	plt.axvspan(start5s+rRNAstart,end5s+rRNAstart,facecolor='r',alpha=0.5)
	plt.axvspan(start28s+rRNAstart,end28s+rRNAstart,facecolor='b',alpha=0.5)
	
	# Generate histogram for transcribed region
	plt.subplot2grid((3,3),(1,0),colspan=3)
	datarDNAOnly=RTpositions-start
	bins=range((start-start),(end-start+2),1)
	hist,bins=np.histogram(datarDNAOnly,bins=bins)
	width=0.7*(bins[1]-bins[0])
	center=(bins[:-1] + bins[1:])/2
	histPlot=np.array(hist,dtype=float)
	histPlot=np.array(histPlot/float(len(RTpositions)),dtype=float)
	plt.bar(center,histPlot,align='center',width=width,color='blue',alpha=0.45)
	plt.tick_params(axis='x',labelsize=2.5) 
	plt.tick_params(axis='y',labelsize=2.5)  
	plt.xlabel('rRNA locus position (bin=1 base)',fontsize=5)
	plt.ylabel('Normalized RT stop / bin',fontsize=2.5)
	plt.axvspan(start18s,end18s,facecolor='g',alpha=0.5)
	plt.axvspan(start5s,end5s,facecolor='r',alpha=0.5)
	plt.axvspan(start28s,end28s,facecolor='b',alpha=0.5)
	plt.xlim(0,rRNAend)
	
	# Individual regions 
	plt.subplot2grid((3,3),(2,0),colspan=1)
	plt.bar(center,histPlot,align='center',width=width,color='green',alpha=0.75)
	plt.xlim(start18s,end18s)
	plt.xlabel('rRNA locus position (bin=1 base)',fontsize=5)
	plt.ylabel('Normalized RT stop / bin',fontsize=2.5)
	plt.tick_params(axis='x',labelsize=5) 
	plt.tick_params(axis='y',labelsize=5) 
	plt.title('18s Region',fontsize=5)
	plt.subplot2grid((3,3),(2,1),colspan=1)
	plt.bar(center,histPlot,align='center',width=width,color='red',alpha=0.75)
	plt.xlim(start5s,end5s)
	plt.xlabel('rRNA locus position (bin=1 base)',fontsize=5)
	plt.tick_params(axis='x',labelsize=5) 
	plt.tick_params(axis='y',labelsize=5) 
	plt.title('5.8s Region',fontsize=5)
	plt.subplot2grid((3,3),(2,2),colspan=1)
	plt.bar(center,histPlot,align='center',width=width,color='blue',alpha=0.75)
	plt.xlim(start28s,end28s)
	plt.xlabel('rRNA locus position (bin=1 base)',fontsize=5)
	plt.tick_params(axis='x',labelsize=5) 
	plt.tick_params(axis='y',labelsize=5)  
	plt.title('28s Region',fontsize=5)
	plt.tight_layout()

############
## PLOT 5 ##
############

def getBindingFrac(type_specific):
	# 5' position on the negative strand is snoRNA stop coordinate.
	neg_data=type_specific[type_specific['strand_snoRNA']=='-']
	neg_data['diff']=np.abs(neg_data['Stop_snoRNA']-neg_data['Start']) 
	neg_data['frac']=neg_data['diff']/(neg_data['Stop_snoRNA']-neg_data['Start_snoRNA'])
	# 5' position on the positive strand is snoRNA start coordinate.
	pos_data=type_specific[type_specific['strand_snoRNA']=='+']
	pos_data['diff']=np.abs(pos_data['Start_snoRNA']-pos_data['Start'])
	pos_data['frac']=pos_data['diff']/(pos_data['Stop_snoRNA']-pos_data['Start_snoRNA'])
	DF_snoProfile=pd.concat([neg_data,pos_data])
	return DF_snoProfile

def plot_snorna(snorna_file):
	bf_sno=pd.read_table(snorna_file,header=None)
	bf_sno.columns=['Chr','Start','End','name','Q','Strand','Chr_snoRNA','Start_snoRNA','Stop_snoRNA','name_snoRNA','Type','strand_snoRNA']
	snoTypes=pd.DataFrame(bf_sno.groupby('Type').size())
	snoTypes.columns=['Reads']
	snoTypes['Fraction']=snoTypes['Reads']/sum(snoTypes['Reads'])
	outfilepathToSave=cfg.outfilepath+'/PlotData_readsPerSnoRNAType.txt'
	snoTypes.to_csv(outfilepathToSave)

	snoTypesAndGenes=pd.DataFrame(bf_sno.groupby(['Type','name_snoRNA']).size())
	snoTypesAndGenes.columns=['Count_per_gene']
	outfilepathToSave=cfg.outfilepath+'/PlotData_geneStatsPerSnoRNAType.txt'
	snoTypesAndGenes.to_csv(outfilepathToSave)
		  
	ax=plt.subplot(2,2,1)
	pie_wedges=ax.pie(snoTypes['Fraction'],labels=snoTypes.index,labeldistance=1.1,autopct='%1.1f%%')
	plt.rcParams['font.size']=5
	for wedge in pie_wedges[0]:
		wedge.set_edgecolor('black')
		wedge.set_lw(1)

	i=2
	for sType in set(bf_sno['Type']):
		type_specific=bf_sno[bf_sno['Type']==sType]
		sno_profile=getBindingFrac(type_specific)
		
		if sType=='C':
			title="C/D_box"
		elif sType=='H':
			title="H/ACA_box"
		elif sType=='s':
			title="scaRNA"
		else:
			title="Unknown"
		
		outfilepathToSave=cfg.outfilepath+'/PlotData_snoRNAReadDist_%s.txt'%sType
		sno_profile.to_csv(outfilepathToSave)
		
		plt.subplot(2,3,i)
		bins=np.arange(0,1,0.01)
		hist,bins=np.histogram(sno_profile['frac'],bins=bins)
		hist=np.array(hist/float(sno_profile['frac'].shape[0]),dtype=float)
		width=0.7*(bins[1]-bins[0])
		center=(bins[:-1] + bins[1:])/2
		plt.bar(center,hist,align='center',width=width,color='blue',alpha=0.75)
		plt.tick_params(axis='x',labelsize=5) 
		plt.tick_params(axis='y',labelsize=5)  
		plt.xlabel('Fraction of gene body (5p - 3p)',fontsize=5)
		plt.title('Binding profile for %s'%title,fontsize=5)
		i+=1



############
## PLOT 6 ##
############

def getncRNABindingFrac(type_specific):
	# 5' position on the negative strand is snoRNA stop coordinate.
	neg_data=type_specific[type_specific['Strand']=='-']
	neg_data['diff']=np.abs(neg_data['Gene End (bp)']-neg_data['RT_stop']) 
	neg_data['frac']=neg_data['diff']/(neg_data['Gene End (bp)']-neg_data['Gene Start (bp)'])
	# 5' position on the positive strand is snoRNA start coordinate.
	pos_data=type_specific[type_specific['Strand']=='+']
	pos_data['diff']=np.abs(pos_data['Gene Start (bp)']-pos_data['RT_stop'])
	pos_data['frac']=pos_data['diff']/(pos_data['Gene End (bp)']-pos_data['Gene Start (bp)'])
	DF_ncRNAProfile=pd.concat([neg_data,pos_data])
	return DF_ncRNAProfile  
	
def plot_ncrnas(st_stopFiles, expand):
	plotDim=math.ceil(math.sqrt(len(st_stopFiles)))
	i=1
	for st_file in st_stopFiles:
		name=st_file.split('clipGenes_')[1].split('_reads')[0]
		tmp=pd.read_csv(st_file)
		tmp['RT_stop']=tmp['Start']+expand
		tmp_profile=getncRNABindingFrac(tmp)
		plt.subplot(plotDim,plotDim,i)
		bins=np.arange(0,1,0.01)
		hist,bins=np.histogram(tmp_profile['frac'],bins=bins)
		hist=np.array(hist/float(tmp_profile['frac'].shape[0]),dtype=float)
		width=0.7*(bins[1]-bins[0])
		center=(bins[:-1] + bins[1:])/2
		plt.bar(center,hist,align='center',width=width,color='blue',alpha=0.75)
		plt.tick_params(axis='x',labelsize=5) 
		plt.tick_params(axis='y',labelsize=5)  
		plt.xlabel('Fraction of gene body (5p - 3p)',fontsize=5)
		plt.title('Binding profile for %s'%name,fontsize=5)
		i+=1
		

def viral_RT_stops(pathtofile, filename, exoV_index):
	pos_reads = []
	neg_reads = []
	filename_split = filename.split("=")
	ID = "7_" + filename_split[1].split("_")[2]
	name = filename_split[0].strip("_threshold") + "_" + filename_split[1].split("_")[2] + ".txt"

	with open(pathtofile, 'r') as f:
		lines = [line.split("\t") for line in f]
		num_lines = len(lines)
	for i in range(num_lines):
		if lines[i][5][0] == '+':
			pos_reads.append(int(lines[i][2]))
		else:
			neg_reads.append(int(lines[i][2]))
	plot_RT_stops(pos_reads, neg_reads, name)
	plt.savefig(cfg.outfilepath + "Figure%s.png" % ID, format='png',bbox_inches='tight',dpi=300,pad_inches=0.5)
	plt.cla()
	plt.clf()

	# Histogram of RT stops across gene body
	log("exoV_index[0]: " + exoV_index[0])
	exoV_genome = np.genfromtxt(exoV_index[0],dtype='string')
	exoV_genome_bases = exoV_genome[1]
	start=0
	end=len(exoV_genome_bases)
	bins = range(start, end + 2, 1)
	width = 0.7*(bins[1]-bins[0])
	(hist_plus, histPlot_plus, sequence_plus, center_plus) = get_histogram(pos_reads, bins, exoV_genome_bases, start, end)
	(hist_minus, histPlot_minus, sequence_minus, center_minus) = get_histogram(neg_reads, bins, exoV_genome_bases, start, end)
	# Record data
	name = filename_split[1].split("_")[2] + ".txt"
	storageDF = pd.DataFrame()
	storageDF['Sequence'] = pd.Series(list(sequence_plus))
	storageDF['RT_stops'] = np.array(list(hist_plus))
	storageDF['RT_stops_minus'] = np.array(list(hist_minus))
	outfilepathToSave = cfg.outfilepath + '/PlotData_RepeatRNAHist_%s'%name
	storageDF.to_csv(outfilepathToSave)

def make_reads_by_region_pie():
	reads_by_type = {}
	read_counts_by_type = {}
	gene_types = ["5pUTR", "3pUTR", "CDS", "Introns"]

	for gene_type in gene_types:
		reads_by_type[gene_type] = set()
		read_counts_by_type[gene_type] = {}
		name = cfg.outfilepath + "PlotData_ReadsPerGene_%s.txt" % gene_type
		if not glob.glob(name):
			continue
		with open(name, 'r') as f:
			lines = list(csv.reader(f))
			for line in lines:
				reads_by_type[gene_type].add(line[0])
				read_counts_by_type[gene_type][line[0]] = int(line[1])

	total_counts = {}
	for gene_type in gene_types:
		total_counts[gene_type] = 0
		for gene in read_counts_by_type[gene_type]:
			total_counts[gene_type] += read_counts_by_type[gene_type][gene]

	plt.figure(1, figsize=(6,6))
	fracs = []
	for gene_type in gene_types:
		fracs.append(total_counts[gene_type])
	plt.pie(fracs, labels=gene_types, autopct='%1.1f%%', startangle=90)
	plt.axis("equal")
	plt.title('Reads per mRNA region', bbox={'facecolor':'0.8', 'pad':5}, loc = "left")
	plt.savefig(cfg.outfilepath + "Figure_Reads_By_Region_Pie.png",format='png',bbox_inches='tight',dpi=300,pad_inches=0.5)
	plt.cla()
	plt.clf()

def make_reads_by_type_pie2(nsamp, reads, threshold_nr, index_tag):

	outfilepathToRead=cfg.outfilepath + '/PlotData_ReadsPerPipeFile.txt'
	readDF = pd.read_csv(outfilepathToRead)
	df=readDF.set_index("File_name")

	num_reads = {}
	num_reads["Genome"]=df.loc["Genome Mapped (R1)","Reads_per_file"] + df.loc["Genome Mapped (R2)","Reads_per_file"]
	num_reads["Repeats"]=df.loc["Repeat Mapped (R1)","Reads_per_file"] + df.loc["Repeat Mapped (R2)","Reads_per_file"]
	num_reads["EndoVirus"]=df.loc["EndoVirus (R1)","Reads_per_file"] + df.loc["EndoVirus (R2)","Reads_per_file"]
	num_reads["tRNA"]=df.loc["tRNA Mapped (R1)","Reads_per_file"] + df.loc["tRNA Mapped (R2)","Reads_per_file"]
	num_reads["DV"]=df.loc["DV (R1)","Reads_per_file"] + df.loc["DV (R2)","Reads_per_file"]
	num_reads["ZV"]=df.loc["ZV (R1)","Reads_per_file"] + df.loc["ZV (R2)","Reads_per_file"]
	num_reads["HCV_JFH1"]=df.loc["HCV_JFH1 (R1)","Reads_per_file"] + df.loc["HCV_JFH1 (R2)","Reads_per_file"]

	labels = ['Genome', 'EndoVirus','Repeats', 'tRNA', 'DV', 'ZV', 'HCV_JFH1']
	fracs = []
	label = []
	readsum=0
	num_reads_frac={}
	for i in range(len(labels)):
		readsum=readsum+num_reads[labels[i]]
	for i in range(len(labels)):
		num_reads_frac[labels[i]]=float(num_reads[labels[i]])/float(readsum)
	
	for i in range(len(labels)):
		if num_reads_frac[labels[i]]>=0.001:
			fracs.append(num_reads[labels[i]])
			label.append(labels[i])

	patches, texts, extra = plt.pie(fracs, autopct='%1.1f%%', startangle=90)
	plt.legend(patches, label, fontsize = "x-small", loc=9, bbox_to_anchor=(0.5, -0.1))
	plt.axis("equal")
	plt.title('Reads', bbox={'facecolor':'0.8', 'pad':5}, loc = "left")


def isolateUniqueReads():
	reads_by_type = {}
	read_counts_by_type = {}
	gene_types = ["5pUTR", "3pUTR", "CDS", "Introns"]

	for gene_type in gene_types:
		reads_by_type[gene_type] = set()
		read_counts_by_type[gene_type] = {}
		name = cfg.outfilepath + "PlotData_ReadsPerGene_%s.txt" % gene_type
		if not glob.glob(name):
			continue
		with open(name, 'r') as f:
			lines = list(csv.reader(f))
			for line in lines:
				reads_by_type[gene_type].add(line[0])
				read_counts_by_type[gene_type][line[0]] = line[1]

	reads_by_type_reduced = {}
	for gene_type in gene_types:
		reads_by_type_reduced[gene_type] = reads_by_type[gene_type]
		for gene_type2 in gene_types:
			if gene_type != gene_type2:
				reads_by_type_reduced[gene_type] = reads_by_type_reduced[gene_type] - reads_by_type[gene_type2]

	for gene_type in gene_types:
		with open(cfg.outfilepath + "PlotData_ReadsPerGene_exc%s.txt" % gene_type, 'w') as f:
			for gene in reads_by_type_reduced[gene_type]:
				f.write(gene + ",")
				f.write(str(read_counts_by_type[gene_type][gene]) + "\n")

def plot_snorna_type():
	with open(cfg.outfilepath + "PlotData_readsPerSnoRNAType.txt", 'r') as f:
		lines = list(csv.reader(f, delimiter = ","))
		lines.pop(0)
		counts = [0 for i in range(3)]
		types = ["C", "H", "sca"]

		for line in lines:
			if line[0] == "C":
				counts[0] += int(line[1])
			elif line[0] == "H":
				counts[1] += int(line[1])
			else:
				counts[2] += int(line[1])

	plt.figure(1, figsize=(6,6))
	patches, texts, extra= plt.pie(counts, autopct='%1.1f%%',  startangle=90)
	plt.legend(patches, types, fontsize = "x-small", loc=9, bbox_to_anchor=(0.5, -0.1))
	plt.axis("equal")
	plt.title('Reads per snoRNA type', bbox={'facecolor':'0.8', 'pad':5}, loc = "left")
	plt.savefig(cfg.outfilepath + "Figure4b.png",format='png',bbox_inches='tight',dpi=300,pad_inches=0.5)
	plt.cla()
	plt.clf()

def plotTopEndo():
	a=np.random.random(30)
	cs=cm.Set1(np.arange(30)/30.)

	with open(cfg.outfilepath + "endoVirus_numReads.txt", 'r') as f:
		lines = list(csv.reader(f, delimiter = '\t'))
		lines.pop(0)

	names = []
	read_counts = []
	total = 0
	for i in range(25):
		names.append(lines[i][0])
		read_counts.append(float(lines[i][3]))
		total += float(lines[i][3])

	labels = []
	for i in range(25):
		labels.append(names[i] +  ", " + str(round(read_counts[i]/total * 100, 4)) + "%" + ", " + str(int(read_counts[i])))

	plt.figure(1, figsize=(6,6))
	patches, texts = plt.pie(read_counts, startangle=90, colors = cs)
	plt.legend(patches, labels, fontsize = "x-small", loc=9, bbox_to_anchor=(0.5, -0.1))
	plt.axis("equal")
	plt.title('Top 25 EndoVirus', bbox={'facecolor':'0.8', 'pad':5})
	plt.savefig(cfg.outfilepath + "Figure6.png",format='png',bbox_inches='tight',dpi=300,pad_inches=0.5)
	plt.cla()
	plt.clf()

	 
def clean_up():
	os.chdir(cfg.outfilepath)
	os.system("mkdir figures")	
	os.system("mkdir RepeatRNA")
	os.system("mv PlotData_RepeatRNA* RepeatRNA")
	os.system("mkdir ReadsPerGene")
	os.system("mv PlotData_ReadsPerGene* ReadsPerGene")
	os.system("mkdir .PlotData_Other")
	os.system("mv PlotData* .PlotData_Other")
	os.system("mv .PlotData_Other PlotData_Other")
	os.system("mkdir ProteinCoding")
	os.system("mv clipGenes_proteinCoding* ProteinCoding")
	os.system("mkdir rawdata_and_statistics")
	os.system("mv *stats* *.bam rawdata_and_statistics")
	os.system("mv rawdata_and_statistics rawdata_and_stats")
	os.system("mkdir bedfiles")
	os.system("mv *.mergedRT.bed bedfiles")
	os.system("mv *.bw *.bedgraph *_cleaned_sorted.bed *_centerCoord.bed bedfiles")
	os.system("mv *_ens_annotated.bed bedfiles")
	os.system("mkdir tRNA")
	os.system("mv *isotypes* *trnaReads.txt *histo.txt tRNA")
	os.system("mkdir todelete")
	os.system("mv *.* todelete")
	os.system("mv clipGenes_* todelete")

	
