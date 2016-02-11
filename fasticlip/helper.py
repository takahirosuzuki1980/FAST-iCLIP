# helper.py

import os, cmath, math, sys, glob, subprocess, re, argparse, shutil, datetime, csv, commands
import numpy as np
import cfg
from matplotlib_venn import venn2
import pandas as pd
from collections import defaultdict
from operator import itemgetter
import matplotlib as mpl
import matplotlib.pyplot as plt
from optparse import OptionParser
mpl.rcParams['savefig.dpi'] = 2 * mpl.rcParams['savefig.dpi']
mpl.rcParams['path.simplify'] = True
csv.register_dialect("textdialect",delimiter='\t')

################
### TRIMMING ###
################

def log(str):
	print str
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
			uniq_reads=uniq_reads+[outread]
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

def runBowtie(processed_reads, repeat_index, retro_index, trna_index, star_index, star_ratio):
	# Usage: Read mapping.
	# Input: Fastq files of replicate trimmed read files.
	# Output: Path to samfile for each read.

	rep_sam = []
	retro_sam = []
	trna_sam = []
	genome_sam = []
	for infastq in processed_reads:
		rep_mapped = infastq.replace(".fastq", "_mappedToRepeat.sam")
		rep_mapped = cfg.outfilepath + os.path.basename(rep_mapped)
		rep_sam.append(rep_mapped)

		rep_unmapped = infastq.replace(".fastq", "_notMappedToRepeat.fastq")
		rep_unmapped = cfg.outfilepath + os.path.basename(rep_unmapped)
		retro_mapped = rep_unmapped.replace("_notMappedToRepeat.fastq", "_mappedToRetro.sam")
		retro_sam.append(retro_mapped)
		
		retro_unmapped = rep_unmapped.replace("_notMappedToRepeat.fastq", "_notMappedToRetro.fastq")
		trna_mapped = retro_unmapped.replace("_notMappedToRetro.fastq", "_mappedToTrna.sam")
		trna_sam.append(trna_mapped)

		trna_unmapped = retro_unmapped.replace("_notMappedToRetro.fastq", "_notMappedToTrna.fastq")
		genome_star_prefix = trna_unmapped.replace("_notMappedToTrna.fastq", "")
		genome_star_output = trna_unmapped.replace("_notMappedToTrna.fastq", "Aligned.out.sam")
		genome_mapped = trna_unmapped.replace("_notMappedToTrna.fastq", "_mappedToGenome.sam")
		genome_sam.append(genome_mapped)

		if glob.glob(genome_mapped):
			log("Bowtie/STAR already done.")
			continue
			
		log("Mapping {} to repeat".format(infastq))
		cmd = "bowtie2 -p 8 -x {} {} --un {} -S {} > {} 2>&1".format(repeat_index, infastq, rep_unmapped, rep_mapped, rep_mapped + '_stats.txt')
		if cfg.verbose: log(cmd)
		os.system(cmd)
		
		log("Mapping {} to retroviral".format(infastq))
		cmd = "bowtie2 -p 8 -x {} {} --un {} -S {} > {} 2>&1".format(retro_index, rep_unmapped, retro_unmapped, retro_mapped, retro_mapped + '_stats.txt')
		if cfg.verbose: log(cmd)
		os.system(cmd)
				
		log("Mapping {} to tRNA".format(infastq))
		cmd = "bowtie2 -p 8 -x {} {} --un {} -S {} > {} 2>&1".format(trna_index, retro_unmapped, trna_unmapped, trna_mapped, trna_mapped + '_stats.txt')
		if cfg.verbose: log(cmd)
		os.system(cmd)
		
		#cmd = "bowtie2 -p 8 -x {} {} -S {} > {} 2>&1".format(genome_index, trna_unmapped, genome_mapped, genome_mapped + '_stats.txt')
		log("Mapping {} to genome".format(infastq))
		cmd = "STAR --genomeDir {} --runThreadN 8 --genomeLoad LoadAndKeep --readFilesIn {} --outFileNamePrefix {} --alignEndsType EndToEnd --outFilterMismatchNoverLmax {}".format(star_index, trna_unmapped, genome_star_prefix, star_ratio)
		if cfg.verbose: log(cmd)
		star_res = os.system(cmd)
		if star_res != 0:
			log("STAR must be version 2.4.0 or higher; please check that this is the case.")
			exit()
		os.system("mv {} {}".format(genome_star_output, genome_mapped))
		
	return rep_sam, retro_sam, trna_sam, genome_sam

def run_samtools(samfiles, mapq):
	# Usage: Samfile processing (also takes unique mappers only)
	# Input: Sam files from Bowtie mapping.
	# Output: Sorted bedFiles.
	program = 'samtools'
	program2 = 'bamToBed'
	out_bedfiles = []
	for samfile in samfiles:
		bamfile_sort = samfile.replace('.sam','_sorted') 
		bedfile = bamfile_sort.replace('_sorted', '_withDupes.bed') 
		out_bedfiles.append(bedfile)
			
		cmd_1 = "cat {} | samtools view -q {} -Suo - - | samtools sort - {}".format(samfile, mapq, bamfile_sort)
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

def remove_blacklist_retro(gen_bed, blacklistregions, repeatregions):
	# Usage: Remove repeat regions from bedfile following mapping.
	# Input: .bed file after mapping (duplicates removed by samtools) and blastlist regions removed.
	# Output: Bedfile with repeat regions removed.

	masked=[]
	in_no_repeat_RTstop_collapsed=[]
	for bedIn in gen_bed:
		no_blacklist = bedIn.replace('.bed', '_noBlacklist.bed')
		no_repeat = bedIn.replace('.bed', '_noRepeat.bed')
		repeat = bedIn.replace('.bed', '_repeat.bed')
			
		cmd_1 = "bedtools intersect -a {} -b {} -wa -v -sorted > {}".format(bedIn, blacklistregions, no_blacklist)
		cmd_2 = "bedtools intersect -a {} -b {} -wa -v -sorted -s > {}".format(no_blacklist, repeatregions, no_repeat)
		cmd_3 = "bedtools intersect -a {} -b {} -wa -wb -sorted -s > {}".format(no_blacklist, repeatregions, repeat)

		if cfg.verbose: log(cmd_1)
		os.system(cmd_1)
		if cfg.verbose: log(cmd_2)
		os.system(cmd_2)
		if cfg.verbose: log(cmd_3)
		os.system(cmd_3)
		masked.append(no_repeat)
		
		# #For some extra testing & running threshold vs correlation plot
		# #Get 5' ends (3' for neg strand) and count table of collapsed RT stop positions
		# no_repeat_RTstop = bedIn.replace('.bed', '_noRepeat_fiveprime.bed')		
		# no_repeat_RTstop_collapsed = bedIn.replace('.bed', '_noRepeat_fiveprime_collapsed.txt')	
		# in_no_repeat_RTstop_collapsed.append(no_repeat_RTstop_collapsed)
		# cmd_4 =	"cat {} | awk -F\"\\t\" ' BEGIN {{ OFS=\"\\t\" }} $6==\"+\"{{ print $1, $2, $2+1, $4, $5, $6 }} $6==\"-\" {{ print $1, $3-1, $3, $4, $5, $6 }} ' | cut -f1,2,3,4,5,6 > {} ".format(no_repeat, no_repeat_RTstop)
		# cmd_5 =	"cat {} | cut -f1,2,6 | awk ' BEGIN {{ OFS=\"\\t\" }} {{ a[$0]++ }} END {{ for(i in a) {{ print i, a[i] }} }} ' > {} ".format(no_repeat_RTstop, no_repeat_RTstop_collapsed)
		
		# print cmd_4
		# os.system(cmd_4)
		# print cmd_5
		# os.system(cmd_5)
	
	# #Run correlation table with R script
	# in_no_repeat_RTstop_collapsed = ' '.join(in_no_repeat_RTstop_collapsed)
	# cmd_6 = "Rscript correlation_filter.R {]".format(in_no_repeat_RTstop_collapsed)
	# print cmd_6
	# os.system(cmd_6)
			
	return masked
	

def separateStrands(mappedReads):
	# Usage: separate positive and negative strands.
	# Input: Paths to two bed files from Samtools.
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
				neg_edit.write('\t'.join((chrom, str(int(end) - 1), str(int(end) + 29), name, quality, strand)) + '\n') # end-1 since it's [start, end)
	return negativeStrandEdit

def isolate5prime(strandedReads):
	# Usage: Isolate only the Chr, 5' position (RT stop), and strand.
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
				f.write('\t'.join((chrom,start,strand))+'\n')
	return RTstops

def fileCat(destinationFile,fileList):
	f = open(destinationFile, "w")
	for tempfile in fileList:
		readfile = open(tempfile, "r")
		f.write(readfile.read())
		readfile.close()
	f.close()

def RTcounts(RTfile):
	posRT_R1=pd.DataFrame(pd.read_table(RTfile,index_col=None,header=None,sep='\t'))
	posRT_R1.columns=['Chr','Start','Strand']
	cts=posRT_R1.groupby(['Chr','Start']).size()
	return cts

def countPassed(x, n):
	ct = 0
	for i in x:
		if i >= n: ct += 1
	return ct
	
def mergeRT(RTstopFiles, outfilename, statsfilename, minpass, threshold, expand, strand):
	# Usage: Merge RT stops between replicates and keep only those positions that exceed threshold.
	# Input: Files with RT stops for each replicate, outfile, threshold, strand, and bases to expand around RT stop.
	# Output: None. Writes merged RT stop file.
	
	cts = [RTcounts(fn) for fn in RTstopFiles]
	m = pd.concat(cts, axis=1, join='inner')
	numpass = m.apply(lambda x: countPassed(x, threshold), axis=1)
	m = m[numpass >= minpass]
	m_filter = m.copy(deep=True)
	m_filter['sum'] = m.apply(sum, axis=1)
	m_filter['mean'] = m.apply(np.mean, axis=1)
	m_filter['stdev'] = m.apply(np.std, axis=1)
	
	f = open(outfilename, 'w')
	fs = open(statsfilename, 'w')
	for i in m_filter.index:
		chrom = i[0]
		RT = i[1]
		count = m_filter.loc[i,'sum']
		mean = str(m_filter.loc[i, 'mean'])
		stdev = str(m_filter.loc[i, 'stdev'])
		if RT > expand:
			read = '\t'.join((chrom, str(int(RT)-expand), str(int(RT)+expand), 'CLIPread', '255', strand))+'\n'
			sread = '\t'.join((chrom, str(int(RT)-expand), str(int(RT)+expand), 'CLIPread', '255', strand, mean, stdev))+'\n'
		else:
			read = '\t'.join((chrom,str(int(RT)), str(int(RT)+expand), 'CLIPread', '255', strand))+'\n'
			sread = '\t'.join((chrom,str(int(RT)), str(int(RT)+expand), 'CLIPread', '255', strand, mean, stdev))+'\n'
		f.write(read*(count))
		fs.write(sread)
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

################
### CLIPPER  ###
################

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
	
def make_gene_list_from_annotation(CLIPPERlowFDR):
	outfile = CLIPPERlowFDR + '.geneNames'
	cmd = "cut -f10 {} > {}".format(CLIPPERlowFDR, outfile)
	os.system(cmd)
	return outfile
	
def runCLIPPER(RTclusterfile,genome,genomeFile):
	# Usege: Process the mergedRT file and pass through CLIPper FDR script.
	# Input: Merged RT file.
	# Output: CLIPper input (.bed) file and output file.
	
	bamfile_sorted = RTclusterfile.replace('.bed','.srt')  
	cmd1 = "bedToBam -i {} -g {} | samtools sort - {}".format(RTclusterfile, genomeFile, bamfile_sorted)
	os.system(cmd1)
	
	bamfile_sorted += ".bam"
	mapStats=bamfile_sorted.replace('.srt.bam','.mapStats.txt') 
	cmd2 = "samtools flagstat {} > {}".format(bamfile_sorted, mapStats)
	cmd3 = "samtools index {}".format(bamfile_sorted)
	os.system(cmd2)
	os.system(cmd3)
	
	CLIPPERout_dup = RTclusterfile.replace('.bed','_CLIP_clusters_dupl') 
	cmd4 = "clipper --bam {} {} --outfile={} > /dev/null 2>&1".format(bamfile_sorted, genome, CLIPPERout_dup)
	os.system(cmd4)
	
	# added by BD 4/12/15 to merge adjacent clip clusters and remove duplicates
	CLIPPERout = CLIPPERout_dup.replace('_CLIP_clusters_dupl','_CLIP_clusters') 
	with open(CLIPPERout_dup,'r') as ifile, open(CLIPPERout,'w') as ofile:
		reader = csv.reader(ifile, 'textdialect')
		writer = csv.writer(ofile, 'textdialect')
		currRow = ['chr1',0,0,0,0,'+']
		for row in reader:
			currStart = int(currRow[1])
			currEnd = int(currRow[2])
			newStart = int(row[1])
			newEnd = int(row[2])
			if currStart==newStart and currEnd==newEnd: continue #duplicates
			if math.fabs(newStart-currEnd) <= 15 and currRow[5]==row[5]: #overlap and same strand
				if int(currRow[1]) != 0: #not the first one
					currRow[2]=newEnd #merge the two adjacent clusters
			else: #not overlap
				if int(currRow[1]) != 0:
					writer.writerow(currRow)
				currRow = row #cycle continues
		writer.writerow(currRow) #fencepost
	
	return CLIPPERout

def modCLIPPERout(CLIPPERin, CLIPPERout):
	# Usage: Process the CLIPper output and isolate lowFDR reads based upon CLIPper windows.
	# Input: .bed file passed into CLIPper and the CLIPper windows file.
	# Output: Low FDR reads recovered using the CLIPer windows file, genes per cluster, gene list of CLIPper clusters, and CLIPper windows as .bed.
	
	CLIPperOutBed = CLIPPERout + '.bed' # CLIPper windows as a bed file
	CLIPperReadsPerCluster = CLIPPERout + '.readsPerCluster' # Number of reads per CLIPper cluster
	CLIPperGeneList = CLIPPERout + '.geneNames' # Gene names returned from the CLIPper file
	CLIPPERlowFDR = CLIPperOutBed.replace('.bed','_lowFDRreads.bed') # Low FDR reads returned filtered through CLIPper windows
	
	with open(CLIPPERout,'r') as infile, open(CLIPperOutBed,'w') as f, open(CLIPperReadsPerCluster,'w') as g, open(CLIPperGeneList,'w') as h:
		for line in infile:	
			try:
				# *** Old CLIPper includes a header that cannot be parsed. Handle this. ***
				# *** Old CLIPper: Ensembl genes are parsed with <name>_<cluster>_<count>. ***
				chrom,start,end,name,stats,strand,start_2,end_2 = line.strip().split('\t')
				readPerCluster=name.strip().split('_')[2]
				geneName=name.strip().split('_')[0].split('.')[0]
				f.write('\t'.join((chrom,start,end,name,stats,strand))+'\n')
				g.write((readPerCluster+'\n'))
				h.write((geneName+'\n'))
			except:
				continue
				
	# Intersect input reads with the CLIPper windows and report full result for both.
	cmd = "bedtools intersect -a {} -b {} -wa -wb -s > {}".format(CLIPPERin, CLIPperOutBed, CLIPPERlowFDR)
	os.system(cmd)
	
	return (CLIPPERlowFDR, CLIPperReadsPerCluster, CLIPperGeneList, CLIPperOutBed)

def getBedCenterPoints(inBed, expand):
	# Usage: Obtain center coordinates of bedFile.
	# Input: BedFile.
	# Output: Center coordinates returned.
	outBed = inBed.replace('.bed','_centerCoord.bed')	
	with open(inBed, 'r') as ifile, open(outBed, 'w') as ofile:
		reader = csv.reader(ifile, 'textdialect')
		writer = csv.writer(ofile, 'textdialect')
		for row in reader:	
			writer.writerow([row[0], int(row[1]) + expand, int(row[1]) + expand + 1, row[9], row[4], row[5]])
	return outBed

def cleanBedFile(inBed):
	# Usage: Sort and recover only first 6 fields from a bed file.
	# Input: BedFile.
	# Output: Sorted bedFile with correct number of fields.
	sortedBed = inBed.replace('.bed','_cleaned_sorted.bed')	
	cmd = "cut -f1-6 {} | sort -k1,1 -k2,2n > {}".format(inBed, sortedBed)
	os.system(cmd)
	return sortedBed

def makeBedGraph(cleanBed,sizesFile):
	# Usage: From a bedFile, generate a bedGraph and bigWig.
	# Input: BedFile.
	# Output: BedGraph file.
	outname = cleanBed.replace('.bed','.bedgraph')
	outname2 = cleanBed.replace('.bed','.bw')
	cmd1 = "bedtools genomecov -bg -split -i {} -g {} > {}".format(cleanBed, sizesFile, outname)
	cmd2 = cfg.home + "/bin/bedGraphToBigWig {} {} {}".format(outname, sizesFile, outname2)
	os.system(cmd1)
	os.system(cmd2)
	return outname
		
def makeClusterCenter(windowsFile):
	# Usage: Generate a file of cluster centers.
	# Input: Raw CLIPper output file.
	# Output: File with coordinates for the center of each CLIPper cluster.
	cleanBed = cleanBedFile(windowsFile)
	centers=cleanBed.replace('.bed','.clusterCenter')
	f = open(centers, 'w')
	with open(cleanBed, 'r') as infile:
		for line in infile:
			elementList = line.strip().split('\t')
			diff=abs(int((int(elementList[1])-int(elementList[2]))/2))
			f.write(elementList[0]+'\t'+str(int(elementList[1])+diff)+'\t'+str(int(elementList[1])+diff+1)+'\n')
	f.close()
	return centers

def getClusterIntensity(bedGraph,centerCoordinates):
	# Usage: Generate a matrix of read itensity values around CLIPper cluster center.
	# Input: BedGraph and cluster center file.
	# Output: Generates a matrix, which is passed into R.
	program=cfg.home + '/bin/grep_chip-seq_intensity.pl'
	cmd = "perl {} {} {} > /dev/null".format(program, centerCoordinates, bedGraph)
	os.system(cmd)
	
#####################
### PARTITIONING  ###
#####################

def getLowFDRReadTypes(CLIPPERlowFDR,pathToGeneLists):
	# Usage: Given a list of genes, return all reads for the associated genes.
	# Input: Gene list and the path to lowFDR read file.
	# Output: List of reads assocaited with the given genes.
	lowFDRgenelist = []
	for path in pathToGeneLists:
		outfile = path + '_LowFDRreads.bed'
		cmd = "grep -F -f {} {} > {}".format(path, CLIPPERlowFDR, outfile)
		os.system(cmd)
		lowFDRgenelist=lowFDRgenelist+[outfile]
	return lowFDRgenelist

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

def getLowFDRGeneTypes(CLIPperGeneList,geneAnnot):
	# Usage: Get all genes listed under each type, compare to CLIPper targets. (--clipper option)
	# Input: .bed file passed into CLIPper and the CLIPper windows file.
	# Output: Path to file containing all CLIPper genes of each type.
	geneTypes=[]
	for genepath in geneAnnot:
		if 'snoRNA' not in genepath:
			lowFDRgenes = compareLists(CLIPperGeneList, genepath, os.path.split(genepath)[1])
			geneTypes.append(lowFDRgenes)
	return geneTypes
	
def get_gene_counts(bedFile):
	bf=pd.DataFrame(pd.read_table(bedFile,header=None))
	bf.columns=['Chr','Start','Stop','CLIPper_name','Q','Strand']
	bf['geneName']=bf['CLIPper_name'].apply(lambda x: x.split('_')[0])
	geneCounts=bf.groupby('geneName').size()
	geneCounts.sort(ascending=False)
	return geneCounts

def getSnoRNAreads(negAndPosMerged,snoRNAindex):
	program='intersectBed'		
	bedFile=cfg.outfilepath+'clipGenes_snoRNA_LowFDRreads.bed'
	cmd = "bedtools intersect -a {} -b {} -s -wa -wb -sorted > {}".format(negAndPosMerged, snoRNAindex, bedFile)
	os.system(cmd)
	return bedFile

def countSnoRNAs(bedFile_sno):
	bf=pd.DataFrame(pd.read_table(bedFile_sno,header=None))
	bf.columns=['Chr','Start','End','CLIPper_name','Q','Strand','Chr_snoRNA','Start_snoRNA','Stop_snoRNA','name_snoRNA','Type','strand_snoRNA']
	geneCounts=bf.groupby('name_snoRNA').size()
	geneCounts.sort(ascending=False)
	return geneCounts

def countRemainingGeneTypes(remaining):
	for bedFile in remaining:
		try:
			bf=pd.DataFrame(pd.read_table(bedFile,header=None))
			bf.columns=['Chr','Start','End','ReadName','Q','Strand','CLIPper_winChr','CLIPper_winStart','CLIPper_winEmd','CLIPper_winaName','CLIPper_winP','CLIPper_winStrand']
			# *** THIS MAY DEPEND UPON THE VERSION OF CLIPPER USED ***
			bf['geneName']=bf['CLIPper_winaName'].apply(lambda x: x.split('_')[0])
			geneCounts=bf.groupby('geneName').size()
			geneCounts.sort(ascending=False) 
						
			head,fname=os.path.split(bedFile)
			geneType=fname.split("_")[1]
			cfg.outfilepathToSave=cfg.outfilepath+'/PlotData_ReadsPerGene_%s'%geneType
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
	fivePfile = cfg.outfilepath + '/PlotData_ReadsPerGene_5pUTR'
	threePfile = cfg.outfilepath + '/PlotData_ReadsPerGene_3pUTR'
	CDSfile = cfg.outfilepath + '/PlotData_ReadsPerGene_CDS'
	exonfile = cfg.outfilepath + '/PlotData_ReadsPerGene_Exons'
	intronfile = cfg.outfilepath + '/PlotData_ReadsPerGene_Introns'
	allfile = cfg.outfilepath + '/PlotData_ReadsPerGene_proteinCoding'
	
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
	
	ofn = cfg.outfilepath + '/PlotData_proteinCoding_byRegion'
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

def makeAvgGraph(bedGraph,utrFile,genesFile,sizesFile):
	# Usage: Generate a matrix of read intensity values across gene body.
	# Input: BedGraph.
	# Output: Generates two matrices.
	program1 = cfg.home + '/bin/bedGraph2tab.pl'
	tabFile = bedGraph.replace('.bedgraph','.tab')
	cmd = "perl {} {} {} {} {} > /dev/null".format(program1, genesFile, sizesFile, bedGraph, tabFile)
	os.system(cmd)

	program2 = cfg.home + "/bin/averageGraph_scaled_tab.pl"
	outhandle = tabFile.replace('.tab','_UTRs')
	cmd = "perl {} {} {} {} {} > /dev/null".format(program2, utrFile, tabFile, tabFile, outhandle)
	os.system(cmd)	

def getGeneStartStop(bedFile,geneRef):
	try:
		bf=pd.DataFrame(pd.read_table(bedFile,header=None))
		bf.columns=['Chr','Start','End','ReadName','Q','Strand','CLIPper_winChr','CLIPper_winStart','CLIPper_winEmd','CLIPper_winaName','CLIPper_winP','CLIPper_winStrand']
		bf['geneName']=bf['CLIPper_winaName'].apply(lambda x: x.split('_')[0].split('.')[0])
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
	filesToCount.extend([base + '_trimmed_mappedToRepeat_withDupes.bed' for base in bases])
	filesToCount.extend([base + '_trimmed_mappedToTrna_withDupes.bed' for base in bases])
	filesToCount.extend([base + '_trimmed_mappedToGenome_withDupes.bed' for base in bases])
	filesToCount.extend([base + '_trimmed_mappedToGenome_withDupes_noRepeat.bed' for base in bases])
	
	if cfg.run_clipper:
		clipperIN = cfg.outfilepath+cfg.sampleName+'_threshold=%s_%s_allreads.mergedRT_snoRNAremoved_miRNAremoved.bed'%(threshold_nr,index_tag)
		clipperOUT = cfg.outfilepath+cfg.sampleName+'_threshold=%s_%s_allreads.mergedRT_snoRNAremoved_miRNAremoved_CLIP_clusters_lowFDRreads.bed'%(threshold_nr,index_tag)
		filesToCount.extend([clipperIN, clipperOUT])
	
	fileNameBases = ['Raw', 'No dupes', 'Repeat Mapped', 'tRNA Mapped', 'Genome Mapped', 'Blacklist Masked']
	fileNames = [b + " (R{})".format(i + 1) for b in fileNameBases for i in range(nsamp)]
	if cfg.run_clipper: fileNames.extend(['ClipperIn', 'ClipperOut'])
	
	counts = []
	counter = 0
	for fileString in filesToCount:
		temp = lineCount(fileString)
		if counter < 4:
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
	if cfg.verbose: log(counts)
	readDF=pd.DataFrame()
	readDF['File_name'] = fileNames
	readDF['Reads_per_file'] = counts
	outfilepathToSave=cfg.outfilepath + '/PlotData_ReadsPerPipeFile'
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
	record.sort('genesBound',inplace=True)
	outfilepathToSave=cfg.outfilepath + '/PlotData_ReadAndGeneCountsPerGenetype'
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
	
def plot_ReadsPerCluster(threshold_nr, index_tag):
	readPerCluster=cfg.outfilepath+cfg.sampleName+'_threshold=%s_%s_allreads.mergedRT_snoRNAremoved_miRNAremoved_CLIP_clusters.readsPerCluster'%(threshold_nr,index_tag)
	clust=pd.DataFrame(pd.read_table(readPerCluster,header=None))
	clust.columns=['ReadsPerCluster']
	clust=clust['ReadsPerCluster']
	interval=10
	bins=range(min(clust)-10,max(clust)+10,interval)
	hist,bins=np.histogram(clust,bins=bins)
	width=0.7*(bins[1]-bins[0]) 
	center=(bins[:-1] + bins[1:])/2 
	plt.bar(center, hist,align='center',width=width)
	locs,pltlabels = plt.yticks(fontsize=5)
	locs,pltlabels = plt.xticks(center,center,fontsize=5)
	plt.setp(pltlabels, rotation=90, fontsize=3.5)
	plt.tick_params(axis='yticks',labelsize=5) 
	plt.xlabel('Reads per cluster (bin=%s)'%interval,fontsize=5)
	plt.ylabel('Frequency (RT stop count)',fontsize=5)
	plt.title('Reads per cluster',fontsize=5)
	plt.xlim(0,100) # Make the histogram easy to view.
	# plt.xlim(-interval,np.max(center)+interval)
	
def plot_ClusterSizes(threshold_nr, index_tag):
	clipClusters=cfg.outfilepath+cfg.sampleName+'_threshold=%s_%s_allreads.mergedRT_snoRNAremoved_miRNAremoved_CLIP_clusters'%(threshold_nr,index_tag)
	clust=pd.DataFrame(pd.read_table(clipClusters,header=None,skiprows=1))
	clust.columns=['chr','start','end','name','score','strand','m1','m2']
	clust['clusterSize']=clust['start']-clust['end']
	clust['clusterSize']=clust['clusterSize'].apply(lambda x: math.fabs(x))
	plt.boxplot(clust['clusterSize'])
	plt.tick_params(axis='x',labelbottom='off') 
	ax=plt.gca()
	for line in ax.get_xticklines():
		line.set_markersize(0)
	plt.ylabel('Cluster length (bases)',fontsize=5)
	locs,pltlabels = plt.yticks(fontsize=5)
	plt.title('Cluster size',fontsize=5)

def plot_clusterBindingIntensity(threshold_nr, index_tag):
	clusterCenterHeatmap=cfg.outfilepath+cfg.sampleName+'_threshold=%s_%s_allreads.mergedRT_snoRNAremoved_miRNAremoved_CLIP_clusters_cleaned_sorted.clusterCenter_heatmap.txt'%(threshold_nr,index_tag)
	hmap=pd.DataFrame(pd.read_table(clusterCenterHeatmap,header=None,skiprows=1))
	hmap_vals=hmap.ix[:,1:]
	sums=hmap_vals.sum(axis=1)
	hmap_vals=hmap_vals.loc[np.argsort(sums),:]
	plt.ylim(0,hmap_vals.shape[0])
	p=plt.pcolormesh(np.array(hmap_vals),cmap='Blues')
	plt.tick_params(axis='x',labelbottom='off') 
	plt.xlabel('Cluster position',fontsize=5)
	locs,pltlabels = plt.yticks(fontsize=5)
	plt.ylabel('Cluster number',fontsize=5)
	plt.title('Read distribution',fontsize=5)

def readUTRfile(path):
	geneCounts=pd.read_csv(path,header=None)
	geneCounts.columns=['Gene_name','Count']
	return geneCounts

def plot_readsBymRNAregion(ax): 
	fivePfile = cfg.outfilepath+'/PlotData_ReadsPerGene_5pUTR'
	threePfile = cfg.outfilepath+'/PlotData_ReadsPerGene_3pUTR'
	CDSfile = cfg.outfilepath+'/PlotData_ReadsPerGene_CDS'
	pc_5pReads = readUTRfile(fivePfile)['Count'].sum() if glob.glob(fivePfile) else 0
	pc_3pReads = readUTRfile(threePfile)['Count'].sum() if glob.glob(threePfile) else 0
	pc_CDSReads = readUTRfile(CDSfile)['Count'].sum() if glob.glob(CDSfile) else 0
	non_intronic=pc_5pReads+pc_3pReads+pc_CDSReads
	allProteinCoding=cfg.outfilepath +'clipGenes_proteinCoding_LowFDRreads_centerCoord.bed'
	all_pc=pd.DataFrame(pd.read_table(allProteinCoding,header=None))
	pc_allReads=all_pc.shape[0]
	v=[float(pc_allReads-non_intronic)/pc_allReads,float(pc_5pReads)/pc_allReads,float(pc_CDSReads)/pc_allReads,float(pc_3pReads)/pc_allReads]
	pie_wedges=ax.pie(v,labels=["Intronic","5p UTR","CDS","3pUTR"],labeldistance=1.1,autopct='%1.1f%%')
	plt.rcParams['font.size']=5
	for wedge in pie_wedges[0]:
		wedge.set_edgecolor('black')
		wedge.set_lw(1)
		
def plot_figure1(nsamp, reads, threshold_nr, index_tag):
	if cfg.run_clipper:
		plt.subplot(2,3,1) 
		plot_ReadAccounting(nsamp, reads, threshold_nr, index_tag)
		plt.subplot(2,3,2)
		plot_ReadsPerCluster(threshold_nr, index_tag)
		plt.subplot(2,3,3)
		plot_ClusterSizes(threshold_nr, index_tag)
		plt.subplot(2,3,4)
		plot_clusterBindingIntensity(threshold_nr, index_tag)
		ax = plt.subplot(2,3,5)
		plot_readsBymRNAregion(ax)
		plt.subplot(2,3,6)
		plot_BoundGeneTypes()
	else:
		plt.subplot(2,3,1) 
		plot_ReadAccounting(nsamp, reads, threshold_nr, index_tag)
		ax = plt.subplot(2,3,2)
		plot_readsBymRNAregion(ax)
		plt.subplot(2,3,3)
		plot_BoundGeneTypes()
	plt.tight_layout()
	
############
## PLOT 2 ##
############


def plot_mRNAgeneBodyDist():
	averageGraph=cfg.outfilepath+'clipGenes_proteinCoding_LowFDRreads_centerCoord_UTRs_scaled_cds200_abt0_averageGraph.txt'
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
		screen=readUTRfile(cfg.outfilepath+'/PlotData_ReadsPerGene_5pUTR')
	elif utr=='3p':
		ix=(hmap_m[range(1,401)].sum(axis=1)==0)&(hmap_m[range(401,601)].sum(axis=1)>0)
		screen=readUTRfile(cfg.outfilepath+'/PlotData_ReadsPerGene_3pUTR')
	elif utr=='5p3p':
		ix=(hmap_m[range(201,401)].sum(axis=1)==0)&(hmap_m[range(401,601)].sum(axis=1)>0)&(hmap_m[range(1,201)].sum(axis=1)>0)
		screen=readUTRfile(cfg.outfilepath+'/PlotData_ReadsPerGene_5p3pUTR')
	else:  # utr=='CDS'
		ix=(hmap_m[range(1,201)].sum(axis=1)==0)&(hmap_m[range(401,601)].sum(axis=1)==0)&(hmap_m[range(201,401)].sum(axis=1)>0)
		screen=readUTRfile(cfg.outfilepath+'/PlotData_ReadsPerGene_CDS')
		
	# Ensure all genes are also identified in pre-allocated gene lists.
	screen['Gene_name'] = screen['Gene_name'].map(lambda x: x.split('.')[0]) # remove transcript number from gene name
	hmap_m_utrSpec=hmap_m.ix[ix,:]
	hmap_m_utrSpec_filter=pd.merge(hmap_m_utrSpec,screen,left_on='ENSG_ID',right_on='Gene_name',how='inner')
	sums=hmap_m_utrSpec_filter[range(1,601)].sum(axis=1)
	hmap_m_utrSpec_filter=hmap_m_utrSpec_filter.loc[np.argsort(sums),:]
	return hmap_m_utrSpec_filter

def plot_geneBodyPartition(ensemblGeneAnnot):
	treatMatrix=cfg.outfilepath+'clipGenes_proteinCoding_LowFDRreads_centerCoord_UTRs_scaled_cds200_abt0_treatmatrix.txt'
	hmap=pd.DataFrame(pd.read_table(treatMatrix,header=None,skiprows=1))
	
	# Ensure genes recovered from this analysis are independently identified using partitioning of CLIPper cluster data.
	hmap['ENSG_ID']=hmap.ix[:,0].apply(convertENBLids, args=(ensemblGeneAnnot,))
	bound_pc = cfg.outfilepath+'clipGenes_proteinCoding'
	pc_genes=pd.DataFrame(pd.read_table(bound_pc,header=None,))
	pc_genes.columns=['ENSG_ID']
	hmap_m=pd.merge(hmap,pc_genes,left_on='ENSG_ID',right_on='ENSG_ID',how='inner') 
	
	# Isolate intronic bound genes.
	tosave=cfg.outfilepath+'PlotData_ExclusiveBound_Intronic' 
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
		fn = cfg.outfilepath+'/PlotData_ReadsPerGene_'+name
		if not glob.glob(fn): continue
		
		plt.subplot2grid((2,4),(1,i),colspan=1)
		utrMatrix=getUTRbindingProfile(geneTypes[i],hmap_m)
		tosave=cfg.outfilepath+'PlotData_ExclusiveBound_%s'%geneTypes[i] 
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
		start = int(row[9])
		end = int(row[8])
	return np.array(RTpositions_plus), np.array(RTpositions_minus), start, end
	
def get_histogram(RTpositions, bins, repeat_genome_bases, start, end):
	hist, bins = np.histogram(RTpositions, bins=bins)
	center = (bins[:-1] + bins[1:])/2
	
	# Normalize
	histPlot = np.array(hist, dtype=float)
	histPlot = np.array(histPlot/float(len(RTpositions)), dtype=float)
		
	# Make sure same length
	sequence = repeat_genome_bases[start:end+1]
	hist = hist[:len(sequence)]
	histPlot = hist[:len(sequence)]
	center = center[:len(sequence)]
	
	return hist, histPlot, sequence, center	

def plot_repeatRNA(repeatGenomeBuild):
	repeat_genome = np.genfromtxt(repeatGenomeBuild,dtype='string')
	repeat_genome_bases = repeat_genome[1]
	
	repFiles=glob.glob(cfg.outfilepath + '/PlotData_RepeatRNAreads_*')
	
	plotDim=math.ceil(math.sqrt(len(repFiles)))
	i=0
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
		plt.bar(center_plus, histPlot_plus, align='center', width=width, color='blue', alpha=0.45)
		plt.tick_params(axis='x', labelsize=2.5) 
		plt.tick_params(axis='y', labelsize=2.5)  
		plt.title('RT stops for %s: %s'%(name, len(RTpositions_plus)), fontsize=5)
		plt.xlim(start, end)  
		
		# Record data
		storageDF = pd.DataFrame()
		storageDF['Sequence'] = pd.Series(list(sequence_plus))
		storageDF['RT_stops'] = np.array(list(hist_plus))
		storageDF['RT_stops_norm'] = np.array(list(histPlot_plus))	
		storageDF['RT_stops_minus'] = np.array(list(hist_minus))
		storageDF['RT_stops_minus_norm'] = np.array(list(histPlot_minus))	  
		outfilepathToSave = cfg.outfilepath + '/PlotData_RepeatRNAHist_%s'%name
		storageDF.to_csv(outfilepathToSave)
		i += 1
	plt.tight_layout()

############
## PLOT 4 ##
############

def plot_rDNA(start18s, end18s, start5s, end5s, start28s, end28s, rRNAend):
	plt.subplot2grid((3,3),(0,0),colspan=3)
	name='rDNA'
	rDNA=glob.glob(cfg.outfilepath + 'PlotData_RepeatRNAreads_rDNA')
	
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
	bf_sno.columns=['Chr','Start','End','CLIPper_name','Q','Strand','Chr_snoRNA','Start_snoRNA','Stop_snoRNA','name_snoRNA','Type','strand_snoRNA']
	snoTypes=pd.DataFrame(bf_sno.groupby('Type').size())
	snoTypes.columns=['Reads']
	snoTypes['Fraction']=snoTypes['Reads']/sum(snoTypes['Reads'])
	outfilepathToSave=cfg.outfilepath+'/PlotData_readsPerSnoRNAType'
	snoTypes.to_csv(outfilepathToSave)

	snoTypesAndGenes=pd.DataFrame(bf_sno.groupby(['Type','name_snoRNA']).size())
	snoTypesAndGenes.columns=['Count_per_gene']
	outfilepathToSave=cfg.outfilepath+'/PlotData_geneStatsPerSnoRNAType'
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
		
		outfilepathToSave=cfg.outfilepath+'/PlotData_snoRNAReadDist_%s'%sType
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
		name=st_file.split('clipGenes_')[1].split('_LowFDRreads')[0]
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

def clean_up():
	os.chdir(cfg.outfilepath)
	os.system("mkdir figures")
	os.system("mv Figure* figures")
	
	os.system("mkdir RepeatRNA")
	os.system("mv PlotData_RepeatRNA* RepeatRNA")
	
	os.system("mkdir ReadsPerGene")
	os.system("mv PlotData_ReadsPerGene* ReadsPerGene")
	
	os.system("mkdir PlotData_Other")
	os.system("mv PlotData* PlotData_Other")
	
	os.system("mkdir ProteinCoding")
	os.system("mv clipGenes_proteinCoding* ProteinCoding")
	
	os.system("mkdir rawdata_and_stats")
	os.system("mv *stats* *.bam *Log.final.out rawdata_and_stats")
	os.system("mv runLog rawdata_and_stats")
	
	os.system("mkdir bedfiles")
	os.system("mv *.mergedRT.bed bedfiles")
	os.system("mv *.bw *.bedGraph *_cleaned_sorted.bed *_centerCoord.bed bedfiles")
	if cfg.run_clipper:
		os.system("mv *_allreads.mergedRT_CLIP_clusters_lowFDRreads* rawdata_and_stats")
		os.system("mv *.mergedRT_CLIP_clusters.bed rawdata_and_stats")
	else:
		os.system("mv *_ens_annotated.bed bedfiles")

	os.system("mkdir tRNA")
	os.system("mv *iclipro/ tRNA")
	os.system("mv *isotypes* *trnaReads.txt *histo.txt tRNA")

	os.system("mkdir todelete")
	os.system("mv *.* todelete")
	os.system("mv clipGenes_* todelete")
	os.system("mv rawdata_and_stats/runLog ./")
	