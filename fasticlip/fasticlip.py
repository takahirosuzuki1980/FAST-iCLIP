#!/usr/bin/env python
# FAST-iCLIP

import os, cmath, math, sys, glob, subprocess, re, argparse, shutil, datetime, csv
from helper import *
import cfg
import numpy as np
import pandas as pd
from collections import defaultdict
from operator import itemgetter
import matplotlib as mpl
import matplotlib.pyplot as plt
from optparse import OptionParser
mpl.rcParams['savefig.dpi'] = 2 * mpl.rcParams['savefig.dpi']
mpl.rcParams['path.simplify'] = True
csv.register_dialect("textdialect",delimiter='\t')

### Checking environment variables ###
cfg.home = os.environ.get("FASTICLIP_PATH")
if not cfg.home and glob.glob("docs/"):
	cfg.home = os.getcwd()
if not cfg.home:
	print "Error: Environment variable FASTICLIP_PATH not found. Please set it to the FAST-iCLIP installation directory."
	exit()
if not glob.glob(cfg.home + "/docs/"):
	print "Error: docs folder not inside FASTICLIP_PATH. Make sure this environment variable is set correctly, or \
	that you have run ./configure inside the installation directory."
	exit()

### Parsing arguments ###
parser = argparse.ArgumentParser(description="FAST-iCLIP: a pipeline to process iCLIP data", epilog="Example: fasticlip -i rawdata/example_MMhur_R1.fastq rawdata/example_MMhur_R2.fastq --mm9 -n MMhur -o results")
parser.add_argument('-i', metavar='INPUT', nargs='+', help="Up to 4 input fastq or fastq.gz files separated by a space", required=True)
parser.add_argument('--trimmed', action='store_true', help="flag if files are already trimmed")
group = parser.add_mutually_exclusive_group()
group.add_argument('--GRCh38', action='store_true', help="required if your CLIP is from human")
group.add_argument('--GRCm38', action='store_true', help="required if your CLIP is from mouse")
parser.add_argument('-s', metavar="STAR_INDEX", help="Path to STAR index for your organism", required=True)
parser.add_argument('-n', metavar='NAME', help="Name of output directory", required=True)
parser.add_argument('-o', metavar='OUTPUT', help="Name of directory where output directory will be made", required=True)
parser.add_argument('-f', metavar='N', type=int, help="First base to keep on 5' end of each read. Default is 18.", default=18)
parser.add_argument('-a', metavar='ADAPTER', help="3' adapter to trim from the end of each read. Default is AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG.", default='AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG')
parser.add_argument('-tr', metavar='REPEAT_THRESHOLD_RULE', type=str, help="m,n: at least m samples must each have at least n RT stops mapped to repeat RNAs. Default is 1,4 (1 sample); 2,3 (2 samples); x,2 (x>2 samples)")
parser.add_argument('-tv', metavar='VIRAL_THRESHOLD_RULE', type=str, help="m,n: at least m samples must each have at least n RT stops mapped to viral RNAs. Default is 1,4 (1 sample); 2,3 (2 samples); x,2 (x>2 samples)")
parser.add_argument('-tn', metavar='NONREPEAT_THRESHOLD_RULE', type=str, help="m,n: at least m samples must each have at least n RT stops mapped to nonrepeat RNAs. Default is 1,4 (1 sample); 2,3 (2 samples); x,2 (x>2 samples)")
parser.add_argument('-sr', metavar='STAR_RATIO', type=float, help="Maximum mismatches per base allowed for STAR genome mapping (corresponds to outFilterMismatchNoverLmax). Default is 0.08 (2 mismatches per 25 mapped bases).", default=0.08)
parser.add_argument('-bm', metavar='BOWTIE_MAPQ', type=int, help="Minimum MAPQ (Bowtie alignment to repeat/tRNA/retroviral indexes) score allowed. Default is 42.", default=42)
parser.add_argument('-q', metavar='Q', type=int, help="Minimum average quality score allowed during read filtering. Default is 25.", default=25)
parser.add_argument('-p', metavar='P', type=int, help="Percentage of bases that must have quality > q during read filtering. Default is 80.", default=80)
parser.add_argument('-l', metavar='L', type=int, help="Minimum length of read. Default is 15.", default=15)
parser.add_argument('--verbose', action='store_true', help="Print everything")

# organism
args = parser.parse_args()
if not (args.GRCh38 or args.GRCm38):
	print "Error: must include --GRCh38 or --GRCm38. Exiting."
	exit()
if args.GRCh38: org = 'human' 
else: org = 'mouse'

# input files
reads = args.i
reads_files = []
for fn in reads:
	reads_files.append(fn)
	if not glob.glob(fn):
		print "Error: input file " + fn + " not accessible. Exiting."
		exit()
		
# sample name and output directory and STAR index
cfg.sampleName = args.n
cfg.outfilepath = args.o
if not glob.glob(cfg.outfilepath):
	print "Error: output directory " + cfg.outfilepath + " not accessible. Exiting."
	exit()
cfg.outfilepath = cfg.outfilepath + '/%s/'%cfg.sampleName
if not glob.glob(cfg.outfilepath): os.system("mkdir " + cfg.outfilepath)

star_index = args.s  # Location of STAR index
if not glob.glob(star_index):
	print "Error: STAR index at " + star_index + " not accessible. Exiting."
	exit()
	
# exoViruses
exoViruses = []
exoViruses = ['DV', 'ZV', 'HCV_JFH1']
for v in exoViruses:
	viral_pref = cfg.home + "/docs/viral/{}".format(v)
	if not glob.glob(viral_pref + ".*"):
		print "Error: Viral index with prefix {} not accessible. Exiting.".format(viral_pref)
		exit()

# Create log and start pipeline
cfg.logFile=cfg.outfilepath + "runLog"
cfg.logOpen=open(cfg.logFile, 'w')

### Parameters ###

iCLIP3pBarcode = args.a # Barcode sequence to trim from reads.
star_ratio = args.sr  # Maximum mismatch/base ratio allowed for STAR
mapq = args.bm # Minimum MAPQ score allowed allowed for Bowtie2
q = args.q # Minimum quality score to keep during filtering.
p = args.p # Percentage of bases that must have quality > q during filtering.
l = args.l +  args.f - 1 # Minimum length of read + 5' adapter
iCLIP5pBasesToTrim=args.f # Number of reads to trim from 5' end of clip reads + 1 (this number represents the first base kept)
expand = 15 # Bases to expand around RT position after RT stops are merged.
cfg.verbose = True if args.verbose else False
NameDelim = '_' # Delimiter that for splitting gene name in the CLIPper windows file.

nsamp = len(reads)
if not args.tn:
	if nsamp == 1: threshold_nr = 4
	elif nsamp == 2: threshold_nr = 3
	else: threshold_nr = 2
	minpass_nr = nsamp
else:
	try: [minpass_nr, threshold_nr] = [int(x) for x in args.tn.split(',')]
	except: 
		print "Nonrepeat threshold rule must be in form m,n where m,n are integers"
		exit()
	minpass_nr = min(minpass_nr, nsamp) 

if not args.tr:
	if nsamp == 1: threshold_rep = 4
	elif nsamp == 2: threshold_rep = 3
	else: threshold_rep = 2
	minpass_rep = nsamp
else:
	try: [minpass_rep, threshold_rep] = [int(x) for x in args.tr.split(',')]
	except: 
		print "Repeat threshold rule must be in form m,n where m,n are integers."
		exit()
	minpass_rep = min(minpass_rep, nsamp) 

if not args.tv:
	if nsamp == 1: threshold_viral = 4
	elif nsamp == 2: threshold_viral = 3
	else: threshold_viral = 2
	minpass_viral = nsamp
else:
	try: [minpass_viral, threshold_viral] = [int(x) for x in args.tv.split(',')]
	except: 
		print "Viral threshold rule must be in form m,n where m,n are integers."
		exit()
	minpass_viral = min(minpass_viral, nsamp) 

if org == 'human':
	repeat_index=cfg.home + '/docs/GRCh38/repeat/rep_spaced' # bt2 index for repeat RNA.
	repeatGenomeBuild=cfg.home+'/docs/GRCh38/repeat/repeatRNA_spaced.fa' # Sequence of repeat index.
	repeatAnnotation=cfg.home+'/docs/GRCh38/repeat/Hs_repeatIndex_spaced_positions.txt' # Repeat annotation file.
	endoVirus_index=cfg.home+'/docs/GRCh38/retroviral/'
	exoVirus_index=cfg.home+'/docs/viral/'
	start18s=3655
	end18s=5523
	start5s=6601
	end5s=6757
	start28s=7925
	end28s=12994
	rRNAend=13357
	index=cfg.home + '/docs/GRCh38/GRCh38_STAR' # bt2 index for mapping.
	index_tag='GRCh38' # Name of bt2 index.
	genomeFile=cfg.home+'/docs/GRCh38/GRCh38.sizes' # Genome file for bedGraph, etc.
	#blacklistregions=cfg.home+'/docs/GRCh38/wgEncodeDukeMapabilityRegionsExcludable.bed' # Blacklist masker.
	repeatregions=cfg.home+'/docs/GRCh38/GRCh38_repeatMasker.bed' # Repeat masker.
	geneAnnot=glob.glob(cfg.home+'/docs/GRCh38/gene_types/*') # List of genes by type.
	miRNAmasker=cfg.home+'/docs/GRCh38/miR_sort_clean.bed' # miRNA masker file.
	fivePUTRBed=cfg.home+'/docs/GRCh38/5pUTRs_Ensbl_sort_clean_uniq.bed' # UTR annotation file.
	threePUTRBed=cfg.home+'/docs/GRCh38/3pUTRs_Ensbl_sort_clean_uniq.bed' # UTR annotation file.
	exonBed=cfg.home+'/docs/GRCh38/Exons_Ensbl_sort_clean_uniq.bed' # UTR annotation file.
	#utrFile=cfg.home+'/docs/GRCh38/UTR_annotation.txt' # UTR annotation file.
	#genesFile=cfg.home+'/docs/GRCh38/GRCh38_ensGene_76.txt' # Gene annotation file.
	sizesFile=cfg.home+'/docs/GRCh38/GRCh38.sizes' # Genome sizes file. 
	snoRNAindex=cfg.home+'/docs/GRCh38/snoRNA_coordinates.bed' # snoRNA coordinate file.
	tRNAindex=cfg.home+'/docs/GRCh38/trna/tRNA_hg19'  # this is now CCA tailed
	geneStartStopRepo=cfg.home+'/docs/GRCh38/all_genes.txt'
	geneStartStopRepoBed = cfg.home+'/docs/GRCh38/genes_BED6.bed'
elif org == 'mouse':
	repeat_index=cfg.home + '/docs/GRCm38/repeat/rep_spaced' # bt2 index for repeat RNA.
	repeatGenomeBuild=cfg.home+'/docs/GRCm38/repeat/Mm_repeatRNA_spaced.fa' # Sequence of repeat index.
	repeatAnnotation=cfg.home+'/docs/GRCm38/repeat/Mm_repeatIndex_spaced_positions.txt' # Repeat annotation file.
	endoVirus_index=cfg.home+'/docs/GRCm38/retroviral/'
	exoVirus_index=cfg.home+'/docs/viral/'
	start18s=4007
	end18s=5876
	start5s=6877
	end5s=7033
	start28s=8123
	end28s=12836
	rRNAend=13401
	index=cfg.home + '/docs/GRCm38/GRCm38_STAR' # bt2 index for mapping.
	index_tag='GRCm38' # Name of bt2 index.
	genomeFile=cfg.home+'/docs/GRCm38/GRCm38.sizes' # Genome file for bedGraph, etc.
	#blacklistregions=cfg.home+'/docs/GRCm38/mm9-blacklist.bed' # Blacklist masker.
	repeatregions=cfg.home+'/docs/GRCm38/GRCm38_repeatMasker.bed' # Repeat masker.
	geneAnnot=glob.glob(cfg.home+'/docs/GRCm38/gene_types/*') # List of genes by type. 
	miRNAmasker=cfg.home+'/docs/GRCm38/miRNA_coordinates.bed' # miRNA masker file.
	fivePUTRBed=cfg.home+'/docs/GRCm38/5UTR.bed' # UTR annotation file.
	threePUTRBed=cfg.home+'/docs/GRCm38/3UTR.bed' # UTR annotation file. 
	exonBed=cfg.home+'/docs/GRCm38/exons.bed' # UTR annotation file. 
	#utrFile=cfg.home+'/docs/GRCm38/UTR_annotation.txt' # UTR annotation file. 
	#genesFile=cfg.home+'/docs/GRCm38/GRCm38_ensGene.txt' # Gene annotation file. 
	sizesFile=cfg.home+'/docs/GRCm38/GRCm38.sizes' # Genome sizes file. 
	snoRNAindex=cfg.home+'/docs/GRCm38/snoRNA_coordinates.bed' # snoRNA coordinate file. 
	tRNAindex=cfg.home+'/docs/GRCm38/tRNA/tRNA_mm9'  # this is now CCA tailed
	geneStartStopRepo=cfg.home+'/docs/GRCm38/all_genes.txt'
	geneStartStopRepoBed = cfg.home+'/docs/GRCm38/genes_BED6.bed'

### start running pipeline ###
now=datetime.datetime.now()
cfg.logOpen.write("Timestamp: %s\n"%str(now))
cfg.logOpen.write("\n###Parameters used###\n")
cfg.logOpen.write("3' barcode: %s\n'"%iCLIP3pBarcode)
cfg.logOpen.write("Minimum quality score (q): %s\n"%q)
cfg.logOpen.write("Percentage of bases with > q: %s\n"%p)
cfg.logOpen.write("5' bases to trim: %s\n'"%iCLIP5pBasesToTrim)
cfg.logOpen.write("Threshold for minimum number of RT stops (repeat): %s samples with >= %s RT stops\n"%(minpass_rep, threshold_rep))
cfg.logOpen.write("Threshold for minimum number of RT stops (nonrepeat): %s samples with >= %s RT stops\n"%(minpass_nr, threshold_nr))
cfg.logOpen.write("\n\n\n")

def main():

	# 1. Trim and map
	
	log("\nProcessing THis sample {}".format(cfg.sampleName))

	if not args.trimmed: 
		log("\nRemoving duplicates")
		dup_removed_reads = remove_dup(reads, q, p)

	if not args.trimmed: 
		log("\nTrimming 5' and 3'")
		processed_reads = trim(dup_removed_reads, iCLIP3pBarcode, l, iCLIP5pBasesToTrim)
	else: processed_reads = reads

	log("\nRun mapping to indexes.")
	(viral_sam, rep_sam, endoVirus_sam, trna_sam, gen_sam) = run_mapping(processed_reads, exoViruses, repeat_index, endoVirus_index, tRNAindex, star_index, star_ratio)

	log("\nRun samtools.")
	viral_bed = run_samtools(viral_sam, "-q {}".format(mapq))
	rep_bed = run_samtools(rep_sam, "-q {}".format(mapq))
	endoVirus_bed = run_samtools(endoVirus_sam, "-q {}".format(mapq))  # do more with this later, so no flags yet
	trna_bed = run_samtools(trna_sam, "-q {}".format(mapq))
	gen_bed = run_samtools(gen_sam, "-q 255") # we're using STAR here so 255 signifies unique mapping

	# 2.1 Process viral RT stops
	if exoViruses:
		log("\nViral RT stop isolation.")
		exoVirus_to_beds = defaultdict(lambda: [])
		exoVirus_to_fa = defaultdict(lambda: [])
		for v in exoViruses:
			exoVirus_to_fa[v].append(exoVirus_index+v+".fasta")
			log("virus: " + v)
			for bed in viral_bed:
				log("bed: " + bed)
				if 'mappedTo{}'.format(v) in bed:
					log("ok") 
					exoVirus_to_beds[v].append(bed)
		
		for v in exoVirus_to_beds:
			log("check: " + v  )
			viral_bedfiles = exoVirus_to_beds[v]
			readsByStrand_v = separateStrands(viral_bedfiles)
			negativeRTstop_v = isolate5prime(modifyNegativeStrand(readsByStrand_v[0])) 
			positiveRTstop_v = isolate5prime(readsByStrand_v[1]) 

			posMerged = cfg.outfilepath + cfg.sampleName + '_viral_{}_positivereads.mergedRT'.format(v)
			negMerged = cfg.outfilepath + cfg.sampleName + '_viral_{}_negativereads.mergedRT'.format(v)
			negAndPosMerged = cfg.outfilepath + cfg.sampleName + '_threshold={}_viral_{}_allreads.mergedRT.bed'.format(threshold_viral, v)
		
			mergeRT(positiveRTstop_v, posMerged, posMerged + '_stats', minpass_viral, threshold_viral, '+')
			mergeRT(negativeRTstop_v, negMerged, negMerged + '_stats', minpass_viral, threshold_viral, '-')
			fileCat(negAndPosMerged, [posMerged, negMerged])
			fileCat(negAndPosMerged + '_stats', [posMerged + '_stats', negMerged + '_stats'])

		"""log("Making RT stop histograms")
		filename = cfg.sampleName + '_threshold={}_viral_{}_allreads.mergedRT.bed'.format(threshold_viral, v)
		viral_RT_stops(negAndPosMerged,filename)"""

	# 2.2 Process repeat RT stops	
	log("\nRepeat RT stop isolation.")
	readsByStrand_rep=separateStrands(rep_bed)
	negativeRTstop_rep=isolate5prime(modifyNegativeStrand(readsByStrand_rep[0])) 
	positiveRTstop_rep=isolate5prime(readsByStrand_rep[1]) 

	log("Merge Repeat RT stops.")
	posMerged = cfg.outfilepath+cfg.sampleName+'_repeat_positivereads.mergedRT'
	negMerged = cfg.outfilepath+cfg.sampleName+'_repeat_negativereads.mergedRT'
	negAndPosMerged = cfg.outfilepath+cfg.sampleName+'_threshold=%s'%threshold_rep+'_repeat_allreads.mergedRT.bed'
	
	mergeRT(positiveRTstop_rep, posMerged, posMerged + '_stats', minpass_rep, threshold_rep, '+')
	mergeRT(negativeRTstop_rep, negMerged, negMerged + '_stats', minpass_rep, threshold_rep, '-')
	fileCat(negAndPosMerged, [posMerged, negMerged])
	fileCat(negAndPosMerged + '_stats', [posMerged + '_stats', negMerged + '_stats'])

	# 2.3 Process nonrepeat RT stops
	log("Nonrepeat RT stop isolation.")
	gen_norepeat_bed = remove_RepeatMaskerRegions(gen_bed, repeatregions)
	#gen_norepeat_bed = remove_RepeatMaskerRegions(gen_bed, blacklistregions, repeatregions)
	readsByStrand = separateStrands(gen_norepeat_bed)
	negativeRTstop = isolate5prime(modifyNegativeStrand(readsByStrand[0])) 
	positiveRTstop = isolate5prime(readsByStrand[1]) 

	log("Merge Nonrepeat RT stops.")
	posMerged = cfg.outfilepath+cfg.sampleName+'_%s_positivereads.mergedRT'%index_tag
	negMerged = cfg.outfilepath+cfg.sampleName+'_%s_negativereads.mergedRT'%index_tag
	negAndPosMerged = cfg.outfilepath+cfg.sampleName+'_threshold=%s'%threshold_nr+'_%s_allreads.mergedRT.bed'%index_tag
	mergeRT(positiveRTstop, posMerged, posMerged + '_stats', minpass_nr, threshold_nr, '+')
	mergeRT(negativeRTstop, negMerged, negMerged + '_stats', minpass_nr, threshold_nr, '-')
	fileCat(negAndPosMerged,[posMerged,negMerged])
	fileCat(negAndPosMerged + '_stats', [posMerged + '_stats', negMerged + '_stats'])

	# 3. Process genic RT stops
	log("\nGetting list of snoRNAs")
	bedFile_sno = getSnoRNAreads(negAndPosMerged, snoRNAindex)
	if os.stat(bedFile_sno).st_size > 0:
		geneCounts_sno = countSnoRNAs(bedFile_sno) 
		outfilepathToSave = cfg.outfilepath + '/PlotData_ReadsPerGene_snoRNA.txt'
		geneCounts_sno.to_csv(outfilepathToSave)
		
	log("\nFiltering out snoRNAs and miRNAs")
	sno_mirna_filtered_reads = filter_snoRNAs(negAndPosMerged, snoRNAindex, miRNAmasker)

	log("\nAnnotating reads by gene")
	GeneReads = annotate_genes(sno_mirna_filtered_reads, geneStartStopRepoBed)
	GeneList = make_gene_list_from_annotation(GeneReads) # list of genes; commented out all things making this

	log("Make bedGraphs")
	# pre-masking, gene annotation bedgraph
	bedGraphout = makeBedGraph(cleanBedFile(negAndPosMerged),genomeFile)
	allGeneCentersBedGraph = makeBedGraph(cleanBedFile(negAndPosMerged), genomeFile)

	# - all
	log("\nPartition reads by type.")
	pathToGeneLists = getGeneTypes(GeneList,geneAnnot)
	pathToReadLists = getReadTypes(GeneReads,pathToGeneLists)
		
	# - protein coding
	geneRef=pd.DataFrame(pd.read_table(geneStartStopRepo))
	proteinCodingReads = cfg.outfilepath + 'clipGenes_proteinCoding_reads.bed'
	proteinCodingReads_centered = getBedCenterPoints(proteinCodingReads, namecol=9)
	geneCounts_pc = get_gene_counts(proteinCodingReads_centered) 
	cfg.outfilepathToSave = cfg.outfilepath + '/PlotData_ReadsPerGene_proteinCoding.txt'
	geneCounts_pc.to_csv(cfg.outfilepathToSave)
	
	# - lncRNA
	lincRNAReads = cfg.outfilepath + 'clipGenes_lincRNA_reads.bed'
	lincRNAReads_centered = getBedCenterPoints(lincRNAReads, namecol=9)
	if os.stat(lincRNAReads_centered).st_size > 0:
		geneCounts_linc = get_gene_counts(lincRNAReads_centered)
		outfilepathToSave = cfg.outfilepath + '/PlotData_ReadsPerGene_lincRNA.txt'
		geneCounts_linc.to_csv(outfilepathToSave)

		bf=pd.DataFrame(pd.read_table(lincRNAReads_centered,header=None))
		bf.columns=['Chr','Start','Stop','name','Q','Strand']
		bf['geneName']=bf['name'].apply(lambda x: x.split('_')[0])
		merge=pd.merge(geneRef,bf,left_on='Ensembl Gene ID',right_on='geneName')
		ncRNA_startStop=merge[['Ensembl Gene ID','Gene Start (bp)','Gene End (bp)','Start','Stop','Strand']]
		outfilepathToSave = lincRNAReads_centered.replace(".bed",".geneStartStop")
		ncRNA_startStop.to_csv(outfilepathToSave)

	# - other
	remaining = [f for f in glob.glob(cfg.outfilepath+"*_reads.bed") if 'lincRNA' not in f and 'proteinCoding' not in f and 'snoRNA' not in f]
	countRemainingGeneTypes(remaining)


	# 5. Process introns and UTRs
	log("\nIntron and UTR analysis.")
	exonreads, intronreads, fivePreads, threePreads, cdsreads = extract_regions(proteinCodingReads_centered, fivePUTRBed, threePUTRBed, exonBed)

	if os.stat(fivePreads).st_size > 0:
		geneCounts_5p=get_gene_counts(fivePreads) 
		outfilepathToSave=cfg.outfilepath+'/PlotData_ReadsPerGene_5pUTR.txt'
		geneCounts_5p.to_csv(outfilepathToSave)

	if os.stat(threePreads).st_size > 0:
		geneCounts_3p=get_gene_counts(threePreads) 
		outfilepathToSave=cfg.outfilepath+'/PlotData_ReadsPerGene_3pUTR.txt'
		geneCounts_3p.to_csv(outfilepathToSave)

	if os.stat(cdsreads).st_size > 0:
		geneCounts_cds=get_gene_counts(cdsreads) 
		outfilepathToSave=cfg.outfilepath+'/PlotData_ReadsPerGene_CDS.txt'
		geneCounts_cds.to_csv(outfilepathToSave) 

	if os.stat(exonreads).st_size > 0:
		geneCounts_5p=get_gene_counts(exonreads) 
		outfilepathToSave=cfg.outfilepath+'/PlotData_ReadsPerGene_Exons.txt'
		geneCounts_5p.to_csv(outfilepathToSave)

	if os.stat(intronreads).st_size > 0:
		geneCounts_3p=get_gene_counts(intronreads) 
		outfilepathToSave=cfg.outfilepath+'/PlotData_ReadsPerGene_Introns.txt'
		geneCounts_3p.to_csv(outfilepathToSave)

	gene_binding_by_region()
	isolateUniqueReads()

	# 6. Analysis of gene bodies, CLIP binding sites (iCLIPro), tRNAs
	###log("\nncRNA gene body analysis.")
	###remaining=[f for f in glob.glob(cfg.outfilepath+"*_reads.bed") if 'lincRNA' not in f and 'proteinCoding' not in f and 'snoRNA' not in f]
	###for bedFile in remaining:
	###	st_stop=getGeneStartStop(bedFile,geneRef)
        ###
	###log("\nRun iCLIPro.")
	###iclipro = run_iclipro(gen_sam)

	log("\nRun tRNA isotype counting.")
	trna_readlist = trna_isotype_count(trna_sam, minpass_rep, threshold_rep)

	# 7. Repeat RNAs
	log("\nRecord repeat RNA.")
	repeat_genome_bases,repeatAnnotDF=makeRepeatAnnotation(repeatGenomeBuild,repeatAnnotation)
	repeatAnnotDF.set_index('Name',inplace=True,drop=False)
	repeatMerged=glob.glob(cfg.outfilepath+"*repeat_allreads.mergedRT.bed") # Get merged data for repeat index.
	rep=pd.read_table(repeatMerged[0],dtype=str,header=None)
	rep.columns=['Rep_index','Start','Stop','Read_name','Q','Strand']
	rep['RT_stop']=rep['Start'].astype(int)+expand
	for ix in repeatAnnotDF.index:
		end=repeatAnnotDF.loc[ix,'IndexEnd']
		repName=repeatAnnotDF.loc[ix,'Name']
		gene_hits=rep[(rep['RT_stop']<int(repeatAnnotDF.loc[ix,'IndexEnd']))&(rep['RT_stop']>int(repeatAnnotDF.loc[ix,'IndexStart']))]
		gene_hits['Repeat_End']=repeatAnnotDF.loc[ix,'IndexEnd']
		gene_hits['Repeat_Start']=repeatAnnotDF.loc[ix,'IndexStart']
		outfilepathToSave=cfg.outfilepath + '/PlotData_RepeatRNAreads_%s.txt'%repName
		gene_hits.to_csv(outfilepathToSave)

	# 8. Plots
	log("Running Retroviral Mapping")
	organism = "GRCh38" if org == "human" else "GRCm38"
	cmd = "python fasticlip/retroviralMapping.py -n %s --%s" % (cfg.sampleName, organism)
	os.system(cmd)

	log("\nMake plots.")
	import matplotlib
	import commands
	matplotlib.rcParams['savefig.dpi'] = 2 * matplotlib.rcParams['savefig.dpi']

	log("Making Figure 1")
	fig1 = plt.figure(1)
	plot_figure1(nsamp, reads, threshold_nr, index_tag, exoViruses, reads_files)
	fig1.tight_layout()
	fig1.savefig(cfg.outfilepath+'Figure1.png',format='png',bbox_inches='tight',dpi=300,pad_inches=0.5)

	log("Making Figure 2")
	fig2 = plt.figure(2)
	plot_repeatRNA(repeatGenomeBuild)
	fig2.tight_layout()
	fig2.savefig(cfg.outfilepath+'Figure2a.png',format='png',bbox_inches='tight',dpi=300,pad_inches=0.5)
       	plt.cla()
        plt.clf()

	log("Making Figure 3")
	fig3 = plt.figure(3)
	plot_rDNA(start18s, end18s, start5s, end5s, start28s, end28s, rRNAend)
	fig3.tight_layout()
	fig3.savefig(cfg.outfilepath+'Figure3a.png',format='png',bbox_inches='tight',dpi=300,pad_inches=0.5)
       	plt.cla()
        plt.clf()

	log("Making Figure 4")
	if os.stat(bedFile_sno).st_size > 0:
		fig4 = plt.figure(4)
		plot_snorna(bedFile_sno)
		fig4.tight_layout()
		fig4.savefig(cfg.outfilepath+'Figure4a.png',format='png',bbox_inches='tight',dpi=300,pad_inches=0.5, fontsize = 5)
		plot_snorna_type()
	else:
		log("No snoRNA reads; not making Figure 4")
       	plt.cla()
        plt.clf()

	log("Making Figure 5")
	st_stopFiles = glob.glob(cfg.outfilepath+"*.geneStartStop")
	st_stopFiles = [f for f in st_stopFiles if 'rRNA' not in f]
	fig5 = plt.figure(5)
	plot_ncrnas(st_stopFiles, expand)
	fig5.tight_layout()
	fig5.savefig(cfg.outfilepath+'Figure5.png',format='png',bbox_inches='tight',dpi=300,pad_inches=0.5)

	log("Making Figure 6")
	plotTopEndo()

        log("Making Figure 7")
        if exoViruses:
	        for v in exoVirus_to_beds:
			negAndPosMerged = cfg.outfilepath + cfg.sampleName + '_threshold={}_viral_{}_allreads.mergedRT.bed'.format(threshold_viral, v)
                        filename = cfg.sampleName + '_threshold={}_viral_{}_allreads.mergedRT.bed'.format(threshold_viral, v)
			viral_RT_stops(negAndPosMerged,filename, exoVirus_to_fa[v])
        clean_up()
        cfg.logOpen.close()

if __name__ == "__init__":
	main()

	
