# motifAnalysis.py
# 12/30/14
# Motif_analysis.ipynb, in script form

import sys, os, re, cmath, math, glob, subprocess, csv, matplotlib, argparse
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
plt.style.use('ggplot')

matplotlib.rcParams['savefig.dpi'] = 2 * matplotlib.rcParams['savefig.dpi']
csv.register_dialect("textdialect", delimiter='\t')

### Parsing arguments ###

parser = argparse.ArgumentParser(description="motifAnalysis.py: use this to look for motifs in CLIP data", epilog="Example: python motifAnalysis.py -n ddx3wt --GRCh38 -l 12,14 -f 15 -m 25")
group = parser.add_mutually_exclusive_group()
parser.add_argument('-n', metavar='SAMPLENAME', help="Sample name; name of folder inside results/", required=True)
parser.add_argument('-p', metavar='PREFIX', help="Prefix for folders where motifs will be stored", default="homer")
group.add_argument('--GRCh38', action='store_true', help="required if your CLIP is from human")
group.add_argument('--GRCm38', action='store_true', help="required if your CLIP is from mouse")
parser.add_argument('-l', metavar='LENGTHS', help="Length(s) of motifs desired; separated by a ','. Default is '6,8,10,12'.", default="6,8,10,12")
parser.add_argument('-f', metavar='FILTER_WINDOW', help="Distance over which RT stops cannot overlap. Default is 10.", default=10)
parser.add_argument('-m', metavar='MOTIF_WINDOW', help="Size of window to find motifs around RT stop. Default is 20.", default=20)
parser.add_argument('--window', action='store_true', help="look at genes with reads in a certain window")
parser.add_argument('-s', metavar='START', type=int, help="Left end of window. 0 = left, 200 = 5UTR/CDS boundary, 400 = CDS/3UTR boundary, 600 = right. Default is 0.", default=0)
parser.add_argument('-e', metavar='END', type=int, help="Right end of window. Default is 600.", default=600)


# organism
args = parser.parse_args()
if not (args.GRCh38 or args.GRCm38):
	print "Error: must include --GRCh38 or --GRCm38. Exiting."
	exit()
if args.GRCh38: 
	org = 'GRCh38' 
	orgname='human'
else: 
	org = 'GRCm38'
	orgname = 'mouse'

sampleName=args.n #name of folder
if not glob.glob("results/" + sampleName):
	print "Folder %s does not exist. Exiting." %sampleName
	exit()
	
lengths=args.l #desired motif lengths, separated by commas with NO spaces
filter_win=int(args.f) #window around RT stop over which no other RT stops can be found
homer_win=int(args.m) #window around RT stop to find motifs within. if homer_win = 20 then we do pos +- 10
# Start/end of region to extract, between 0 and 300
if args.window:
	start=args.s
	end=args.e



# - Functions for extracting cluster data associated with gene list - 
# -------------------------------------------------------------------
def grep(pattern,filein):
	# Usage: Open a file and search all lines for a pattern.
	# Input: Pattern to search (gene name) and file name.
	# Output: List of lines in a file that have the pattern.
	r = []
	filein_open = open(filein, 'r')
	for line in filein_open:
		 # print line
		 if re.search(pattern,line):
			r.append(line)
	filein_open.close()
	return r

def extractClusters(geneList,allRT, filter_win):
	# Usage: Extract a set of  RT stops based on gene name
	# Input: Gene list
	# Output: RT stops in those genes
	# BD Update 12/13/2014: only unique RT stops
	# BD Update 6/11/2015: filterSummit.pl AND no more unique RT stops
	extractedRT = geneList+'_rtstops.bed'
	outfh = open(extractedRT, 'w')
	writer = csv.writer(outfh, 'textdialect')
	
	# List of allowed gene names (genes that have RT stops only in the defined region)
	namesToQuery = np.genfromtxt(geneList,usecols=(0,),delimiter='\t',dtype='string')
	
	# Iterate through all RT stops in clusters
	ifile = open(allRT, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	for row in reader:
		geneName = row[3].split('_')[0].split('.')[0]
		if geneName in namesToQuery: writer.writerow(row)
		
	ifile.close()
	outfh.close()
	
	# Run filtering
	proc = subprocess.Popen(["dos2unix",extractedRT])
	proc.communicate()
	filteredRT = geneList+'_rtstops_summits.bed'	
	proc = subprocess.Popen(["perl","filterSummit.pl",extractedRT, filteredRT, str(filter_win)])
	proc.communicate()	
	
	return filteredRT

# - Functions for running HOMER analysis - 
def shuffleBedFile(inBed):
	# Usage: Shuffle a bed file (usually used a background file for HOMER analysis) 
	# Input: Bedfile
	# Output: SHuffled bedFile
	program = 'shuffleBed'
	referenceFile = os.getcwd()+'/docs/%s/%s_transcriptome_collapse_exon.bed'%(org,org)
	genomeFile = os.getcwd()+'/docs/%s/%s.%s.genome'%(org, orgname, org)
	try:
		shuffledBed = inBed.replace('.bed','_shuffled.bed')
		outfh = open(shuffledBed, 'w')
		proc = subprocess.Popen([program,'-i',inBed,'-incl',referenceFile,'-g',genomeFile],stdout=outfh)
		proc.communicate()
		return shuffledBed
	except:
		print "Problem generating shuffled bedfile."
		
def makeBedForHOMER(inBed, homer_win):
	# Usage: This modified the bedfile for processing homer by making the first field a concatenation of chr_start_end
	# Input: Clean bedfile with first 5 field properly assigned
	# Output: Modified bed file with first field chr_start_end, and name excluded.
	# BD 6/11/2015: adds a +/-win to each RT stop identified as the summit of each cluster
	# Make sure bedfile only has 5 fields
	print inBed
	bedForHOMER=inBed.replace('.bed','_forHOMER.bed')	
	print bedForHOMER
	# Open new file
	f = open(bedForHOMER, 'w')
	with open(inBed, 'r') as infile:
		for line in infile:	
			elementList = line.strip().split('\t')
			if len(elementList) != 6: continue
			# Re-write the bed file with chr replaced
			f.write('\t'.join((elementList[0]+'_'+elementList[1]+'_'+elementList[2],elementList[0],str(int(elementList[1])-homer_win/2),str(int(elementList[2])+homer_win/2),elementList[5],'\n')))
	f.close()
	return bedForHOMER

def runHOMER(inBed,lengths,homer_win,outDirName):
	# Usage: Run the HOMER motif finding algorithm 
	# Input: Bedfile properly modified for HOMER
	# Output: A directory containing the HOMER output files
	program='findMotifsGenome.pl'
	homerReferenceFile = os.getcwd()+'/docs/%s/%s_transcriptome_collapse_exon.bed'%(org, org)
	# Convert the input bedFile into HOMER compatible format
	inBedForHOMER=makeBedForHOMER(inBed, homer_win)
	# Get the path of the input file 
	path,filename=os.path.split(inBedForHOMER)
	outDir=path+'/'+outDirName
	
	# print path, filename
	# return

	# Call HOMER, which will generate a directory of files
	proc = subprocess.Popen([program,inBedForHOMER,org,outDir,'-rna','-bg',homerReferenceFile, '-len', lengths, '-S', '10'])
	proc.communicate()
	

def convertENBLids(inNames):
	# Usage: Convert ENST to ENSG (unique ID) using ENSEMBL annotation file
	# Input: List of ENST IDs
	# Output: List of ENSG IDs
	genesFile = os.getcwd() + '/docs/%s/%s_ensembl_genes.txt'%(org, org)
	ensemblIDfile=np.genfromtxt(genesFile,usecols=(1,12,),delimiter='\t',dtype='string') # Note that column lookup is zero indexed
	temp=[]
	for name in inNames:
		outName=ensemblIDfile[ensemblIDfile[:,0]==name,1]
		temp=temp+[outName]
	temp=np.array(temp)
	return temp

def customGeneExtraction(outfilepath,startIndex,endIndex):
	# Extract genes with reads that fall within a specified window 
	# Input: File path, and coordinates for data selection
	# Output: None
	
	# Extract data from treat matrix
	treatMatrixCols=600
	treatMatrix=glob.glob(outfilepath+'/rawdata/*treatmatrix.txt')[0]
	treatMatrixData=np.genfromtxt(treatMatrix,skip_header=1,usecols=range(1,treatMatrixCols+1),delimiter='\t',dtype='float')
	geneNames=np.loadtxt(treatMatrix,dtype='string',skiprows=1,usecols=(0,),delimiter='\t')
	#print geneNames[0:10]
	# pathToNameConversion=os.getcwd() + '/docs/refSeq_to_Ensl_all.txt'
	# nameConversionToEnsembl=np.genfromtxt(pathToNameConversion,usecols=(0,1,),delimiter='\t',dtype='string') # Note that column lookup is zero indexed
	
	# Convert to ENSG IDs and check for genes in the initial list 
	geneNames=convertENBLids(geneNames)
	masterList = outfilepath+'/rawdata/clipGenes_proteinCoding'
	masterNames = np.genfromtxt(masterList,usecols=(0,),delimiter='\t',dtype='string') # Gene names isolated from Figure 1e
	indexer=[]
	for geneName in geneNames:
		if geneName in masterNames:
			indexer=indexer+[1]
		else:
			indexer=indexer+[0]
	indexer=np.array(indexer,dtype=bool)
	geneNames=geneNames[indexer]
	treatMatrixData=treatMatrixData[indexer,:]
	
	customData=treatMatrixData[treatMatrixData[:,range(startIndex,endIndex+1)].sum(axis=1) > 0,:]
	customNames=geneNames[treatMatrixData[:,range(startIndex,endIndex+1)].sum(axis=1) > 0]
	
	tosave=outfilepath+'TEMP'
	np.savetxt(tosave,customData,fmt="%s")
	
	treatSums=customData.sum(axis=1)
	sortedIndex=list(reversed([i[0] for i in sorted(enumerate(treatSums),key=lambda x:x[1])]))
	sortedData=customData[sortedIndex,:]
	sortedNames=customNames[sortedIndex]
	tosave=outfilepath+'SourceData_CustomGeneExtraction'
	np.savetxt(tosave,np.unique(sortedNames),fmt="%s")
	return tosave

def plotExtractedRegion(sample,plotNum,start,end):
	outfilepath=os.getcwd() + '/results/%s/'%sample
	# Bed file with protein coding reads
	filteredProteinCoding = outfilepath+'/rawdata/clipGenes_proteinCoding_reads_centerCoord_snoRNAremoved_miRNAremoved.bed'
	averageGraph=outfilepath+'/rawdata/clipGenes_proteinCoding_reads_centerCoord_snoRNAremoved_miRNAremoved_cleaned_sorted_UTRs_scaled_cds200_abt0_averageGraph.txt'	
	# Number of columns 
	avgGraphCols=600
	avgGraphData=np.loadtxt(averageGraph,skiprows=1,dtype='float',usecols=range(1,avgGraphCols+1))
	plt.subplot(1,2,plotNum)
	ylimit=max(avgGraphData[1,:])*1.1
	a= plt.plot(avgGraphData[1,:],color='blue',linewidth='1')
	plt.ylim(0,ylimit)
	plt.vlines(200,0,ylimit,linestyles='dashed')
	plt.vlines(400,0,175,linestyles='dashed')
	plt.tick_params(axis='x',labelbottom='off') 
	plt.axvspan(start,end,facecolor='r',alpha=0.5)
	plt.tick_params(axis='y',labelsize=5) 
	plt.title(sample)
	
	b=PdfPages("extractedRegion.pdf")
	b.savefig(a)
	b.close()
	
	
	
	
	
	
	
# Running things
# --------------
outfilepath=os.getcwd()+'/results/%s/'%sampleName

Clusters = glob.glob(outfilepath + '/bedfiles/*threshold*mergedRT_snoRNAremoved_miRNAremoved_CLIP_clusters_reads_cleaned_sorted.bed')
if not Clusters:
	Clusters = glob.glob(outfilepath + '/bedfiles/*threshold*mergedRT_snoRNAremoved_miRNAremoved_ens_annotated.bed')[0]
else:
	Clusters = Clusters[0]
	
print "Cluster file to process: {}".format(Clusters)

if not args.window:
	# regions we want to look at
	UTRcluster_names = ['5p', '3p', 'cds', 'introns', 'exons']
	UTRclusters = ["clipGenes_proteinCoding_reads_centerCoord_" + x + '.bed' for x in UTRcluster_names]

	# Run HOMER on all clusters
	filteredRT = Clusters[:-4] + "_rtstops_summits.bed"
	proc = subprocess.Popen(["perl","bin/filterSummit.pl", Clusters, filteredRT, str(filter_win)])
	proc.communicate()
	runHOMER(filteredRT, lengths, homer_win, args.p + '_allReads')
	
	# shuffledReads=shuffleBedFile(filteredRT)
	# runHOMER(shuffledReads,lengths,'homer_allReads_shuffle')

	# Run HOMER on clusters for specific regions
	folderNames = [args.p + '_' + x for x in UTRcluster_names]
	i=0
	for clusterFile in UTRclusters:
		infile = outfilepath + '/ProteinCoding/' + clusterFile
		print "{} file to process: {}".format(UTRcluster_names[i], infile)

		filteredRT = infile[:-4] + "_rtstops_summits.bed"
		proc = subprocess.Popen(["perl","bin/filterSummit.pl", infile, filteredRT, str(filter_win)])
		proc.communicate()
		
		runHOMER(filteredRT, lengths, homer_win, folderNames[i])
		i+=1
else:
	# Input: Extract all reads in the region of interest
	pathToCustomFile=customGeneExtraction(outfilepath,start,end)

	# Plot the region
	plotExtractedRegion(sampleName,1,start,end)

	# Run HOMER
	outfilepath=os.getcwd()+'/results/%s/'%sampleName
	clusters=extractClusters(pathToCustomFile,Clusters)
	runHOMER(clusters,lengths,'rawdata/homer_custom')