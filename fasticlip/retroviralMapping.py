
# endoVirusMapping.py
# 1/8/14
# endoVirus_mapping.ipynb, in script form

import sys, os, re, cmath, math, glob, subprocess, csv, matplotlib, argparse
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

#matplotlib.rcParams['savefig.dpi'] = 2 * matplotlib.rcParams['savefig.dpi']
csv.register_dialect("textdialect", delimiter='\t')

# *** INPUT : Sample to process. Please fill these in
# ---------------------------------------------------

### Parsing arguments ###

parser = argparse.ArgumentParser(description="endoVirusMapping.py: a script to map raw CLIP data to endoVirus sequences", epilog="Example: python endoVirusMapping.py -n ddx3wt --GRCh38")
group = parser.add_mutually_exclusive_group()
parser.add_argument('-n', metavar='SAMPLENAME', help="Sample name; name of folder inside results/", required=True)
group.add_argument('--GRCh38', action='store_true', help="required if your CLIP is from human")
group.add_argument('--GRCm38', action='store_true', help="required if your CLIP is from mouse")
parser.add_argument('-m', metavar='MAPQ', type=int, help="Minimum MAPQ score allowed. Default is 42.", default=42)

# organism
args = parser.parse_args()
if not (args.GRCh38 or args.GRCm38):
	print "Error: must include --GRCh38 or --GRCm38. Exiting."
	exit()
if args.GRCh38: 
	org = 'GRCh38' 
	orgname='human'
	retroIndex='docs/GRCh38/retroviral/retro_spaced'
	retro_pos='docs/GRCh38/retroviral/Hs_retroviralIndex_spaced_positions.txt'
	retro_genome='docs/GRCh38/retroviral/Hs_retroviralIndex_spaced.fa'
else: 
	org = 'GRCm38'
	orgname = 'mouse'
	retroIndex='docs/GRCm38/retroviral/retro_spaced'
	retro_pos='docs/GRCm38/retroviral/Mm_retroviralIndex_spaced_positions.txt'
	retro_genome='docs/GRCm38/retroviral/Mm_retroviralIndex_spaced.fa'

mapq = str(args.m)
sampleName=args.n #name of folder
if not glob.glob("results/" + sampleName):
	print "Folder %s does not exist. Exiting." %sampleName
	exit()

# Get raw data
global outfilepath
outfilepath=os.getcwd() + '/results/%s'%sampleName

def modifyName(filepath,newTag):
	# Useage: Modifies the filepath name. 
	# Input: File path of format <path>/<name>.fastq and a string to add to the name.
	# Output: Returns the modified path of type <old path>_<new modifier>.fastq
    head, tail = os.path.split(filepath)
    oldname = tail.split('.')[0]
    newName = head+"/"+oldname+"_"+newTag
    return newName

def runBowtie(fastqFiles):
    # Useage: Short read mapping to reference (GRCh38).
    # Input: Fastq files of replicates (not paired end).
    # Output: Path to samfile for each read.
    program = 'bowtie2'
    mappedReads=[]
    unMappedReads=[]
    print "Performing Bowtie..."
    for infastq in fastqFiles:

        outfile = modifyName(infastq,"mappedToendoVirus.sam")
        
        proc = subprocess.Popen([program,'-x',retroIndex,'-U',infastq,'-S',outfile],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out, err = proc.communicate()
        result = out.decode()
        error = err.decode()
        print "Result : ",result 
        print "Error : ",error
        mappedReads = mappedReads + [outfile]
    return mappedReads

# Run Bowtie
readsForendoVirusmapping=glob.glob(outfilepath+"/*_notMappedToTrna.fastq")
mappedReads = runBowtie(readsForendoVirusmapping)

def runSamtools(samfiles, mapq):
    # Useage: Samfile processing.
    # Input: Sam files from Bowtie mapping.
    # Output: Duplicate removed, sorted bedFiles.
    program = 'samtools'
    program2 = 'bamToBed'
    outBedFiles=[]
    
    for samfile in samfiles:

        # Convert to bamfile
        bamfile = samfile.replace('.sam', '.bam')  
        proc = subprocess.Popen( [program, 'view', '-q', mapq, '-bS', '-o', bamfile, samfile])
        proc.communicate()
        
        # Sort the bamfile and note that samtools sort adds the .bam handle
        bamfile_sort = bamfile.replace('.bam', '_sorted') 
        proc2 = subprocess.Popen( [program, 'sort', bamfile, bamfile_sort])
        proc2.communicate()
        
        # Convert to bedFile
        bedFile = bamfile_sort.replace('_sorted', '_withDupes.bed')
        outfh = open(bedFile, 'w')
        proc3 = subprocess.Popen( [program2,'-i', bamfile_sort+'.bam'],stdout=outfh)
        proc3.communicate()
        
        outBedFiles=outBedFiles+[bedFile]
        
    return outBedFiles

# Run Samtools
mappedReads=glob.glob(outfilepath+"/*mappedToendoVirus.sam")
print "Process mapped data"  
mappedBedFiles=runSamtools(mappedReads, mapq)

def makeRepeatAnnotation():    
    # Repeat index sequence 
    repeatGenomeBuild=retro_genome
    repeat_genome=np.genfromtxt(repeatGenomeBuild,dtype='string')
    repeat_genome_bases=repeat_genome[1]
    repeat_genome_size=len(repeat_genome[1])
    # Repeat index positions
    repeatAnnotation=retro_pos
    repeatAnnotDF=pd.DataFrame(pd.read_table(repeatAnnotation,header=None))
    repeatAnnotDF.columns=['Name','Length','IndexStart','IndexEnd']
    # Python list extraction is not end index inclusive; to extract sequence, use end + 1.
    repeatAnnotDF['End_for_extraction']=repeatAnnotDF['IndexEnd']+1 
    repeatAnnotDF=repeatAnnotDF.set_index('Name',drop=False)
    return (repeat_genome_bases,repeatAnnotDF)
# - Repeat annotation - 
repeat_genome_bases,repeatAnnotDF=makeRepeatAnnotation()

# - Get mapped bedfiles - 
def readBed(path):
    bedFile = pd.read_table(path,dtype=str,header=None)
    bedFile.columns=['Index','Start','Stop','Name','QS','Strand']
    bedFile['Start']=bedFile['Start'].astype(int)
    return bedFile

mappedBed=glob.glob(outfilepath+"/*mappedToendoVirus_withDupes.bed")
bedR1=readBed(mappedBed[0])
bedR2=readBed(mappedBed[1])

# - Record gene position of RT stop - 
recordHits=pd.DataFrame()
for ix in repeatAnnotDF.index:
    end=repeatAnnotDF.loc[ix,'IndexEnd']
    repName=repeatAnnotDF.loc[ix,'Name']
    hits_r1=bedR1[(bedR1['Start']<int(repeatAnnotDF.loc[ix,'IndexEnd'])) & (bedR1['Start']>int(repeatAnnotDF.loc[ix,'IndexStart']))].shape[0]
    hits_r2=bedR2[(bedR2['Start']<int(repeatAnnotDF.loc[ix,'IndexEnd'])) & (bedR2['Start']>int(repeatAnnotDF.loc[ix,'IndexStart']))].shape[0]
    recordHits.loc[repName,'hits_r1']=hits_r1
    recordHits.loc[repName,'hits_r2']=hits_r2
recordHits['sum']=recordHits['hits_r1']+recordHits['hits_r2']
recordHits.fillna(0,inplace=True)
recordHits.sort(['sum'],inplace=True,ascending=False)


grZero=recordHits[recordHits['sum']>0]
pathToSave=outfilepath + '/endoVirus_numReads.txt'
recordHits.to_csv(pathToSave,sep='\t')

# - Evaluate consistency between replicates - 
def plotRepGraph():
      fig = plt.figure()
      plt.subplot(2,2,1)
      plt.scatter(recordHits['hits_r1']/bedR1.shape[0],recordHits['hits_r2']/bedR2.shape[0])
      plt.xticks(fontsize=5)
      plt.yticks(fontsize=5)
      plt.xlabel('Normalized Hits / endoVirus RNA (Rep1)',fontsize=5)
      plt.ylabel('Normalized Hits / endoVirus RNA (Rep2)',fontsize=5)

      plt.subplot(2,2,2)
      plt.scatter(np.log10(recordHits['hits_r1']),np.log10(recordHits['hits_r2']))
      plt.xticks(fontsize=5)
      plt.yticks(fontsize=5)
      plt.xlabel('Log10 (Raw hits / endoVirus RNA, Rep1)',fontsize=5)
      plt.ylabel('Log10 (Raw hits / endoVirus RNA, Rep2)',fontsize=5)
      return fig

pp = PdfPages(outfilepath + '/endoVirus_RNA_rep_comparison.pdf')
fig = plotRepGraph()
pp.savefig(fig)
pp.close()
# Binning and a matrix
ofile = open(outfilepath + '/endoVirus_RNA_bins.txt', 'w')
writer = csv.writer(ofile, 'textdialect')

#for repName in grZero.index:
for repName in recordHits.index:
	# Hits
    	hits_r1=bedR1[(bedR1['Start']<int(repeatAnnotDF.loc[repName,'IndexEnd'])) & (bedR1['Start']>int(repeatAnnotDF.loc[repName,'IndexStart']))]
	hits_r2=bedR2[(bedR2['Start']<int(repeatAnnotDF.loc[repName,'IndexEnd'])) & (bedR2['Start']>int(repeatAnnotDF.loc[repName,'IndexStart']))]

	binSize=1
	bins=range(repeatAnnotDF.loc[repName,'IndexStart'],repeatAnnotDF.loc[repName,'IndexEnd']+2,binSize) # Make sure bins are end coordinate inclusive
	# amin
	histr1,bins=np.histogram(hits_r1['Start'],bins=bins)
	histr2,bins=np.histogram(hits_r2['Start'],bins=bins)
	
	# Histogram
	numBins=100
	a=map(lambda x: float(x)*len(histr1)/numBins, range(numBins)) # convert my desired scale to the current scale		
	hits_r1_bin = np.interp(a, range(len(histr1)), histr1)
	a=map(lambda x: float(x)*len(histr2)/numBins, range(numBins)) # convert my desired scale to the current scale		
	hits_r2_bin = np.interp(a, range(len(histr2)), histr2)

	hits_bin = list(hits_r2_bin)
	outputRow = [repName]
	outputRow.extend(hits_bin)
	writer.writerow(outputRow)
	
ofile.close()
#os.system(outfilepath + '/*.sam')


