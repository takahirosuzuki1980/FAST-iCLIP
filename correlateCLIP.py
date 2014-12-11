# correlateCLIP.py
# 11/23/14
# correlates RT stops from two CLIP files (*noBlacklist_noRepeat.bed)
# uses RT stop maker from FAST-iCLIP

import csv, sys, collections, glob, os
csv.register_dialect("textdialect", delimiter='\t')

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
			if glob.glob(neg_strand) and glob.glob(pos_strand): continue
				
			neg = open(neg_strand, 'w')
			pos = open(pos_strand, 'w')
			for line in infile:	
				if str(line.strip().split('\t')[5]) == '-':
					neg.write(line)
				elif str(line.strip().split('\t')[5]) == '+':
					pos.write(line)
	return [negativeStrand,positiveStrand]

def modifyNegativeStrand(negativeStrandReads):
	# Usage: For negative stranded reads, ensure 5' position (RT stop) is listed first.
	# Input: Bed file paths to all negative stranded.
	# Output: Paths to modified bed files.
	negativeStrandEdit=[]
	for negativeRead in negativeStrandReads:
		neg_strand_edited=negativeRead.replace('_neg.bed','_negEdit.bed')
		negativeStrandEdit=negativeStrandEdit+[neg_strand_edited]
		if glob.glob(neg_strand_edited): continue
		
		neg_edit = open(neg_strand_edited, 'w')
		with open(negativeRead, 'r') as infile:
			for line in infile:	
				chrom,start,end,name,quality,strand=line.strip().split('\t')
				neg_edit.write('\t'.join((chrom,end,str(int(end)+30),name,quality,strand))+'\n')
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
			if glob.glob(RTstop): continue
			
			f = open(RTstop, 'w')
			for line in infile:	
				chrom,start,end,name,quality,strand=line.strip().split('\t')
				f.write('\t'.join((chrom,start,strand))+'\n')
	return RTstops

def getHashmap(filePair):
	fn1 = filePair[0]
	fn2 = filePair[1]
	regionToCount = collections.defaultdict(lambda: [0,0])
	
	# read in file 1
	ifile = open(fn1, 'r')
	reader = csv.reader(ifile, 'textdialect')
	for row in reader:
		key = '__'.join(row)
		regionToCount[key][0] += 1
	ifile.close()
	
	# read in file 2
	ifile = open(fn2, 'r')
	reader = csv.reader(ifile, 'textdialect')
	for row in reader:
		key = '__'.join(row)
		regionToCount[key][1] += 1
	ifile.close()
	
	return regionToCount
	
def main():
	fn1, fn2, ofn = sys.argv[1], sys.argv[2], sys.argv[3]
	if len(sys.argv) < 4: 
		print "fn1, fn2, ofn required. Exiting."
		exit()
		
	print "RT stop isolation."
	mappedBedFiles_rep = [fn1, fn2]
	readsByStrand_rep=separateStrands(mappedBedFiles_rep)
	modified = modifyNegativeStrand(readsByStrand_rep[0])
	negativeRTstop_rep=isolate5prime(modified) 
	positiveRTstop_rep=isolate5prime(readsByStrand_rep[1])

	print "Making hashmap of RT stop to count"
	regionToCount = getHashmap(positiveRTstop_rep)
	regionToCount_neg = getHashmap(negativeRTstop_rep)
	regionToCount.update(regionToCount_neg)
	
	for fn in modified + negativeRTstop_rep + positiveRTstop_rep + readsByStrand_rep[0] + readsByStrand_rep[1]:
		os.system("rm -f " + fn)	
	
	print "Writing list of RT stops"
	regionKeys = regionToCount.keys()
	list1 = [regionToCount[key][0] for key in regionKeys]
	list2 = [regionToCount[key][1] for key in regionKeys]
	
	ofile = open(ofn, 'w')
	writer = csv.writer(ofile, 'textdialect')
	for region in regionToCount:
		'''
		if regionToCount[region][0] * regionToCount[region][1] == 0: continue #both need to have RT stops
		'''
		if regionToCount[region][0] + regionToCount[region][1] < 1: continue #this is the threshold that was set earlier
		regions = region.split('__')
		regions.extend(regionToCount[region])
		writer.writerow(regions)
	ofile.close()
		
	print "Running R script"
	print "Rscript ~/CLIP/correlateCLIP_plot.r " + ofn
	os.system("Rscript ~/CLIP/correlateCLIP_plot.r " + ofn)
			
if __name__ == '__main__':
	main()
