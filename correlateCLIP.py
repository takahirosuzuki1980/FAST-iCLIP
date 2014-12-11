# correlateCLIP.py
# 11/23/14
# correlates RT stops from two CLIP files (*_allreads.mergedRT.bed)

import csv, sys, collections
csv.register_dialect("textdialect", delimiter='\t')

def main():
	fn1, fn2, ofn = sys.argv[1], sys.argv[2], sys.argv[3]
	if len(sys.argv) < 4: 
		print "fn1, fn2, ofn required. Exiting."
		exit()
		
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
	
	# get lists
	regionKeys = regionToCount.keys()
	list1 = [regionToCount[key][0] for key in regionKeys]
	list2 = [regionToCount[key][1] for key in regionKeys]
	
	# write file
	ofile = open(ofn, 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	for region in regionToCount:
		regions = region.split('__')
		regions.extend(regionToCount[region])
		writer.writerow(regions)
		
	# run R script
	os.system("Rscript ~/CLIP/correlateCLIP_plot.r " + ofn)
	
		
if __name__ == '__main__':
	main()