# processUTR.py
# 8/18/14
# cleans up UTRs (get one 5' and one 3' UTR for each gene)

import csv, collections
csv.register_dialect("textdialect", delimiter='\t')

def main():
	ifile = open("/home/lmartin/CLIP/docs/mm9/mm9_utrs_raw.txt", 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	geneTo5 = {}
	geneTo3 = {}
	reader.next()
	for row in reader:
		geneName = row[0]
		chrom = 'chr' + row[7]
		strand = '-' if row[2] == '-1' else '+'
		name = '__'.join([geneName, chrom, strand])

		if row[3] != '': 
			length = int(row[4]) - int(row[3])
			if name not in geneTo5 or geneTo5[name][1] - geneTo5[name][0] < length:
				geneTo5[name] = (int(row[3]), int(row[4]))
		if row[5] != '': 
			length = int(row[6]) - int(row[5])
			if name not in geneTo3 or geneTo3[name][1] - geneTo3[name][0] < length:
				geneTo3[name] = (int(row[5]), int(row[6]))
	ifile.close()
	
	ofile = open("/home/lmartin/CLIP/docs/mm9/processing/mm9_5pUTR.bed", 'w')
	writer = csv.writer(ofile, 'textdialect')
	for name in geneTo5:
		[geneName, chrom, strand] = name.split('__')
		outputRow = [chrom, geneTo5[name][0], geneTo5[name][1], geneName, 0, strand]
		writer.writerow(outputRow)			
	ofile.close()
	
	ofile = open("/home/lmartin/CLIP/docs/mm9/processing/mm9_3pUTR.bed", 'w')
	writer = csv.writer(ofile, 'textdialect')
	for name in geneTo3:
		[geneName, chrom, strand] = name.split('__')
		outputRow = [chrom, geneTo3[name][0], geneTo3[name][1], geneName, 0, strand]
		writer.writerow(outputRow)			
	ofile.close()
	
if __name__ == '__main__':
	main()
