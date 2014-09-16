# processExons.py
# 8/18/14
# cleans up exons (gets all unique exons for a gene) for mm9

import csv, collections
csv.register_dialect("textdialect", delimiter='\t')

def main():
	ifile = open("/home/lmartin/CLIP/docs/mm9/processing/mm9_exons_raw.txt", 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	geneToExons = collections.defaultdict(lambda: set())
	reader.next()
	for row in reader:
		geneName = row[0]
		chrom = 'chr' + row[3]
		strand = '-' if row[2] == '-1' else '+'
		name = '__'.join([geneName, chrom, strand])
		geneToExons[name].add((row[4], row[5]))
	
	ifile.close()
	
	ofile = open("/home/lmartin/CLIP/docs/mm9/processing/mm9_exons_formatted.bed", 'w')
	writer = csv.writer(ofile, 'textdialect')
	for name in geneToExons:
		[geneName, chrom, strand] = name.split('__')
		for exon in geneToExons[name]:
			outputRow = [chrom, exon[0], exon[1], geneName, 0, strand]
			writer.writerow(outputRow)
			
	ofile.close()
	
if __name__ == '__main__':
	main()
