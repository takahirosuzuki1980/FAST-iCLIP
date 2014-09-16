# processExons.py
# 8/18/14
# cleans up exons (gets all unique exons for a gene) for mm9

import csv, collections
csv.register_dialect("textdialect", delimiter='\t')

def main():
	ifile = open("/home/lmartin/CLIP/docs/mm9/processing/mm9_ensembl_genes.txt", 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	ofile = open("/home/lmartin/CLIP/docs/mm9/mm9_annotation.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	reader.next()
	for row in reader:
		geneName = row[1]
		strand = row[3]
		
		[ts, te, cs, ce] = [int(x) for x in row[4:8]]
		length = te - ts + 1
		writer.writerow([geneName, 'Transcript', 1, length])
		if strand == '+':
			writer.writerow([geneName, '5UTR', 1, cs - ts])
			writer.writerow([geneName, 'CDS', cs - ts + 1, ce - ts])
			writer.writerow([geneName, '3UTR', ce - ts + 1, length])
		else:
			writer.writerow([geneName, '5UTR', 1, te - ce])
			writer.writerow([geneName, 'CDS', te - ce + 1, te - cs])
			writer.writerow([geneName, '3UTR', te - cs + 1, length])		
	
	ifile.close()
			
	ofile.close()
	
if __name__ == '__main__':
	main()
