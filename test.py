import argparse, glob

parser = argparse.ArgumentParser(description="FASTCLIP: a pipeline to process CLIP data")
parser.add_argument('-i', metavar='INPUT', nargs='+', help="input FASTQ file(s); list with spaces in between", required=True)
group = parser.add_mutually_exclusive_group()
group.add_argument('--hg19', action='store_true', help="required if your CLIP is from human")
group.add_argument('--mm9', action='store_true', help="required if your CLIP is from mouse")
parser.add_argument('-n', metavar='NAME', help="Name of output directory", required=True)
parser.add_argument('-o', metavar='OUTPUT', help="Name of directory where output directory will be made", required=True)
parser.add_argument('-f', metavar='N', type=int, help="Number of bases to trim from 5' end of each read. Default is 13.", default=13)
parser.add_argument('-e', metavar='SEQUENCE', help="3' adapter to trim from the end of each read. Default is AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG.", default='AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG')

args = parser.parse_args()
if not args.hg19 or args.mm9:
	print "Error: must include --hg19 or --mm9. Exiting."
	exit()
	
infiles = args.i
for fn in infiles:
	if not glob.glob(fn):
		print "Error: input file " + fn + " not acccessible. Exiting."
		exit()
		
print args.n
outfilepath = args.o
if not glob.glob(outfilepath):
	print "Error: output directory " + outfilepath + " not acccessible. Exiting."
	exit()
	
print args.f