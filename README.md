FAST-iCLIP
==========

Fully Automated and Standardized iCLIP (FAST-iCLIP) is a fully automated tool to process iCLIP data. Please cite the following paper:

`Zarnegar B, Flynn RA, Shen Y, Do BT, Chang HY, Khavari PA. Ultraefficient irCLIP pipeline for characterization for protein-RNA interactions. Nature Methods (2016)`

This package contains two main sets of tools: an executable called `fasticlip` to run iCLIP on human and mouse data, and several (possibly deprecated) iPython notebooks to process iCLIP data from viral genomes.

The following README will focus mainly on `fasticlip`. The pdf in the repository contains further instructions for using the iPython notebooks.

**Table of Contents**

- [Usage](#user-content-usage)
	- [Required arguments](#user-content-required-arguments)
	- [Optional arguments](#user-content-optional-arguments)
- [Installation instructions](#user-content-installation-instructions)
- [Dependencies](#user-content-dependencies)
- [Input](#user-content-input)
- [Output](#user-content-output)
- [How the pipeline works](#user-content-how-the-pipeline-works)

Usage [*** UPDATE ***]
-----

`fasticlip [-h] -i INPUT [INPUT ...] [--trimmed] [--GRCh38 | --GRCm38] -s STAR_INDEX -n NAME -o OUTPUT [-f N] [-a ADAPTER] [-tr REPEAT_THRESHOLD_RULE] [-tn NONREPEAT_THRESHOLD_RULE] [-tv EXOVIRUS_THRESHOLD_RULE] [-sr STAR_RATIO] [-bm BOWTIE_MAPQ] [-q Q] [-p P] [-l L] [--verbose]`

Example: `fasticlip -i rawdata/example_MMhur_R1.fastq rawdata/example_MMhur_R2.fastq --GRCm38 -s docs/GRCm38/GRCm38_STAR/ -n MMhur -o results`

Example: `fasticlip -i rawdata/example_Hmhur_R1.fastq rawdata/example_Hmhur_R2.fastq --GRCh38 -s docs/GRCh38/GRCh38_STAR/ -n Hmhur -o results`

### Required arguments

  flag | description
  ------------------|------------------------------------------------
  -h, --help   |      show this help message and exit
  -i INPUT(s) | At least one input FASTQ (or fastq.gz) files; separated by spaces
  --GRCh38        |    required if your CLIP is from human
  --GRCm38          |   required if your CLIP is from mouse
  -s STAR_INDEX  |   Path to STAR index for your organism
  -n NAME         |  Name of output directory
  -o OUTPUT        | Name of directory where output directory will be made
  
### Optional arguments

  flag | description
  ------------------|------------------------------------------------
  --trimmed       |      flag if files are already trimmed
  -f N          |    Number of bases to trim from 5' end of each read. Default is 14. If using irCLIP RT primers, this value should be 18.
  -a ADAPTER     |   3' adapter to trim from the end of each read. Default is AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG.
  -tr REPEAT_THRESHOLD_RULE | m,n: at least m samples must each have at least n RT stops mapped to repeat RNAs. Default is 1,4 (1 sample); 2,3 (2 samples); x,2 (x>2 samples)
  -tv EXOVIRAL_THRESHOLD_RULE | m,n: at least m samples must each have at least n RT stops mapped to viral genome. Default is 1,4 (1 sample); 2,3 (2 samples); x,2 (x>2 samples)
  -tn NONREPEAT_THRESHOLD_RULE | m,n: at least m samples must each have at least n RT stops mapped to nonrepeat RNAs. Default is 1,4 (1 sample); 2,3 (2 samples); x,2 (x>2 samples)
  -sr STAR_RATIO   |     Maximum mismatches per base allowed for STAR genome mapping (corresponds to outFilterMismatchNoverLmax). Default is 0.08 (2 mismatches per 25 mapped bases).
  -bm BOWTIE_MAPQ   |    Minimum MAPQ (Bowtie alignment to repeat/tRNA/retroviral indexes) score allowed. Default is 42.  
  -q Q         |         Minimum average quality score allowed during read filtering. Default is 25.
  -p P          |    Percentage of bases that must have quality > q during filtering. Default is 80.
  -l L           |       Minimum length of read. Default is 15.
  --verbose	| Prints out lots of things :)


Installation instructions
------------

1. Clone this repository by running one of the following:
	- `git clone git@github.com:ChangLab/FAST-iCLIP.git` if you use ssh authentication
	- `git clone https://github.com/ChangLab/FAST-iCLIP.git` otherwise
2. Type `cd FAST-iCLIP` to enter the folder.
3. Run `./configure`. This will check for dependencies (below) and download necessary files (bowtie indices, gene lists and genomes, and example iCLIP data).
4. Run `sudo python setup.py install`. If you do not have sudo privileges, run `python setup.py install --user` or `python setup.py install --prefix=<desired directory>`.
5. You should see three new folders inside `FAST-iCLIP`: `docs`, `rawdata`, and `results`.
6. Try running the following command: 
  `fasticlip -i rawdata/example_MMhur_R1.fastq rawdata/example_MMhur_R2.fastq --GRCm38 -s docs/GRCm38/GRCm38_STAR/ -n MMhur -o results`. It should run in ~1 hour. Look inside `results/MMhur` for output files.

You can `fasticlip` from outside its installation directory. To do this, add the following lines to the end of your `.bash_profile` script:

	export FASTICLIP_PATH=<your absolute FAST-iCLIP installation path>
	export PATH=$FASTICLIP_PATH:$PATH
Save the file, then run `source ~/.bash_profile`.

Dependencies
------------

The version numbers listed have been tested successfully. There can be difficulties if you choose to run updated versions of some of these dependencies.

- Python 2.7: https://www.python.org/download/releases/2.7/ (Important: does not work with Python3)
- matplotlib 1.5: http://matplotlib.org/
- Pandas 0.18.1: http://pandas.pydata.org/
- Bowtie 2.1: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
- STAR 2.4.0: https://github.com/alexdobin/STAR
- bedtools 2.25.0: http://bedtools.readthedocs.org/en/latest/
- FASTX-Tookit 0.0.13: http://hannonlab.cshl.edu/fastx_toolkit/
- matplotlib-venn 0.11.4: https://pypi.python.org/pypi/matplotlib-venn
- iCLIPro 0.1.1: http://www.biolab.si/iCLIPro/doc/

Input
-----

At least one FASTQ or compressed FASTQ (fastq.gz). Use the `--trimmed` flag if trimming has already been done.

Output
------

Three subdirectories inside the named directory within `results`. 
- `figures` has 6 figures in pdf and png format.
  - *Figure 1* visualizes the some of the relevant summary data.
    - A. Read count summary per pipeline step. The source data is: PlotData_ReadsPerPipeFile
    - B. Bar graph of gene count per RNA type. The source data is: PlotData_ReadAndGeneCountsPerGenetype
    - C. Pie chart of RT stops mapping to known features of mRNAS including 5'UTR, Introns, CDS, and 3'UTR.
    	- This uses reads obtrained from intersection with ENSEMBL-derived UTR coordinates. The source data is: PlotData_ReadsPerGene_*UTR or CDS
    - D. Pie chart of RT stops mapped to all indexes included in the FAST-iCLIP pipeline.

  - *Figure 2* provides coverage histograms of binding across each repeat RNA element, both sense and antisense strands. 
    - Source data: PlotData_RepeatRNAHist_*
    - RT stops mapping to the positive and negative strands are shown in blue and red, respectively.

  - *Figure 3* provides coverage histograms of binding across the rRNA, highlighting mature rRNA regions. 
    - Source data: PlotData_RepeatRNAHist_*
    - RT stops mapping to the positive and negative strands are shown in blue and red, respectively.

  - *Figures 4a and 4b* provide a summary of snoRNA binding data.
    - Histograms display RT stop position within an average snoRNA gene body.
    - The pie chart provides a summary of reads per snoRNA type.

  - *Figure 5* provides histograms of RT stop position within gene body for all remaining ncRNA types.

  - *Figure 6* provides a pie chart composed of RT stops from the top 15 best bound endoVirus elements.
    - Total RT stop counts per element and percentage of the total endoVirus mapped reads are included for each element in the legend.

  - *Figure 7* provides histograms of RT stop position across the genome for any exoViruses (DV, ZV, or HCV).
    - RT stops mapping to the positive and negative strands are shown in blue and red, respectively.

- `rawdata` has all the PlotData files used to make the figures, as well as intermediate files that can be useful in generating other plots.
- `todelete` has files that are unnecessary to keep.

How the pipeline works
----------------------

1. Prepare 2 FASTQ or FASTQ.gz files corresponding to the two replicates for an iCLIP/irCLIP experiment.

2. Duplicate removal, quality filter, and trim adapter from the 3' end from the reads
  - Duplicate removal:
	  - This step takes advantage of the fact that 5' end of each read has a random barcode.
	  - Each initial starting molecule that was RT'd will have a unique barcode.
	  - Therefore, PCR duplicates are removed by collapsing molecules with identical 5' barcode sequences.
  - Remove adapter region from the 3' end of the read.
  - The adapter is an optional input parameter.
  - Default is to remove sequences less than N=33 nucleotides. 
  - Q33 specifies the quality score encoding format.

3. After duplicate removal, remove the 5' barcode sequence. Default removes 13 nts.
  - For libraries made with sequences found in the FAST-iCLIP manuscript, remove 13 nts (-f 14)	
  - For libraries made with sequences found in the irCLIP manuscript, remove 17 nts (-f 18)	

4. We then map the reads to indexes in the following order:
  - A. exoVirus index
  - B. repeat RNA index
  - C. endoVirus index
  - D. tRNA index
  - E. Human or Mouse non-repetitive genome
  	- By default, we only keep reads that are unique and perfectly aligned.
  	- We then remove reads that map to blacklist or repeat regions.

5. After mapping, we isolate the 5' position (RT) stop for both positive and negative strand reads.
  - This represents the cross-link site in the initial experiment.

6. For each replicate, we analyze the RT stop position and read length using iCLIPro. 
  - See Haberman et al 2017 (https://www.ncbi.nlm.nih.gov/pubmed/28093074)
  - iCLIPro analyzes the RT stops vs cDNA length. It has been observed that RBPs and the complexes they bind RNAs in can occupy small or large areas on target RNAs. If RNase digestion is too stringent the identified RT stops via cDNA truncation may be shifted.

7. For each strand, we merge RT stops between replicates. 
  - This means that at RT stop position must be 'conserved' between replicates.
  - If conserved, we count the total number of instances of the RT position for both replicates.
  - If the total counts exceed a specified threshold, then we record these RT stops.
  - We then make bedGraph and BigWig files, allowing visualization in genome browsers. 

8. Partition RT stops by gene type.
  - We intersect merged RT stops with Ensembl annotation of gene name by RNA type.
  - Some protein coding and lincRNA genes can have embedded snoRNAs or miRNAs in their introns, making this more challenging.
  	- In turn, we re-generate the initial RT stop and intersect this with two different filters.
  	- One filter is a snoRNA mask and the other filter is a miR mask: both derived from Ensembl annotations.
  	- These masks allow us to remove all "protein coding" RT stops that fall within annotated snoRNA/miR regions.

9. Quantification of reads per gene.
  - For each gene type, we quantify the number of reads per gene.
  - For all but snoRNAs, this is computed using the bed files obtained above.
  - For snoRNAs, we intersect the initial pool of RT stops with a custom annotation file of snoRNAs built from Ensembl annotations of noncoding RNAs. 
  - Collectively, this gives us reads per gene for each gene type.

10. Partition protein coding reads by functional mRNA elements.
  - We intersect snoRNA/miR filtered reads with Ensembl-derived UTR coordinates. 
  - We perform this such output files generated are both "exclusive" and "inclusive".
  	- Exclusive = a gene will have RT stops only in the specific element, without any RT stops elsewhere in that gene.
  	- Inclusive = a gene will have RT stops in the specific element, but can also have RT stops elsewhere in that gene.

11. Partition reads by ncRNA binding region.
  - For non-coding RNAs, we simply annotate reads with the start and stop position for each ncRNA.
  - This allows us to determine the position of each RT stop with respect to the full length of the gene.

12. Partition repeat-mapped RT stops by region.
  - The repeat RNA mapped RT stops are paritioned using the repeat custom index annotation. 
  - As with the ncRNAs, this is later used for visualization. Both sense and antisense mapped RT stops are visualized.
