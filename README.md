FAST-iCLIP
==========

Fully Automated and Standardized iCLIP (FAST-iCLIP) is a fully automated tool to process iCLIP data.

This package contains two main sets of tools: an executable called `fasticlip` to run iCLIP on human and mouse data, and several iPython notebooks to process iCLIP data from viral genomes.

The following README will focus mainly on `fasticlip`. The pdf in the repository contains further instructions for using the iPython notebooks.

**Table of Contents**

- [Usage](#user-content-usage)
	- [Required arguments](#user-content-required-arguments)
	- [Optional arguments](#user-content-optional-arguments)
- [Installation instructions](#user-content-installation-instructions)
- [Dependencies](#user-content-dependencies)
- [Input](#user-content-input)
- [Output](#user-content-output)
- [Debugging](#user-content-debugging)
- [How the pipeline works](#user-content-how-the-pipeline-works)

Usage
-----

`fasticlip [-h] -i INPUT1 INPUT2 [--hg19 | --mm9] -n NAME -o OUTPUT [-f N] [-a ADAPTER] [-t THRESHOLD] [-q Q] [-p P]`

Example: `fasticlip -i rawdata/example_MMhur_R1.fastq rawdata/example_MMhur_R2.fastq --mm9 -n MMhur -o results`

### Required arguments

  flag | description
  ------------------|------------------------------------------------
  -h, --help   |      show this help message and exit
  -i INPUT1 INPUT2  | 2 input FASTQ files separated by a space
  --hg19        |    required if your CLIP is from human
  --mm9          |   required if your CLIP is from mouse
  -n NAME         |  Name of output directory
  -o OUTPUT        | Name of directory where output directory will be made
  
### Optional arguments

  flag | description
  ------------------|------------------------------------------------
  -f N          |    Number of bases to trim from 5' end of each read. Default is 13.
  -a ADAPTER     |   3' adapter to trim from the end of each read. Default is A            GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG.
  -t THRESHOLD    |  Stringency of RT stop filtering, where higher numbers indicate greater stringency. Default is 3.
  -q Q             | Minimum quality score to keep during filtering. Default is 25.
  -p P          |    Percentage of bases that must have quality > q during filtering. Default is 80.


Installation instructions
------------

1. Clone this repository by running `git clone git@github.com:ChangLab/FAST-iCLIP.git` if you use ssh authentication or `git clone https://github.com/ChangLab/FAST-iCLIP.git` otherwise.
2. `cd` into the `FAST-iCLIP` folder.
3. To install, run `./configure`. This will check for dependencies (below) and download necessary files (bowtie indices, gene lists and genomes, and example iCLIP data).
4. You should see three new folders inside `FAST-iCLIP`: `docs`, `rawdata`, and `results`.
5. Try running the following command: 
  `fasticlip -i rawdata/example_MMhur_R1.fastq rawdata/example_MMhur_R2.fastq --mm9 -n MMhur -o results`. It should run fairly quickly. Look inside `results/MMhur` for output files.

Dependencies
------------

1. Python 2.7: https://www.python.org/download/releases/2.7/
2. iPython: http://ipython.org/install.html (optional)
3. iPython notebook: http://ipython.org/notebook (optional)
4. Matplotlib for Python (plotting): http://matplotlib.org/
5. Pandas for Python (data): http://pandas.pydata.org/
6. Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
7. bedtools: http://bedtools.readthedocs.org/en/latest/
8. CLIPper: https://github.com/YeoLab/clipper/wiki/CLIPper-Home
9. FASTX-Tookit: http://hannonlab.cshl.edu/fastx_toolkit/

Input
-----

There must be two replicates of your iCLIP data in FASTQ format. These FASTQ files should NOT be trimmed or processed because trimming and processing will occur as part of the pipeline.

Output
------

Three subdirectories inside the named directory within `results`. 
- `figures` has 6 figures in pdf and png format.
  - *Figure 1* visualizes the some of the relevant summary data.
    - It includes a read count summary per pipeline step. The source data is: PlotData_ReadsPerPipeFile
    - It includes a pie chart of UTR binding.
    - This uses reads obtrained from intersection with ENSEMBL-derived UTR coordinates. 
    - The source data is: PlotData_ReadsPerGene_*UTR or CDS
    - It also includes a bar graph of gene count per RNA type.
    - The source data is: PlotData_ReadAndGeneCountsPerGenetype
 
  - *Figure 2* provides a richer summary of the UTR data.
    - The upper panel is an aggregate trace of binding derived from a custom perl script.
    - The lower panels provide a heatmap of binding intensity for gene exclusivly bound in each UTR or CDS.
    - This allows us to isolate genes with exclusive UTR,CDS,or intronic binding.
    - The source data is: PlotData_ExclusiveBound_*

  - *Figures 3 and 4* provide coverage histograms of binding across each repeat RNA and rRNA, respectivly. 
    - Source data: PlotData_RepeatRNAHist_*

  - *Figure 5* provides a summary of snoRNA binding data.
    - The pie chart provides a summary of reads per snoRNA type.
    - These are complimented by histograms of RT stop position within the snoRNA gene body.

  - *Figure 6* provides histograms of RT stop position within gene body for all remaining ncRNA types.

- `rawdata` has all the PlotData files used to make the figures, as well as intermediate files that can be useful in generating other plots.
- `todelete` has files that are unnecessary to keep.

Debugging
----------

1. _Mapping:_ Ensure that bowtie2 is in the $PATH and executable.
2. _CLIPper:_ Uses Python 2.7. 
3. _Any script in bin/:_
    - The provided BedGraphToBigWig is built for Linux. 
  - This, and related scripts, may be downloaded for other platforms: http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.i386/
4. If there is a problem parsing clusters from CLIPper, consider the version used.
  - The older version of CLIPper gives files with name `<name>_<clusterNum>_<readPerCluster>`
  - The newer version gives files with name `name.val__<clusterNum>_<readPerCluster>`

How the pipeline works
----------------------

1. Prepare 2 FASTQ files corresponding to the two replicates for the iCLIP experiment.

2. Trim adapter from the 3' end from the reads.
  - RT primer is cleaved, leaving adapters. 
  - Remove adapter region from the 3' end of the read.
  - The adapter is an optional input parameter.
  - Default is to remove sequnces less than N=33 reads. 
  - Q33 specifies the quality score encoding format.

3. Quality filter.

4. Remove duplicates.
  - This step takes advantage of the fact that 5' end of each read has a random barcode.
  - Each initial starting molecule that was RT'd will have a unique barcode.
  - Therefore, PCR duplicates are removed by collapsing molecules with identical 5' barcode sequences.

5. After duplicate removal, remove the 5' barcode sequence.

6. We then map to a repeat index.
  - We use k=1, meaning bt2 will search for 1 distinct, valid alignments for each read.
  - This step allows us to both remove reads that are normally blacklisted and also map reads to the repeat index.
  - The repeat index is derived from UCSC Genome Browser.

7. After mapping, we isolate the 5' position (RT) stop for both positive and negative strand reads.
  - This represents the cross-link site in the initial experiment.

8. For each strand, we merge RT stops between replicates. 
  - This means that at RT stop position must be conserved between replicates.
  - If conserved, we count the total number of instances of the RT position for both replicates.
  - If the total counts exceed a specified threshold, then we record these RT stops.
  - Finally, we re-generate a "read" around the RT stop using the passed parameter "expand," 
  - A "read" around the RT stop is required for downstream processing.

9. Reads that do not map to the repeat index are then mapped to hg19.
  - The mapped reads are processed by samtools, repeat masker, and blacklist filter.
  - As with the repeat index, we then merge RT stops.

10. Expanded reads from RT stop merging are passed to CLIPper, a peak calling algorithm.
  - CLIPper returns a bed-like file format with window coordinates, reads counted per window, etc.
  - We use these windows to extract "low FDR" reads from the total set of reads passed to CLIPper.
  - We then make bedGraph and BigWig files from this complete pool of "low FDR" reads, allowing easy visualization. 

11. Partition "low FDR" reads by gene type.
  - We partition the gene names recoved from CLIPper using ENSEMBL annotation of gene name by RNA type.
  - Once this is done, we also split the "low FDR" reads recovered from CLIPper by type using the gene name.
  - Protein coding and lincRNA genes can be embedded snoRNAs or miRNAs that make this more challenging.
  - In turn, we re-generate the initial RT stop and intersect this with two different filters.
  - One filter is a snoRNA mask and the other filter is a miR mask. 
  - Both are derived from Ensembl annotations.
  - These masks allow us to remove all "protein coding" RT stops that fall within annotated sno/mi-RNA regions.

12. Quantification of reads per gene.
  - For each gene type, we quantify the number of reads per gene.
  - For all but snoRNAs, this is computed using the bed files obtained above.
  - For snoRNAs, we intersect the initial pool of "low FDR" reads with custom annotation file. 
  - The custom annotation file is a summation of UCSC, Rfam, and ENSEMBL loci annotated as snoRNAs as no individual annotation entirely captured all known loci.
  - Collectively, this gives us reads per gene for each gene type.
  - All are based upon ENSEMBL annotation except for the snoRNAs.

13. Summary of RT stop intensity around CLIPper cluster centers.
  - We generate a bed file of cluster center positions using the CLIPper cluster file output.
  - We use a custom perl script that generates a heatmap of RT stop intensty per cluster.
  - This allows us to later visualize the distribution of RT stops per cluster.

14. Parition protein coding reads by UTR.
  - We intersect sno/mi-RNA filtered reads with ENSEMBL-derived UTR coordinates. 
  - We perform this such that each read assignment is mutually exclusive.
  - This only isolates reads that fall within each UTR type.
  - Similarly, we use a custom perl script to generate a matrix of read intensity per gene.
  - This provides a complete binding profile per gene.

15. Partition reads by ncRNA binding region.
  - For non-coding RNAs, we simply annotate reads with the start and stop position for each ncRNA.
  - This allows us to determine the position of each RT stop with respect to the full length of the gene.

16. Partition repeat-mapped RT stops by region.
  - The repeat RNA mapped RT stops are paritioned using the repeat custom index annotation. 
  - As with the ncRNAs, this is later used for visualization.
