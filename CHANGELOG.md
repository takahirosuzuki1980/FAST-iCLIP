# Changelog

## Unreleased / Known issues
- can run the other scripts (viralgenome, motifAnalysis) outside the installation directory
- Need to move all global variables to cfg.py
- Add hg38 and mm10 annotation files
- docs.tar.gz needs to be updated
- Add more analysis of retroviral data (currently in another script)

## [0.9.3] - 2016-02-07
### Added
- New flag to control STAR N over L ratio (mismatches per base). Default set at 0.08

### Fixed
- Thresholding now implemented correctly for tRNA, repeat, and genome

## [0.9.2] - 2016-02-02
### Added
- Maps to retroviral index now

### Fixed
- ReadsPerPipeFile now outputs counts

## [0.9.1] - 2016-01-31
### Fixed
- Script runs to completion with --clipper
- Figure 1 is no longer stretched out vertically
- When not using CLIPper, protein-coding genes are now correctly enumerated
- Folder structure is now correct

## [0.9.0] - 2016-01-29
### Added
- Installation procedure using setup.py
- Now runs STAR instead of bowtie2 for genome index mapping. STAR is faster and maps across splice junctions.
- Users now have to supply a STAR index (-s flag required)
- CLIPper is now optional (and no longer recommended). Instead, reads are now overlapped by snoRNAs and miRNAs before overlapped to genes.
- Logging is improved, with a --verbose flag to see all command lines run
- Histograms of tRNA isotypes are now outputted, along with comparison across replicates

### Fixed
- Files are now outputted in a more manageable directory tree
- Code split into several files to improve maintainability

## 2016-01-16
### Added
- tRNA indexes now end with CCA modification
- iCLIPro is now run on output

### Fixed
- The 5' end of RT stops on the minus strand is now correctly reported (off-by-one error)

## 2015-12-26
### Added
- Both + and - signal are now reported for each nucleotide in PlotData_RepeatRNAHist_*
- Mean and standard deviation for each significant RT stop are now reported in files that end with _mergedRT.bed_stats

## 2015-12-25
### Added
- Added the ability to run fasticlip outside its installation directory once the appropriate environment variables are set

## 2015-12-23
### Added
- Can now handle any number of replicates. Added default thresholds for RT stop significance for different numbers of replicates.
- Users can define significance rules for RT stops for repeat and non-repeat RNAs separately with -tr and -tn

## 2015-12-16
### Added
- ViralGenome_analysis.py added; can now map to arbitrary viral genome Bowtie2 indexes

## 2015-09-24
### Added
- Can now specify filename prefix in motifAnalysis.py

### Fixed
- In Figure 5, if PlotData files are empty the plots will not be drawn (rather than crashing)

## 2015-09-23
### Fixed
- Clarified the definition of exons, CDS, introns, and UTRs in both fasticlip and motifAnalysis.py
- Added the file PlotData_ReadsPerGene_proteinCoding_byRegion to count up reads for each of the above regions for each gene
- Removed -f value from the bedtools intersect commands, since default is -f 1e-9 which should work

## 2015-09-20
### Added
- Users can now specify the minimum MAPQ quality they want to keep using the -q flag

## 2015-09-14
### Added
- mapping to tRNA indexes (docs.tar.gz also updated to add tRNA indexes)
- tRNA and repeat masker reads are now reported

### Fixed
- Default lower limit read length is 15
- Duplicate removal of reads is now done before 3'/5' trimming to save time
- Various redundancy fixes to speed up making Figures 1, 3, 4
- Using the -sorted flag in bedtools intersect to speed up the repeat and blacklist masking

## 2015-07-13
### Added
- Adding a list of genes that have 5' and 3' UTR binding, but no CDS binding

## 2015-06-12
### Added
- In motifAnalysis.py, lists of RT stops are now filtered so that there is only one representative per significant "window"

## 2015-05-27
### Fixed
- rDNA histogram is now the correct length

## 2015-04-12
### Fixed
- CLIP clusters are now merged if 1) they are duplicates of each other, or 2) they overlap

## 2015-02-26
### Changed
- changed location of docs.tar.gz and updated ./configure

## 2015-01-09
### Added
- added a script to map to retroviral indexes
- Added the necessary retroviral fasta files to the docs.tar.gz file

## 2014-12-30
### Added
- added motifAnalysis.py to analyze motifs using HOMER

## 2014-12-21
### Fixed
- we now combine all isoforms of each gene during plot making

## 2014-12-04
### Added
- Added a minimum read length post-5' trim

## 2014-11-09
### Fixed
- Changed -f13 to -f14 for the 5' trim so that the first 13 nucleotides are trimmed by default. -f should represent the first nucleotide that is kept











