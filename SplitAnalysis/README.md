# Introduction
Cell free eccDNA codes (cfeccDNACodes) contains scripts for linearized cell-free eccDNA (cf-eccDNA) fragments analysis, cf-eccDNA feature extraction and cf-eccDNA feature-based cancer detection models.
The functions provided in cf-eccDNACodes are illustrated as follows.
![codes](https://github.com/user-attachments/assets/5770489b-7ad9-4814-b872-ed6bbe63fa00)



## Getting Started
## Hardware requirements
Enough RAM to support the in-memory operations.

### OS Requirements
This package is supported for *Linux* and *macOS*. The package has been tested on the following systems:
+ Linux: Ubuntu 16.04.7
+ macOS: ventura 13.2.1

### Prerequisites
Please ensure the following dependencies are installed:

**Python Packages:**
- Python (>=3.7)
- numpy
- pandas
- scikit-learn

**R Packages:**
- R (>=4.0)
- Bioconductor: BSgenome, BSgenome.Hsapiens.UCSC.hg38, GenomeInfoDb, GenomicAlignments, GenomicRanges, IRanges, Rsamtools
- CRAN: dplyr, mgsub, readxl, stringr, tidyr

**Other:**
- Perl
	
 
### Installation
1. **Necessary Step:** Download from Github:
```
git clone https://github.com/xqpeng/cfeccDNA_H-DNA_prediction.git
cd cfeccDNA_H-DNA_prediction
```
2. **Recommended Step:** (Conda users, Conda version: 4.12.0) Create your environment and activate it:
```
conda env create -f environment.yaml
conda activate cfeccDNA_analysis
Rscript setup.R
```

## Usage
1. **sh map.sh *<input_dir>* *<output_dir>* *<genome_index_dir>***
   
   Function: Align the paired fastq files from the file directory *input_dir*, and output the sorted bam files to the file directory *output_dir*.
   
   Make sure you have built indexed reference genome and provided it in *genome_index_dir*.
   
   Input parameters:
   
          (1) <input_dir>: the file directory containing the paired fastq files
   
          (2) <output_dir>: the file directory to store the sorted bam files
   
          (3) <genome_index_dir>: the directory containing the indexed reference genome
   
    Output: Sorted bam files


2. **sh runSplit_read.sh *<input_dir>* *<output_dir>***
   
   Function: Extract split-aligned fragments from sorted bam files.
   
   Input parameters:
   
         (1) <input_dir>: the file directory containing the sorted bam files
   
         (2) <output_dir>: the file directory to store the bam files of split-aligned fragments
   
   Output: The bam files of split-aligned fragments
   
3. **sh runDiscordant_read.sh *<input_dir>* *<output_dir>***

   Function: Extract discordantly-aligned fragments from sorted bam files.
   
   Input parameters:
   
         (1) <input_dir>: the file directory containing the sorted bam files
   
         (2) <output_dir>: the file directory to store the bam files of split-aligned fragments
   
   Output: The bam files of split-aligned fragments
   
4. **sh bam2sam.sh -i *<input_dir>* -o *<output_dir>***

	Function: Convert bam files into sam files
	
	Input parameters:
   
          (1) <input_dir>: the file directory containing the bam files
   
          (2) <output_dir>: the file directory to store the converted sam files
	
	Output: sam files
   
5. **perl SplitAnalysis.pl *<list_file>***

   Function:
   
          (1) Extract two ends of each fragment, following the sequencing direction 5'-3'.
   
          (2) Perform statistic analysis on the split-aligned fragments from each sample.
   
          (3) Perform length analysis on the split-aligned fragments from a group of samples.

          (4) Perform lengths, endpoint distances, and identical 6-mer end motifs analysis on the split-aligned fragments from each group.

   Input parameter:

           (1)  <list_file> is a file contain the split-aligned sam filenames of a group of samples.
   
              For example: list_HCC
   
		      which contains three lines and there is a filename on each line, such as following:
   
				         HCC1.split_reads.sam
   
				         HCC2.split_reads.sam
   
				         HCC3.split_reads.sam
   
                         ...
   Output:

           (1) statics_proportion_length_list_HCC; # corresponding to  function 2
   
           (Belowing are optional)
   
		   (2) File of Fragment length distribution
   
		   (3) File of eccDNA molecules size distribution
   
		   (4) File of endpoint distance distribution of Same-Strand 5' Ends fragments
   
		   (5) File of endpoint distance distribution of Opposite-Strand 5' Ends fragments
   
		   (6) File of Identical 6-mer end motifs profile of Same-Strand 5' Ends fragments
   
		   (7) File of Identical 6-mer end motifs profile of split-aligned fragments fragments
   
	       (8) Files of the 5ends of fragments corresponding to the files in the input file. 
   
6. **perl DiscordantAnalysis.pl  *<list_file>***
   
   Function:
   
         (1) Perform statistic analysis on the discordant-aligned fragments from each sample.
   
         (2) Perform  identical 6-mer end motifs analysis on the discordant-aligned fragments from each group

   Input parameter:
   
         (1) list_file is a file contain the discordantly-aligned sam filenames of a group of samples.
   
             For example: list_HCC_discordant_pairs
   
		      which contains three lines and there is a filename on each line, such as following:
   
				         HCC1_discordant_pairs.sam
   
				         HCC2_discordant_pairs.sam
   
				         HCC3_discordant_pairs.sam
                         ...
	 Output:

             (1) statics_discordant_fragments_list_HCC_discordant_pairs;

             (2) File of Identical 6-mer end motifs profile of discordantly-aligned fragments
   

