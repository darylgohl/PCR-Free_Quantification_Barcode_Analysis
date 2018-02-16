# PCR-Free_Quantification_Barcode_Analysis
Code for analyzing PCR-free quantification barcode sequencing data. This program takes in a FASTQ file and barcode reference FASTA file and outputs a plot and a .txt file of barcode counts.

## Prerequisites
Python 2.7

BioPython

## Usage
PCR_free_barcode_analysis.py [-h] [-i] [-r] [-o] [-m]

optional arguments:

      -h, --help                show this help message and exit

      -i , --input_file         Input path for raw FASTQ file [required].

      -r , --reference_file     Input reference FATSA file [required].

      -o , --output_dir         Output directory for plot and barcode count file (default: same folder as input file)

      -m , --mismatches_allowed  Number of mismatches to barcode reference sequences allowed (default 2)
  
## Usage example
PCR_free_barcode_analysis.py -i <PathToFile/InputFileName> -r <PathToFile/ReferenceFileName>
