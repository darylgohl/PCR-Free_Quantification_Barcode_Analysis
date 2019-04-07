# REcount Barcode Analysis
Code for analyzing REcount barcode sequencing data. The REcount_barcode_analysis.py program takes in a FASTQ file and barcode reference FASTA file and outputs a plot and a .txt file with barcode counts.

## Prerequisites
Python 2.7

BioPython

## Usage
REcount_barcode_analysis.py [-h] [-i] [-r] [-o] [-m]

optional arguments:

      -h, --help                show this help message and exit

      -i , --input_file         Input path for raw FASTQ file [required].

      -r , --reference_file     Input reference FATSA file [required].

      -o , --output_dir         Output directory for plot and barcode count file (default: same folder as input file)

      -m , --mismatches_allowed  Number of mismatches to barcode reference sequences allowed (default 2)
  
## Usage example
REcount_barcode_analysis.py -i <PathToFile/InputFileName> -r <PathToFile/ReferenceFileName>
