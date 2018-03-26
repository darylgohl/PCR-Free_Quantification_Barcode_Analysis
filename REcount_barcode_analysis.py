#!/usr/bin/env python

"""
Copyright 2018  Daryl Gohl

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from Bio import SeqIO
import regex
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
__version__ = "1.0"

def get_args(x):
    x.add_argument("-i", "--input_file",
                   type = str,
                   default = None,
                   metavar = '',
                   help="Input path for raw FASTQ file [required].")
    x.add_argument("-r", "--reference_file",
                   type = str,
                   default = None,
                   metavar = '',
                   help="Input reference FATSA file [required].")
    x.add_argument("-o", "--output_dir",
                   type = str,
                   default = '',
                   metavar = '',
                   help = "Output directory for plot and barcode count file (default: same folder as input file)") 
    x.add_argument("-m", "--mismatches_allowed",
                   type = str,
                   default = '2',
                   metavar = '',
                   help = "Number of mismatches to barcode reference sequences allowed (default 2)")  
    args = x.parse_args()
    return args

# Parses command line arguments
argparser = argparse.ArgumentParser(description = "Barcode counting script (v" + __version__ +")\n" + \
                                              "by Daryl Gohl\n" + \
                                              "This program takes in a FASTQ file and barcode reference FASTA file and outputs plots and a .txt file of barcode counts.",
                                add_help = True, 
                                epilog ='')
args = get_args(argparser)
args = vars(args)

# possible input args
filename = args['input_file']
file = os.path.split(filename)[1]
Ref_filename = args['reference_file']
mismatches = args['mismatches_allowed']
out_folder = args['output_dir']

if out_folder == '':
  out_dir = os.path.dirname(filename)
else:
  out_dir = out_folder

#Count number of records in the file
count = 0
for record in SeqIO.parse(filename, "fastq"):
	count += 1
print("There were " + str(count) + " records in file " + filename)
total_recs = count

#Count number of records in the file
count = 0
for record in SeqIO.parse(Ref_filename, "fasta"):
	count += 1
Ref_rec = count
print("There were " + str(Ref_rec) + " records in the reference database")

#Count up barcodes in fastq file
bc_ID_list = [] #Construct name
bc_seq_list = [] #Barcode sequence
count_list = [] #Barcode counts
bc_all_list = []
for record in SeqIO.parse(Ref_filename, "fasta"):
    count = 0
    bc_sub_list = []
    bc_ID = record.id #Collect standard IDs from reference file
    bc_seq = str(record.seq) #Collect standard barcode sequences from reference file
    for i in SeqIO.parse(filename, "fastq"):
        query = r'(?:' + bc_seq +'){s<=' + mismatches + '}' #fuzzy matching - allow up to <arg> mismatches
        test = regex.findall(query, str(i.seq[:20]))
        if test != []: 
            count += 1
            bc_sub_list.append(str(i.seq[:20]))   
    count_list.append(count)
    bc_all_list.append(bc_sub_list)
    bc_ID_list.append(bc_ID)
    bc_seq_list.append(bc_seq)
    print("done with " + bc_ID) 

account = float(np.sum(count_list))
perc = account/total_recs*100

#Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(len(count_list))
ax1 = ax.bar(ind,count_list, width=0.8, bottom=None)
width=0.8
ax.set_ylabel('counts')
ax.set_title('Barcode counts')
ax.set_xticks(ind+width/2)
xtickNames = ax.set_xticklabels(bc_ID_list)
plt.setp(xtickNames, rotation=90, fontsize=8)
fig_name = os.path.join(out_dir,(file[:-6] + "_barcode_counts.png"))
plt.savefig(fig_name,bbox_inches='tight')

#Write out barcode count file
save_name = os.path.join(out_dir,(file[:-6] + "_barcode_counts.txt"))
save_file = open(save_name, "w")
header = ("Name",'\t',"Barcode", '\t',"Counts",'\n')
save_file.write(''.join(map(str, header)))
newtab = '\t'
newline = '\n'
for i, item in enumerate(count_list):
	save_file.write(bc_ID_list[i])
	save_file.write(newtab)
	save_file.write(bc_seq_list[i])
	save_file.write(newtab)
	save_file.write(str(item))
	save_file.write(newline)
save_file.close()
print("There were " + str(account) + "/" + str(total_recs) + "(" + str(perc) + "%) records accounted for.")
