from Bio import SeqIO
import regex
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import gzip

#Input args
#filename = sys.argv[1] #Input fastq file
#Ref_filename = sys.argv[2] #Barcode reference file

filename = "<PATH_TO_FASTQ_FILE>"

Ref_filename = "<PATH_TO_UMGC_423_Variable.fasta>"

folder = "<OUTPUT_DIRECTORY>"

instrument = "<INSTRUMENT>"

if not os.path.exists(folder):
    os.makedirs(folder)

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
qual_list = []
for record in SeqIO.parse(Ref_filename, "fasta"):
    count = 0
    bc_sub_list = []
    temp_fastq_list = []
    R_temp = []
    bc_ID = record.id #Collect standard IDs from reference file
    bc_seq = str(record.seq) #Collect standard barcode sequences from reference file
    for i in SeqIO.parse(filename, "fastq"):
        query = r'(?:' + bc_seq +'){s<=2}' #fuzzy matching - allow up to 2 mismatches
        test = regex.findall(query, str(i.seq[:20]))
        if test != []: #Comment for exact matching
        #if str(i.seq[:20]) == bc_seq: #Uncomment for exact matching
            count += 1
            #R_temp.append(i.letter_annotations["phred_quality"])
            bc_sub_list.append(str(i.seq[:20])) 
            temp_fastq_list.append(i)  
    out_file_name = folder + bc_ID + ".fastq"
    SeqIO.write(temp_fastq_list, out_file_name, "fastq")  
    #qual_list.append(R_temp)  
    count_list.append(count)
    bc_all_list.append(bc_sub_list)
    bc_ID_list.append(bc_ID)
    bc_seq_list.append(bc_seq)
    print("done with " + bc_ID) 
 
os.chdir(folder)
file_names = os.listdir(folder)

data_file_names = []
for i in file_names:
    if i[-6:] ==".fastq":
        data_file_names.append(i)

for i in data_file_names:
	R1 = i[:-6] + "_trimmed.fastq"
	execute = "cutadapt -l 50 " + i + " > " + R1
	os.system(execute)

file_names = os.listdir(folder)

data_file_names = []
for i in file_names:
    if i[-13:] =="trimmed.fastq":
        data_file_names.append(i)
   
data_file_names.sort() 
#Extract q-score and read number information from each sample

#Make data lists
out_dir = folder
full_name = []
fname = []
read_num = []
n_reads = []
Q_mean_by_base = []
Q_stdev_by_base = []
Q_mean_overall = []
Q_stdev_overall = []
Q_all = []
for i, item in enumerate(data_file_names):
    R_temp = []
    counts = 0
    for j,record in enumerate(SeqIO.parse(item, "fastq")):
        R_temp.append(record.letter_annotations["phred_quality"])
        counts += 1
    full_name.append(item)
    fname.append(item)
    read_num.append('1')
    n_reads.append(counts)    
    a = np.array(R_temp)
    Q_all.append(a)
    Q_mean_bb = np.mean(a, axis=0)
    Q_mean_by_base.append(Q_mean_bb)
    Q_stdev_bb = np.std(a, axis=0)
    Q_stdev_by_base.append(Q_stdev_bb)
    Q_mean_o = np.mean(Q_mean_bb)
    Q_mean_overall.append(Q_mean_o)
    Q_stdev_o = np.std(Q_stdev_bb)
    Q_stdev_overall.append(Q_stdev_o)
    #print "done with %s" % item    

#Make separate lists for R1 and R2
R1_full_name = []
R1_fname = []
R1_read_num = []
R1_n_reads = []
R1_Q_mean_by_base = []
R1_Q_stdev_by_base = []
R1_Q_mean_overall = []
R1_Q_stdev_overall = []
R2_full_name = []
R2_fname = []
R2_read_num = []
R2_n_reads = []
R2_Q_mean_by_base = []
R2_Q_stdev_by_base = []
R2_Q_mean_overall = []
R2_Q_stdev_overall = []
for i, item in enumerate(read_num):
    if item == '1':
        R1_full_name.append(full_name[i])
        R1_fname.append(fname[i][:-14])
        R1_read_num.append(read_num[i])
        R1_n_reads.append(n_reads[i])
        R1_Q_mean_by_base.append(Q_mean_by_base[i])
        R1_Q_stdev_by_base.append(Q_stdev_by_base[i])
        R1_Q_mean_overall.append(Q_mean_overall[i])
        R1_Q_stdev_overall.append(Q_stdev_overall[i])        
    elif item == '2':
        R2_full_name.append(full_name[i])
        R2_fname.append(fname[i])
        R2_read_num.append(read_num[i])
        R2_n_reads.append(n_reads[i])
        R2_Q_mean_by_base.append(Q_mean_by_base[i])
        R2_Q_stdev_by_base.append(Q_stdev_by_base[i])
        R2_Q_mean_overall.append(Q_mean_overall[i])
        R2_Q_stdev_overall.append(Q_stdev_overall[i])

os.chdir(out_dir)

#Plot q-scores by sample
figure_width = (len(R1_fname)/48)*6
if figure_width<12:
    figure_width=12
fig = plt.figure(figsize=(figure_width,8))
#plt.style.use('classic')
ax = fig.add_subplot(111)
#plt.title('Average Q-score by sample')
ax.set_ylim(0, 45)
ax.set_ylabel('Mean Q-score')
x = range(len(R1_fname))
x2 = [y+0.5 for y in x]
x3 = [y+0.25 for y in x]
ax.errorbar(x,R1_Q_mean_overall,yerr=[R1_Q_stdev_overall,R1_Q_stdev_overall], fmt='o', color='black', ecolor='lightgray', elinewidth=3, capsize=0, label = "Read 1")
plt.xticks(x, R1_fname, rotation='vertical', fontsize=12)
if len(R2_full_name) != 0:
    ax.errorbar(x2,R2_Q_mean_overall,yerr=[R2_Q_stdev_overall,R2_Q_stdev_overall], fmt='o', color='red', ecolor='lightgray', elinewidth=3, capsize=0, label = "Read 2")
    plt.xticks(x3, R1_fname, rotation='vertical', fontsize=12)
#plt.margins(figure_width/24000.0)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#plt.subplots_adjust(bottom=0.5)
#plt.show()
plt.tight_layout()
#plt.legend(numpoints=1, frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('Q_score_by_sample', bbox_inches='tight', format='png')

#Box plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylabel('Q-score')
ax.boxplot(Q_all)
xtickNames = plt.setp(ax, xticklabels=R1_fname)
locs, labels = plt.xticks()
plt.setp(labels, rotation=90, fontsize=8)
plt.gcf().subplots_adjust(bottom=0.26)
#plt.show()
plt.savefig('Q_score_Boxplot', format='png')


#Plot read numbers by sample
figure_width = (len(R1_fname)/48)*6
if figure_width<12:
    figure_width=12
fig = plt.figure(figsize=(figure_width,8))#plt.style.use('classic')
ax = fig.add_subplot(111)
#plt.title('Average Q-score by sample')
ax.set_ylabel('Number of reads')
x = range(len(R1_fname))
x2 = [y+0.5 for y in x]
ax.bar(x,R1_n_reads, color='black')
plt.xticks(x2, R1_fname, rotation='vertical', fontsize=8)
#plt.margins(0.05,0)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#plt.show()
plt.tight_layout()
plt.savefig('Read_number_by_sample', bbox_inches='tight', format='png')

#Q-score heatmap
import seaborn as sns

# build the figure instance with the desired height
# Two subplots, unpack the axes array immediately
if len(R2_full_name) != 0:
    sns.set_style("whitegrid")
    grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
    figure_height = (len(R1_fname)/48)*6
    if figure_height<12:
        figure_height=12
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,gridspec_kw=grid_kws, figsize=(12,figure_height))
    ax1 = sns.heatmap(R1_Q_mean_by_base, ax=ax1, cmap="RdYlGn",cbar=True, vmin=0)
    ax2 = sns.heatmap(R2_Q_mean_by_base, ax=ax2, cmap="RdYlGn",cbar=True, vmin=0)
    ax1.set_xlabel('Position (read 1)')
    ax1.set_ylabel('Samples')
    names = R1_fname[::-1]
    ax1.set_yticklabels(names, fontsize=6, rotation="horizontal")
    ax1.axes.xaxis.set_ticklabels([])
    #ax1.axes.yaxis.set_ticklabels([])
    ax2.set_xlabel('Position (read 2)')
    ax2.axes.xaxis.set_ticklabels([])
    # let seaborn do it's thing
    #ax = sns.heatmap(R1_Q_mean_by_base, ax=ax, cmap="RdYlGn")
    #sns.heatmap(R1_Q_mean_by_base)
    plt.savefig('Qscore_heatmap', bbox_inches='tight', format='png', dpi=300)
else:
    sns.set_style("whitegrid")
    grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
    fig, ax1 = plt.subplots(1, 1, sharey=True)#,gridspec_kw=grid_kws)
    ax1 = sns.heatmap(R1_Q_mean_by_base, ax=ax1, cmap="RdYlGn",cbar=True)#, vmin=0)
    ax1.set_xlabel('Position (read 1)')
    ax1.set_ylabel('Samples')
    names = R1_fname[::-1]
    ax1.set_yticklabels(names, fontsize=6, rotation="horizontal")
    ax1.axes.xaxis.set_ticklabels([])
    #ax1.axes.yaxis.set_ticklabels([])
    # let seaborn do it's thing
    #ax = sns.heatmap(R1_Q_mean_by_base, ax=ax, cmap="RdYlGn")
    #sns.heatmap(R1_Q_mean_by_base)
    plt.savefig('Qscore_heatmap', bbox_inches='tight', format='png', dpi=300)

save_name = (instrument + "_indiv_Q_score_data.txt")
save_file = open(save_name, "w")
newtab = '\t'
newline = '\n'
save_file.write(instrument)
save_file.write(newtab)
for i in R1_fname:
    save_file.write(i)
    save_file.write(newtab)
save_file.write(newline)
save_file.write("Mean Q score")
save_file.write(newtab)
for i in R1_Q_mean_overall:
    save_file.write(str(i))
    save_file.write(newtab)
save_file.write(newline)
save_file.write("Standard Deviation Q score")
save_file.write(newtab)
for i in R1_Q_stdev_overall:
    save_file.write(str(i))
    save_file.write(newtab)
save_file.write(newline)
save_file.write("Read count")
save_file.write(newtab)
for i in R1_n_reads:
    save_file.write(str(i))
    save_file.write(newtab)
save_file.close()

#Sum all standards for a given size
size_bins = []
for i in data_file_names:
    temp = i.split("_")[3]
    size_bins.append(temp)

unique_sizes = []
for i in size_bins:
    if i not in unique_sizes:
        unique_sizes.append(i)

#concatenate same size files
for i in unique_sizes:
    temp_concat = []
    size_search = "_" + i + "_"
    for j in data_file_names:
        if j.find(size_search) != -1:
            temp_concat.append(j)
    execute = "cat " + temp_concat[0] + " " + temp_concat[1] + " " + temp_concat[2] + " > " + i + "_concat.fastq" 
    os.system(execute)

file_names = os.listdir(folder)
concat_data_files = []
for i in file_names:
    if i[-12:] == 'concat.fastq':
        concat_data_files.append(i)

concat_files_sorted = []
for i in unique_sizes:
    for j in concat_data_files:
        if j.split("_")[0] == i:
            concat_files_sorted.append(j)            

#Make data lists
out_dir = folder
full_name = []
fname = []
read_num = []
n_reads = []
Q_mean_by_base = []
Q_stdev_by_base = []
Q_mean_overall = []
Q_stdev_overall = []
Q_all = []
for i, item in enumerate(concat_files_sorted):
    R_temp = []
    counts = 0
    for j,record in enumerate(SeqIO.parse(item, "fastq")):
        R_temp.append(record.letter_annotations["phred_quality"])
        counts += 1
    full_name.append(item)
    fname.append(item)
    read_num.append('1')
    n_reads.append(counts)    
    a = np.array(R_temp)
    Q_all.append(a)
    Q_mean_bb = np.mean(a, axis=0)
    Q_mean_by_base.append(Q_mean_bb)
    Q_stdev_bb = np.std(a, axis=0)
    Q_stdev_by_base.append(Q_stdev_bb)
    Q_mean_o = np.mean(Q_mean_bb)
    Q_mean_overall.append(Q_mean_o)
    Q_stdev_o = np.std(Q_stdev_bb)
    Q_stdev_overall.append(Q_stdev_o)
    #print "done with %s" % item    

#Make separate lists for R1 and R2
R1_full_name = []
R1_fname = []
R1_read_num = []
R1_n_reads = []
R1_Q_mean_by_base = []
R1_Q_stdev_by_base = []
R1_Q_mean_overall = []
R1_Q_stdev_overall = []
R2_full_name = []
R2_fname = []
R2_read_num = []
R2_n_reads = []
R2_Q_mean_by_base = []
R2_Q_stdev_by_base = []
R2_Q_mean_overall = []
R2_Q_stdev_overall = []
for i, item in enumerate(read_num):
    if item == '1':
        R1_full_name.append(full_name[i])
        R1_fname.append(fname[i][:-14])
        R1_read_num.append(read_num[i])
        R1_n_reads.append(n_reads[i])
        R1_Q_mean_by_base.append(Q_mean_by_base[i])
        R1_Q_stdev_by_base.append(Q_stdev_by_base[i])
        R1_Q_mean_overall.append(Q_mean_overall[i])
        R1_Q_stdev_overall.append(Q_stdev_overall[i])        
    elif item == '2':
        R2_full_name.append(full_name[i])
        R2_fname.append(fname[i])
        R2_read_num.append(read_num[i])
        R2_n_reads.append(n_reads[i])
        R2_Q_mean_by_base.append(Q_mean_by_base[i])
        R2_Q_stdev_by_base.append(Q_stdev_by_base[i])
        R2_Q_mean_overall.append(Q_mean_overall[i])
        R2_Q_stdev_overall.append(Q_stdev_overall[i])

os.chdir(out_dir)

#Plot q-scores by sample
figure_width = (len(unique_sizes)/48)*6
if figure_width<12:
    figure_width=12
fig = plt.figure(figsize=(figure_width,8))
#plt.style.use('classic')
ax = fig.add_subplot(111)
#plt.title('Average Q-score by sample')
ax.set_ylim(25, 40)
ax.set_ylabel('Mean Q-score')
x = range(len(unique_sizes))
x2 = [y+0.5 for y in x]
x3 = [y+0.25 for y in x]
ax.errorbar(x,R1_Q_mean_overall,yerr=[R1_Q_stdev_overall,R1_Q_stdev_overall], fmt='o', color='black', ecolor='lightgray', elinewidth=3, capsize=0, label = "Read 1")
plt.xticks(x, unique_sizes, rotation='vertical', fontsize=12)
if len(R2_full_name) != 0:
    ax.errorbar(x2,R2_Q_mean_overall,yerr=[R2_Q_stdev_overall,R2_Q_stdev_overall], fmt='o', color='red', ecolor='lightgray', elinewidth=3, capsize=0, label = "Read 2")
    plt.xticks(x3, unique_sizes, rotation='vertical', fontsize=12)
#plt.margins(figure_width/24000.0)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#plt.subplots_adjust(bottom=0.5)
#plt.show()
plt.tight_layout()
#plt.legend(numpoints=1, frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('Q_score_by_sample_grouped', bbox_inches='tight', format='png')

#Box plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylabel('Q-score')
ax.boxplot(Q_all)
xtickNames = plt.setp(ax, xticklabels=unique_sizes)
locs, labels = plt.xticks()
plt.setp(labels, rotation=90, fontsize=8)
plt.gcf().subplots_adjust(bottom=0.26)
#plt.show()
plt.savefig('Q_score_Boxplot_grouped', format='png')


#Plot read numbers by sample
figure_width = (len(unique_sizes)/48)*6
if figure_width<12:
    figure_width=12
fig = plt.figure(figsize=(figure_width,8))#plt.style.use('classic')
ax = fig.add_subplot(111)
#plt.title('Average Q-score by sample')
ax.set_ylabel('Number of reads')
x = range(len(unique_sizes))
x2 = [y+0.5 for y in x]
ax.bar(x,R1_n_reads, color='black')
plt.xticks(x2, unique_sizes, rotation='vertical', fontsize=8)
#plt.margins(0.05,0)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#plt.show()
plt.tight_layout()
plt.savefig('Read_number_by_sample_grouped', bbox_inches='tight', format='png')

#Q-score heatmap
import seaborn as sns

# build the figure instance with the desired height
# Two subplots, unpack the axes array immediately
if len(R2_full_name) != 0:
    sns.set_style("whitegrid")
    grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
    figure_height = (len(unique_sizes)/48)*6
    if figure_height<12:
        figure_height=12
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,gridspec_kw=grid_kws, figsize=(12,figure_height))
    ax1 = sns.heatmap(R1_Q_mean_by_base, ax=ax1, cmap="RdYlGn",cbar=True, vmin=0)
    ax2 = sns.heatmap(R2_Q_mean_by_base, ax=ax2, cmap="RdYlGn",cbar=True, vmin=0)
    ax1.set_xlabel('Position (read 1)')
    ax1.set_ylabel('Samples')
    names = unique_sizes[::-1]
    ax1.set_yticklabels(names, fontsize=6, rotation="horizontal")
    ax1.axes.xaxis.set_ticklabels([])
    #ax1.axes.yaxis.set_ticklabels([])
    ax2.set_xlabel('Position (read 2)')
    ax2.axes.xaxis.set_ticklabels([])
    # let seaborn do it's thing
    #ax = sns.heatmap(R1_Q_mean_by_base, ax=ax, cmap="RdYlGn")
    #sns.heatmap(R1_Q_mean_by_base)
    plt.savefig('Qscore_heatmap', bbox_inches='tight', format='png', dpi=300)
else:
    sns.set_style("whitegrid")
    grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
    fig, ax1 = plt.subplots(1, 1, sharey=True)#,gridspec_kw=grid_kws)
    ax1 = sns.heatmap(R1_Q_mean_by_base, ax=ax1, cmap="RdYlGn",cbar=True)#, vmin=0)
    ax1.set_xlabel('Position (read 1)')
    ax1.set_ylabel('Samples')
    names = unique_sizes[::-1]
    ax1.set_yticklabels(names, fontsize=6, rotation="horizontal")
    ax1.axes.xaxis.set_ticklabels([])
    #ax1.axes.yaxis.set_ticklabels([])
    # let seaborn do it's thing
    #ax = sns.heatmap(R1_Q_mean_by_base, ax=ax, cmap="RdYlGn")
    #sns.heatmap(R1_Q_mean_by_base)
    plt.savefig('Qscore_heatmap_grouped', bbox_inches='tight', format='png', dpi=300)

save_name = (instrument + "_grouped_Q_score_data.txt")
save_file = open(save_name, "w")
newtab = '\t'
newline = '\n'
save_file.write(instrument)
save_file.write(newtab)
for i in unique_sizes:
    save_file.write(i)
    save_file.write(newtab)
save_file.write(newline)
save_file.write("Mean Q score")
save_file.write(newtab)
for i in R1_Q_mean_overall:
    save_file.write(str(i))
    save_file.write(newtab)
save_file.write(newline)
save_file.write("Standard Deviation Q score")
save_file.write(newtab)
for i in R1_Q_stdev_overall:
    save_file.write(str(i))
    save_file.write(newtab)
save_file.write(newline)
save_file.write("Read count")
save_file.write(newtab)
for i in R1_n_reads:
    save_file.write(str(i))
    save_file.write(newtab)
save_file.close()
