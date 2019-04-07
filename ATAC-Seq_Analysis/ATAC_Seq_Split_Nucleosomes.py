#START HERE
###############################
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

folder_name = "<FOLDER_CONTAINING_SAM_FILES>"
data_file_names = os.listdir(folder_name)
plt.switch_backend('agg')

files = []
for i in data_file_names:
    if i[-3:] == "sam":
        files.append((folder_name + "/" + i))

for i in files:
	fname = i
	save_name1 = (fname[:-4] + "_NO_NUC.sam")
	save_file1 = open(save_name1, "w")
	save_name2 = (fname[:-4] + "_MONO_NUC.sam")
	save_file2 = open(save_name2, "w")
	save_name3 = (fname[:-4] + "_DI_NUC.sam")
	save_file3 = open(save_name3, "w")
	save_name4 = (fname[:-4] + "_TRI_NUC.sam")
	save_file4 = open(save_name4, "w")
	save_name5 = (fname[:-4] + "_ALL_NUC.sam")
	save_file5 = open(save_name5, "w")
	newtab = '\t'
	newline = '\n'
	with open(fname) as input:
		size_dist = []
		for line in input:
			if line[0] == "@": #header
				save_file1.write(line)
				save_file2.write(line)
				save_file3.write(line)
				save_file4.write(line)
				save_file5.write(line) 	
			else:
				xx = line.split('\t')
				#xx[0] - read name
				#xx[2] - ref genome
				#xx[3] - map position
				#xx[8] - size 
				#xx[9] - sequence
				#xx[10] - q-score
				size_dist.append(abs(int(xx[8])))
				if abs(int(xx[8])) < 100:
					save_file1.write(line)
				elif 179 < abs(int(xx[8])) < 248:
					save_file2.write(line)
					save_file5.write(line)
				elif 314 < abs(int(xx[8])) < 474:
					save_file3.write(line)
					save_file5.write(line)
				elif 557 < abs(int(xx[8])) < 616:
					save_file4.write(line)
					save_file5.write(line)
	#plot (size histogram)
		fig = plt.figure()
		plt.title("Insert size distribution")
		ax = fig.add_subplot(111)
		ax.spines['top'].set_color('none')
		ax.spines['bottom'].set_color('none')
		ax.spines['left'].set_color('none')
		ax.spines['right'].set_color('none')
		#ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
		ax1 = fig.add_subplot(111)
		x = size_dist
		numBins = 2000
		ax1.hist(x,numBins,color='green')
		fig_save_name = fname[:-4] + "_size_histogram.png"
		plt.savefig(fig_save_name)

		#plot (size histogram log scale)
		fig = plt.figure()
		plt.title("Insert size distribution")
		ax = fig.add_subplot(111)
		ax.spines['top'].set_color('none')
		ax.spines['bottom'].set_color('none')
		ax.spines['left'].set_color('none')
		ax.spines['right'].set_color('none')
		#ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
		ax1 = fig.add_subplot(111)
		x = size_dist
		numBins = 2000
		ax1.hist(x,numBins,color='green')
		plt.yscale('log', nonposy='clip')
		fig_save_name = fname[:-4] + "_size_histogram_log.png"
		plt.savefig(fig_save_name)
		plt.clf()
		save_file1.close()
		save_file2.close()
		save_file3.close()
		save_file4.close()
		save_file5.close()

		#Make unique size list
		unique_bc = set(size_dist)
		#Make dictionary
		BC_dict = dict.fromkeys(unique_bc,0)

		for k in size_dist:
			BC_dict[k] += 1

		#Write out size count file
		save_name = (fname[:-6] + "_fragment_size_counts.txt")
		save_file = open(save_name, "w")
		header = ("Size", '\t',"Counts",'\n')
		save_file.write(''.join(map(str, header)))
		newtab = '\t'
		newline = '\n'
		ordered_keys = sorted(BC_dict)
		ordered_values = []
		for j in ordered_keys:
		    save_file.write(str(j))
		    save_file.write(newtab)
		    ordered_values.append(BC_dict[j])
		    save_file.write(str(BC_dict[j]))
		    save_file.write(newline)
		save_file.close()
		
#Convert sam to bam and index bam files
data_file_names = os.listdir(folder_name)
files = []
for i in data_file_names:
    if i[-7:] == "NUC.sam":
        files.append((folder_name + "/" + i))

for i in files:
	#Convert sam to bam file
	bam_name = i[:-3] + "bam"
	execute = "samtools view -Sb " + i + " > " + bam_name
	os.system(execute)
	#Sort bam file
	bam_sort_name = i[:-4] + "_sorted.bam"
	execute = "samtools sort " + bam_name + " -T Temp -o " + bam_sort_name
	os.system(execute)
	#Index bam file
	execute = "samtools index " + bam_sort_name
	os.system(execute)


