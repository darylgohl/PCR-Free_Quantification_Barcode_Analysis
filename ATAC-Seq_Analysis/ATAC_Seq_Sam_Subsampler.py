#START HERE
###############################
import os
import sys
import numpy as np
import random
from shutil import copyfile


folder_name = "<FOLDER_CONTAINING_SAM_FILES>" 
os.chdir(folder_name)
os.mkdir("subsampled")
out_dir = os.path.join(folder_name,"subsampled")
data_file_names = os.listdir(folder_name)
#depth = 20000000
depth = 250000

files = []
for i in data_file_names:
    if i[-10:] == "nodupl.sam":
        files.append((folder_name + "/" + i))

read_counts = []
for i in files:
	count = 0
	with open(i) as input:
		for line in input:
			count += 1
	read_counts.append(((count-23))/2)

#subsample
for i, item in enumerate(files):
	fname = item
	save_name1 = (fname[:-4] + "_subsampled.sam")
	if read_counts[i] < depth:
		dst = os.path.join(out_dir,save_name1.split("/")[-1])
		src = item
		copyfile(src, dst)
	else:
		c = list(range(1, int(read_counts[i]-1)))
		inds = random.sample(c, depth)
		inds.sort()
		inds.reverse() #Note: pop from right side is much more time efficient
		fname = item
		save_name = (fname[:-4] + "_subsampled.sam")
		save_name1 = os.path.join(out_dir,save_name.split("/")[-1])
		save_file1 = open(save_name1, "w")
		newtab = '\t'
		newline = '\n'
		with open(fname) as input:
			count = 0
			for line in input:
				count += 1
				if line[0] == "@": #header
					save_file1.write(line)
					count -= 1	
				elif inds == []:
					break
				elif count == (inds[-1]*2-1):
					save_file1.write(line)
				elif count == (inds[-1]*2):
					inds.pop(-1)
					save_file1.write(line)
			save_file1.close()
