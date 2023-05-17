import sys
import numpy as np

inFile = sys.argv[1]	# the input file
numPeaks = int(sys.argv[2])		# how many local peaks should be adde
seed = int(sys.argv[3])		# random seed

np.random.seed(seed)

# generate the peaks
while True:
	peaks = np.random.choice(194481, numPeaks, replace=False)
	# 151380 is the index of the global peak
	if 151380 not in peaks:
		break


# read the input file
with open(inFile, 'r') as f:
	with open("tmp_input_peaks", 'w') as of:
		lines = f.read().split("\n")
		of.write(lines[0]+"\n")
		for i in range(1,len(lines)):
			if i-1 in peaks:		# this is now an artificial local peak
				of.write(lines[i].split("\t")[0] + "\t2.5\n")
			else:
				of.write(lines[i]+"\n")

