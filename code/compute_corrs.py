import sys
import numpy as np

name = sys.argv[1]
file_name = sys.argv[2]
output_name = sys.argv[3]

data = np.genfromtxt(file_name, delimiter="\t") 

corrs = []
for i in range(3, data.shape[1]):
	# if not all the values are same
	if len(set(data[:,i])) > 1:
		# correlation of i-th column with code robustness (3rd column)
		corrs.append(np.corrcoef(data[:,i], data[:,2])[0,1])


with open(output_name, 'a') as f:
	f.write('\t'.join([str(x) for x in [name] + corrs]) + "\n")
