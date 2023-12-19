import sys
import numpy as np

name = sys.argv[1]
file_name = sys.argv[2]
output_name = sys.argv[3]

# find the global peak
input_name = "input/" + name
with open(input_name, 'r') as f:
	input_data = f.read().split("\n")[1:]
global_max = ""
global_max_fit = -100
num_seqs = 0
for item in input_data:
	if item!="":
		num_seqs += 1
		seq = item.split()[0]
		fit = float(item.split()[1])
		if fit > global_max_fit:
			global_max = seq
			global_max_fit = fit

data = np.genfromtxt(file_name, delimiter="\t") 

mean_fit = data[:,3]
mean_steps = data[:,4]
peaks = data[:,5]
print(peaks)
# entropy of the distribution of reached peaks
def entropy(peaks_vec, num_seqs):
	probs = [int(x.split(":")[1])/num_seqs for x in peaks_vec.split(",")[:-1]]
	ent = 0
	for p in probs:
		ent += p*log(p)
	return -1*ent
entropies = [entropy(x, num_seqs) for x in peaks]

# probability of reaching the global peak
def prob_global(peaks_vec, global_max, num_seqs):
	peaks = peaks_vec.split(",")[:-1]
	for p in peaks:
		seq = p.split(":")[0]
		if seq==global_max:
			return p.split(":")[1]/num_seqs
probs_global = [prob_global(x, global_max, num_seqs) for x in peaks]

corrs = [np.corrcoef(data[:,2], mean_fit)[0,1],
	np.corrcoef(data[:,2], mean_steps)[0,1],
	np.corrcoef(data[:,2], entropies)[0,1],
	np.corrcoef(data[:,2], probs_global)[0,1]]

with open(output_name, 'a') as f:
	f.write('\t'.join([str(x) for x in [name] + corrs]) + "\n")