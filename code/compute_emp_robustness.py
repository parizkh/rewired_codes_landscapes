'''
Computes the correlation of landscape ruggedness results / simulation results with robustness of the genetic codes,
in terms of all properties from the aaindex database.
Parameters:
	[1] File with the landscape ruggedness / simulation results
	[2] Name of the output file
	[3] Which columns contain results that depend on the global peak size; for these columns, only codes that preserve the size of the global peak and under which the 
		global peak forms a single connected region in the genotype space will be considered. Numbers (0-based), separated by commas (no spaces).
	[4] Which amino acids constitute the global peak (so that codes where these amino acids occupy the split codon block can be discarded).
'''

import numpy as np
import sys

'''
Reads the genetic code from a file.
Format of the file: 
	First line: header
	Each other line: Amino acid \t Codon
'''
def read_code(fileName):
	code = {}
	with open(fileName, 'r') as f:
		lines = f.read().split("\n")
		# remove the header
		lines.pop(0)
		for l in lines:
			if l!="":
				splitLine = l.split("\t")
				#print(splitLine)
				code[splitLine[1]] = splitLine[0]
	return code

# Generates all codons neighboring (i.e., one nucleotide substituion away) codons of codon c.
def gen_neighbouring_codons(c):
	res = []
	for i in range(3):
		for b in ['A', 'C', 'G', 'U']:
			if c[i]!=b:
				res.append(c[:i] + b + c[(i+1):])
	return res

# Computes the number of single-nucleotide mutations leading between amino acid pairs in a given genetic code.
def number_of_mutations(code, aa_pairs):
	res = np.zeros(shape=len(aa_pairs))

	codons = list(code.keys())

	for c in codons:
		aa = code[c]
		if aa!="*":
			neighs = gen_neighbouring_codons(c)
			for n in neighs:
				n_aa = code[n]
				if n_aa!="*" and aa!=n_aa:
					try:
						ind = aa_pairs.index(aa+n_aa)
					except ValueError:
						ind = aa_pairs.index(n_aa+aa)
					res[ind]+=1
	return res

# Computes the absolute differences between all pairs of values in vector v.
def compute_diff_vector(v):
	res = []
	for ind1 in range(len(v)):
		for ind2 in range(ind1+1, len(v)):
			res += [abs(v[ind1]-v[ind2])]
	return res

''' 
amino acid permutation: 
Keep the blocks of codons, only randomly re-assign amino acids to them.
Position of stop codons is fixed.
parameters: 
	standard_code ... dictionary providing the standard genetic code; keys = codons, values = amino acids
	seed ... random seed (default 0)
'''
def rand_aa_permutation(standard_code, seed = 0):
	np.random.seed(seed)
	# a list of amino acids
	aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
	'T', 'V', 'W', 'Y', '*']
	# dictionary to store the code
	code = {}
	# randomly shuffle the amino acids 
	shuffled_aas = np.random.permutation(20)
	# the last one is stop codon
	shuffled_aas = np.append(shuffled_aas, [20])
	# go over the standard code and always change amino acid i to shuffled_aas[i]
	for codon, aa in standard_code.items():
		code[codon] = aas[shuffled_aas[aas.index(aa)]]
	return code

'''
Amino acid permutation, restricted to only those permutation that preserve the number of codons per aa.
parameters: 
	standard_code ... dictionary providing the standard genetic code; keys = codons, values = amino acids
	seed ... random seed (default 0)
'''
def rand_aa_permutation_restricted(standard_code, seed = 0):
	np.random.seed(seed)
	# dictionary to store the code
	code = {}
	# aas, ordered by number of codons
	aas = ["M", "W", "C", "D", "E", "F", "H", "K", "N", "Q", "Y", "I", "A", "G", "P", "T", "V", "L", "R", "S", "*"]
	# randomly shuffle the amino acids, preserving the size
	shuffled_aas_1 = np.random.permutation(2)
	shuffled_aas_2 = [x+2 for x in np.random.permutation(9)]
	shuffled_aas_3 = [11]
	shuffled_aas_4 = [x+12 for x in np.random.permutation(5)]
	shuffled_aas_6 = [x+17 for x in np.random.permutation(3)]
	shuffled_aas = np.concatenate((shuffled_aas_1, shuffled_aas_2, shuffled_aas_3, shuffled_aas_4, shuffled_aas_6, [20]))
	# the last one is stop codon
	shuffled_aas = np.append(shuffled_aas, [20])
	# go over the standard code and always change amino acid i to shuffled_aas[i]
	for codon, aa in standard_code.items():
		code[codon] = aas[shuffled_aas[aas.index(aa)]]
	return code



def gen_all_neighbors(seq):
	aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
	'T', 'V', 'W', 'Y']
	res = []
	for i in range(len(seq)):
		for aa in aas:
			new_seq = seq[:i] + aa + seq[(i+1):]
			if new_seq != seq:
				res.append(new_seq)
	return res

############## MAIN 

if __name__ == "__main__":
	# number of randomized genetic codes
	n_codes = 100000

	
	# the amino acid pairs
	aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
	aa_pairs = []
	for i1 in range(len(aas)):
		for i2 in range(i1+1, len(aas)):
			aa_pairs += [aas[i1]+aas[i2]]

	# compute the empirical exchangeabilities (mean absolute fitness difference when changing X->Y, in different backgrounds)
	dataFile = "../input/map.tsv"
	data_dict = {}
	with open(dataFile, 'r') as f:
		data = f.read().split("\n")[1:]
	for item in data:
		if item!="":
			seq = item.split()[0]
			fitness = float(item.split()[1])
			data_dict[seq] = fitness
	empirical_diffs = [0]*len(aa_pairs)
	num_contexts = [0]*len(aa_pairs)
	for seq in data_dict.keys():
		for i in range(len(seq)):
			for aa in aas:
				new_seq = seq[:i] + aa + seq[(i+1):]
				if new_seq != seq:
					aa_pair = seq[i]+new_seq[i] if seq[i]+new_seq[i] in aa_pairs else new_seq[i]+seq[i]
					fitness_diff = abs(data_dict[seq] - data_dict[new_seq])
					ind = aa_pairs.index(aa_pair)
					empirical_diffs[ind] += fitness_diff
					num_contexts[ind] += 1
	
	empirical_diffs = [empirical_diffs[i]/num_contexts[i] for i in range(len(aa_pairs))]
	with open("output/empirical_diffs.tsv", 'w') as f:
		for i in range(len(aa_pairs)):
			f.write(aa_pairs[i] + " \t" + str(empirical_diffs[i]) + '\n')
	
	
	# the standard genetic code
	standard_code = read_code("input/code_standard.tsv")

	# create the matrix with numbers of aa-aa mutations allowed by each code
	codes_mutations = np.zeros(shape=(n_codes, len(aa_pairs)))
	
	for seed in range(n_codes):		
		if seed==0:
			# standard			
			code = standard_code
		else:
			code = rand_aa_permutation_restricted(standard_code, seed)
		# compute the number of allowed amino acid substitutions
		codes_mutations[seed,:] = number_of_mutations(code, aa_pairs)
			

	print("codes done")
	np.savetxt("output/codes_muts_matrix.csv", codes_mutations, fmt="%i", delimiter=",", header=",".join(aa_pairs), comments="")

	# empirical diffs
	ermc = np.matmul(codes_mutations, empirical_diffs)
	with open("output/ermc_emp", 'w') as f:
		for i in range(100000):
			f.write('\t'.join([str(i), str(ermc[i])]) + "\n")
