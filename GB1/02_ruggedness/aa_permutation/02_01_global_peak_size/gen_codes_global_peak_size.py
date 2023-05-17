import numpy as np
import math
import sys

# a list of amino acids
aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
	'T', 'V', 'W', 'Y', '*']


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
Random codon assignemnt:
Assign codons to amino acids independently and randomly, no requirements on degeneracy.
The only requirement is that all 20 amino acids have at least one codon assigned to them.
parameters: 
	standard_code ... dictionary with the SGC 
	seed ... random seed (default 0)
'''
def rand_randomAssignment(standard_code, seed = 0):
	np.random.seed(seed)
	# dictionary to store the code
	code = {}
	# list of amino acids
	aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	# fix stop codons
	code['UAA'] = "*"
	code['UAG'] = "*"
	code['UGA'] = "*"
	# shuffle the codons
	# stop codons are fixed - don't include them in the permutation
	codons = list(standard_code.keys())
	perm = np.random.permutation(64)
	codons = [codons[x] for x in perm if codons[x] not in ['UAA', 'UAG', 'UGA']]
	# the first 20 codons are assigned each to the corresponding amino acid
	# to fulfill the requirement that each amino acid has at least one codon
	for i in range(20):
		code[codons[i]] = aas[i]
	# the remaining codons are assigned to amino acids randomly
	for i in range(41):
		codon = codons[20+i]
		aa = aas[np.random.randint(0,20)]
		code[codon] = aa

	return code





#####################
# parameters
randType = sys.argv[1]	# randomization type: "aa_pemutation" or "random"
outFile = sys.argv[2]	# output file
of = open(outFile, 'w')


#################### the standard code
standard_code={}
file = "input/code_standard.tsv"
with open(file, 'r') as f:
	lines = f.read().split("\n")
	# remove the header
	lines.pop(0)
	for l in lines:
		if l!="":
			splitLine = l.split("\t")
			standard_code[splitLine[1]] = splitLine[0]


######## randomized codes
for seed in range(100000):
	if seed==0:
		code = standard_code
	elif randType=="aa_permutation":
		code = rand_aa_permutation(standard_code, seed)
	elif randType=="random":
		code = rand_randomAssignment(standard_code, seed)

	# compute the size of the global peak (WWLA)
	vals = list(code.values())
	global_size = (vals.count("W")**2) * vals.count("L") * vals.count("A")

	# write the statistics to outputFile
	of.write('\t'.join([str(x) for x in [seed, global_size]]) + "\n")


of.close()