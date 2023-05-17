'''
Generates a randomized genetic code (using amino acid permutaiton or random codon assignment), using a specified random seed,
computes its robustness and saves it to a file.
Parameters:
	[1] randomization type (aa_permutation or random)
	[2] seed 
	[3] output file; the value of the code robustness will be saved here
	[4] code output file; the genetic code will be saved here
'''

import numpy as np
import math
import sys

# a list of amino acids
aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
	'T', 'V', 'W', 'Y', '*']

# assignment of amino acids to physicochemical groups, according to Fig. 1A in Pines et al., mBio, 2017
# 0 - basic, 1 - proline, 2 - aromatic, 3 - aliphatic, 4 - polar, 5 - glycine, 6 - acidic, 7 - stop
physchem_groups = [3, 4, 6, 6, 2, 5, 0, 3, 0, 3, 4, 4, 1, 4, 0, 4, 4, 3, 2, 2, 7]


# Generates all codons.
def gen_all_codons():
	res = []
	bases = ['A', 'C', 'G', 'U']
	for c1 in bases:
		for c2 in bases:
			for c3 in bases:
				res = res+[c1+c2+c3]
	return res

# Generates all codons neighboring (i.e., one nucleotide substituion away) codons of codon c.
def gen_neighbouring_codons(c):
	res = []
	for i in range(3):
		for b in ['A', 'C', 'G', 'U']:
			if c[i]!=b:
				res.append(c[:i] + b + c[(i+1):])
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


'''
Computes the robustness of the code.
Robustness is defined as the proportion of single-nucleotide substitutions that do not change the physicochemical 
properties of amino acids.
'''
def robustness(code, groups = None):
	# number of distinct groups
	if groups != None:
		num_groups = len(np.unique(groups))
	# if not specified, equivalent to number of amino acids
	else:
		groups = range(21)
		num_groups = 21
	# all codons
	codons = gen_all_codons()
	# all substitutions
	num_subs = 0
	# conservative substitutions
	num_con = 0

	for i in range(len(codons)):
		codon = codons[i]
		aa = code[codon]
		group = groups[aas.index(aa)]
		neighbours = gen_neighbouring_codons(codon)
		for n in neighbours:
			num_subs += 1
			if groups[aas.index(code[n])] == group:
				num_con += 1

	return num_con/num_subs



#####################
# parameters
randType = sys.argv[1]
seed = int(sys.argv[2])
outFile = sys.argv[3]
of = open(outFile, 'a')
codeFile = sys.argv[4]


# the standard code
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


# generate the randomized code
if seed==0:
	code = standard_code
elif randType=="aa_permutation":
	code = rand_aa_permutation(standard_code, seed)
elif randType=="random":
	code = rand_randomAssignment(standard_code, seed)

# compute robustnes of the code
rob = robustness(code, physchem_groups)	
# for aa permutation: the amino acid occupying the split codon block
if randType=="aa_permutation":
	X = code["UCU"]		

# save the code to file
with open(codeFile, 'w') as f:
	f.write("Aa\tCodon\n")
	for codon,aa in code.items():
		f.write(aa+"\t"+codon+"\n")

# write the statistics to outputFile
if randType=="aa_permutation":
	of.write('\t'.join([str(x) for x in [seed, X, rob]]) + "\t")
else:
	of.write('\t'.join([str(x) for x in [seed, rob]]) + "\t")


of.close()