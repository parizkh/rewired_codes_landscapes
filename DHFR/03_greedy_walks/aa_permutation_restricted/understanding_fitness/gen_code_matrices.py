import numpy as np

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


if __name__ == "__main__":
	# the amino acid pairs
	aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	aa_pairs = []
	for i1 in range(len(aas)):
		for i2 in range(i1+1, len(aas)):
			aa_pairs += [aas[i1]+aas[i2]]

	# create the matrix with numbers of aa-aa mutations allowed by each code
	# high-fitness codes
	f = open("codes_highFitness", 'r')
	ids_high = f.read().split()
	f.close()
	codes_mutations_high = np.zeros(shape=(100, len(aa_pairs)))	
	ind=0
	for i in ids_high:		
		code = read_code("codes_high/code_"+i)
		# compute the number of allowed amino acid substitutions
		codes_mutations_high[ind,:] = number_of_mutations(code, aa_pairs)
		ind += 1

	np.savetxt("matrix_high", codes_mutations_high, fmt="%i", delimiter=",", header=",".join(aa_pairs), comments="")

	# high-fitness codes
	f = open("codes_lowFitness", 'r')
	ids_low = f.read().split()
	f.close()
	codes_mutations_low = np.zeros(shape=(100, len(aa_pairs)))	
	ind=0
	for i in ids_low:		
		code = read_code("codes_low/code_"+i)
		# compute the number of allowed amino acid substitutions
		codes_mutations_low[ind,:] = number_of_mutations(code, aa_pairs)
		ind += 1

	np.savetxt("matrix_low", codes_mutations_low, fmt="%i", delimiter=",", header=",".join(aa_pairs), comments="")
	