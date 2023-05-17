import sys
import numpy as np

int_to_nuc = {0: 'A', 1: 'C', 2: 'G', 3: 'U'}

# Generates all codons neighboring (i.e., one nucleotide substituion away) codons of codon c.
def gen_neighbouring_codons(c):
	res = []
	for i in range(3):
		for b in ['A', 'C', 'G', 'U']:
			if c[i]!=b:
				res.append(c[:i] + b + c[(i+1):])
	return res

# Translates a given RNA seq using given genetic code
def translate(seq, code):
	t = ""
	for i in range(int(len(seq)/3)):
		t += code[seq[3*i:(3*i+3)]]
	return t

# Generates a random RNA sequence of a given length
def rand_RNA_seq(length):
	ints = np.random.choice(4, length, replace=True)
	chars = [int_to_nuc[x] for x in ints]
	return ''.join(chars)



################################
inFile = sys.argv[1]		# input file
numSquares = int(sys.argv[2])		# how many squares exhibiting reciprocal-sign epistasis should be added
seed = int(sys.argv[3])		# random seed

np.random.seed(seed)

# the standard genetic code
code = {}
with open("../input/code_standard.tsv", 'r') as f:
	lines = f.read().split("\n")[1:]
for line in lines:
	if line!="":
		aa = line.split()[0]
		codon = line.split()[1]
		code[codon] = aa


# read the input file
data = {}
with open(inFile, 'r') as f:
	lines = f.read().split("\n")[1:]
for line in lines:
	if line!="":
		seq = line.split()[0]
		fit = line.split()[1]
		data[seq] = float(fit)
seqs = list(data.keys())

# generate given number of rs squares
for i in range(numSquares):
	# sample random sequence
	while True:
		seq0 = rand_RNA_seq(12)
		tseq0 = translate(seq0, code)
		# the sequence should not contain stop codons and not be the global peak
		if "*" not in tseq0 and tseq0!="WWLA":
			score0 = data[tseq0]
			break

	while True:
		# first mutation
		while True:
			pos1 = np.random.randint(12)
			new_char1 = int_to_nuc[np.random.randint(4)]
			seq1 = seq0[:pos1] + new_char1 + seq0[(pos1+1):]
			tseq1 = translate(seq1, code)
			# this is not a synonymous mutation, the sequence is not a global peak and it does not contain stop codons
			if tseq1!=tseq0 and tseq1!="WWLA" and "*" not in tseq1:
				score1 = data[tseq1]
				break

		# second mutation
		while True:
			pos2 = np.random.randint(12)
			new_char2 = int_to_nuc[np.random.randint(4)]
			seq2 = seq0[:pos2] + new_char2 + seq0[(pos2+1):]
			tseq2 = translate(seq2, code)
			# in addition to all the condition above, the mutation must happen in a different position than the first one
			if tseq2!=tseq0 and tseq2!=tseq1 and tseq2!="WWLA" and pos2!=pos1 and "*" not in tseq2:
				score2 = data[tseq2]
				break


		# double mutant
		seq12 = seq1[:pos2] + new_char2 + seq1[(pos2+1):]
		tseq12 = translate(seq12, code)
		if tseq12!=tseq0 and tseq12!=tseq1 and tseq12!=tseq2 and "*" not in tseq12:
			score12 = data[tseq12]
			break

	# permute the socres so that the double mutant has the highest score, the wildtype the second highest, and the two single mutants the two lowest
	scores = [score0, score1, score2, score12]
	scores.sort()
	data[tseq0] = scores[2]
	data[tseq1] = scores[0]
	data[tseq2] = scores[1]
	data[tseq12] = scores[3]

	


# print the new input data
with open("tmp_input_rsEpi", 'w') as f:
	f.write("sequence\tphenoytpe\n")
	for seq, score in data.items():
		f.write(seq + "\t" + str(score) + "\n")