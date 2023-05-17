'''
Generates a summary table of all Ostrov codes: All combinations of amino acids in the four freed codon blocks.
If split codon blocks are created, special characters (X, Z, B, J) are used to encode the "split" part.
The meaning of these special characters (if used) is also saved in the resulting table.
Parameters:
	[1] Name of the output file.
	[2] Name of the file from which the standard genetic code will be read.
'''

import sys
import numpy as np

# a list of amino acids
aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
	'T', 'V', 'W', 'Y', '*']
# assignment of amino acids to physicochemical groups
# 0 - basic, 1 - proline, 2 - aromatic, 3 - aliphatic, 4 - polar, 5 - glycine, 6 - acidic, 7 - stop
physchem_groups = [3, 4, 6, 6, 2, 5, 0, 3, 0, 3, 4, 4, 1, 4, 0, 4, 4, 3, 2, 2, 7]
# the free codon blocks
num_blocks = 4
blocks = {0: ["UUA", "UUG"],
	1: ["UAG"],
	2: ["AGU", "AGC"],
	3: ["AGG", "AGA"]}



# Generates all codons neighboring (i.e., one nucleotide substituion away) codons of codon c.
def gen_neighbouring_codons(c):
	res = []
	for i in range(3):
		for b in ['A', 'C', 'G', 'U']:
			if c[i]!=b:
				res.append(c[:i] + b + c[(i+1):])
	return res

# Checks if a given amino acid (aa) is encoded by any neighbors of a given set of codons (codons), excluding those neighbors in the forbidden_codons list.
def isAmongNeighbours(codons, aa, forbidden_codons, code):
	neighbours = []
	for codon in codons:
		neighbours += gen_neighbouring_codons(codon)
	# remove repeated elements and codons from the codons array from neighbours
	neighbours = set(neighbours) - set(codons) - set(forbidden_codons)
	for n in neighbours:
		if code[n]==aa:
			return True
	return False

# Checks if there is at least one pair of neighboring codons (i.e., in Hamming distance 1) in lists codons1 and codons2.
def areNeighbouring(codons1, codons2):
	for c1 in codons1:
		neighs = gen_neighbouring_codons(c1)
		for c2 in codons2:
			if c2 in neighs:
				return True
	return False

# Generates all possible codons.
def gen_all_codons():
	res = []
	bases = ['A', 'C', 'G', 'U']
	for c1 in bases:
		for c2 in bases:
			for c3 in bases:
				res = res+[c1+c2+c3]
	return res

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


################ MAIN
# parameters
outFile = sys.argv[1]
SGCFile = sys.argv[2]

# the standard code
standard_code = {}
with open(SGCFile, 'r') as f:
	lines = f.read().split("\n")
	# remove the header
	lines.pop(0)
	for l in lines:
		if l!="":
			splitLine = l.split("\t")
			aa = splitLine[0]
			if aa=="X":
				aa="S"
			standard_code[splitLine[1]] = aa



with open(outFile, 'w') as f:
	# write the header
	f.write('\t'.join(["Code", "Block1", "Block2", "Block3", "Block4", "X", "Z", "B", "J", "Robustness"])+"\n")
	
	####### standard
	rob = robustness(standard_code, physchem_groups)
	f.write('\t'.join(["standard", "-", "-", "X", "-", "S", "-", "-", "-", str(rob)])+"\n")

	####### 1 change
	code_num = 0
	# for all blocks
	for i in range(num_blocks):
		block_i = blocks[i]
		# for all amino acids
		for aa_i in aas:
			# this is a change relative to the SGC
			if standard_code[block_i[0]]!=aa_i:
				code = standard_code.copy()

				# reassignments of the codon blocks 
				blocks_reassignments = ["-"]*4
				blocks_reassignments[2] = "X"	# the split serine block
				blocks_reassignments[i] = "Z"

				# meaning of the special characters
				special_chars = {"X": "S", "Z": aa_i, "B": "-", "J": "-"}

				# if i==2, we don't have "X" (the split serine codon block is reassigned)
				if i==2:
					special_chars["X"] = "-"

				# check whether we need the special character
				# if i==3 and aa_i=="S", we can connect block 3 to block 2, both will encode "X"
				if i==3 and aa_i=="S":
					blocks_reassignments[i] = "X"
					special_chars["Z"] = "-"
				elif isAmongNeighbours(blocks[i], aa_i, [], code) or aa_i=="*":
					blocks_reassignments[i] = aa_i
					special_chars["Z"] = "-"

				# compute the robustness of the code
				for codon in blocks[i]:
					code[codon] = aa_i
				rob = robustness(code, physchem_groups)

				f.write("\t".join(["1-"+str(code_num)] + blocks_reassignments + 
					[special_chars["X"], special_chars["Z"], special_chars["B"], special_chars["J"], str(rob)]) + "\n")
				code_num += 1


	####### 2 changes
	code_num = 0
	for i in range(num_blocks):
		for j in range(i+1, num_blocks):
			block_i = blocks[i]
			block_j = blocks[j]

			for aa_i in aas:
				for aa_j in aas:
					# not SGC
					if standard_code[block_i[0]]!=aa_i and standard_code[block_j[0]]!=aa_j:
						code = standard_code.copy()

						# reassignments of the codon blocks
						blocks_reassignments = ["-"]*4
						blocks_reassignments[2] = "X"	# split serine
						blocks_reassignments[i] = "Z"
						blocks_reassignments[j] = "B"

						# meaning of the special characters
						special_chars = {"X": "S", "Z": aa_i, "B": aa_j, "J": "-"}

						### block i
						# we check whether the i-th block has a codon encoding aa_i among its neighbors
						# we need to exclude the j-th block, which is also re-assigned 
						#	(if i==0 and j==1, or i==2 and j==3, the i-th and j-th block are neighbors)
						# if yes, then we do not need a new special character
						if isAmongNeighbours(blocks[i], aa_i, blocks[j], code) or aa_i=="*":
							blocks_reassignments[i] = aa_i
							special_chars["Z"] = "-"
							for codon in blocks[i]:
								code[codon] = aa_i
						else:
							for codon in blocks[i]:
								code[codon] = "Z"
						# if i==2, we don't have the serine split codon block
						if i==2:
							special_chars["X"] = "-"


						### block j
						# Special case: block 2 not reassigned (-> X) and block 3 reassigned to S -> assign to X
						if i!=2 and j==3 and aa_j=="S":
							blocks_reassignments[j] = "X"
							special_chars["B"] = "-"
						# again check whether we need the special character
						elif isAmongNeighbours(blocks[j], aa_j, [], code) or aa_j=="*":
							blocks_reassignments[j] = aa_j
							special_chars["B"] = "-"
							# check if we now need to reassign block i
							if aa_i==aa_j and ((i==0 and j==1) or (i==2 and j==3)):
								blocks_reassignments[i] = aa_i
								special_chars["Z"] = "-"
						else:
							# it's still possible we can use "Z"
							if aa_i==aa_j and ((i==0 and j==1) or (i==2 and j==3)):
								blocks_reassignments[j] = "Z"
								special_chars["B"] = "-"
						# if j==2, we don't have the split serine codon block
						if j==2:
							special_chars["X"] = "-"

						# compute the robustness of the code
						for codon in blocks[i]:
							code[codon] = aa_i
						for codon in blocks[j]:
							code[codon] = aa_j
						rob = robustness(code, physchem_groups)

						# write to output
						f.write("\t".join(["2-"+str(code_num)] + blocks_reassignments + 
							[special_chars["X"], special_chars["Z"], special_chars["B"], special_chars["J"], str(rob)]) + "\n")
						code_num += 1




	########## 3 changes
	code_num = 0
	for i in range(num_blocks):
		for j in range(i+1, num_blocks):
			for k in range(j+1, num_blocks):
				block_i = blocks[i]
				block_j = blocks[j]
				block_k = blocks[k]

				for aa_i in aas:
					for aa_j in aas:
						for aa_k in aas:
							if standard_code[block_i[0]]!=aa_i and standard_code[block_j[0]]!=aa_j and standard_code[block_k[0]]!=aa_k:
								code = standard_code.copy()

								# reassignments of the codon blocks
								blocks_reassignments = ["-"]*4
								blocks_reassignments[2] = "X"	# split serine
								blocks_reassignments[i] = "Z"
								blocks_reassignments[j] = "B"
								blocks_reassignments[k] = "J"

								# meaning of the special characters
								special_chars = {"X": "S", "Z": aa_i, "B": aa_j, "J": aa_k}

								#### the i-th block (block 0 or block 1)
								# we check whether the i-th block has a codon encoding aa_i among its neighbors
								# we need to exclude the j-th block, which is also re-assigned 
								#	(if i==0 and j==1, the i-th and j-th block are neighbors)
								# if yes, then we do not need a new special character
								if isAmongNeighbours(blocks[i], aa_i, blocks[j], code) or aa_i=="*":
									blocks_reassignments[i] = aa_i
									special_chars["Z"] = "-"
									for codon in blocks[i]:
										code[codon] = aa_i
								else:
									for codon in blocks[i]:
										code[codon] = "Z"

								#### the j-th block (block 1 or block 2)
								if j==1:
									# check whether we need the special character
									if isAmongNeighbours(blocks[j], aa_j, [], code) or aa_j=="*":
										blocks_reassignments[j] = aa_j
										special_chars["B"] = "-"
										# check if we now need to re-assign block 0 (block i)
										#	(block j neighbors aa_j and aa_i==aa_j; blocks i (0) and j (1) are neighbors)
										if aa_i==aa_j and special_chars["Z"]!="-":
											special_chars["Z"] = "-"
											blocks_reassignments[i] = aa_i
									else:
										# it is still possible that blocks i and j can be merged
										#	(block i might be reassigned to Z, and thus not found by the isAmongNeighbours function)
										# then we don't need "B", we can assign both i and j to "Z"
										if aa_i==aa_j:
											special_chars["B"] = "-"
											blocks_reassignments[j] = blocks_reassignments[i]
								else:
									assert j==2
									# we don't have the standard split serine codon block now
									special_chars["X"] = "-"
									# check if aa_j is among neighbors of block j
									# exclude block k, because it neighbors block j and is also reassigned
									if isAmongNeighbours(blocks[j], aa_j, blocks[k], code) or aa_j=="*":
										blocks_reassignments[j] = aa_j
										special_chars["B"] = "-"
										for codon in blocks[j]:
											code[codon] = aa_j
									else:
										for codon in blocks[j]:
											code[codon] = "B"

								#### the k-th block
								if k==2:
									assert j==1
									# we don't have the split serine codon block
									special_chars["X"] = "-"
									if isAmongNeighbours(blocks[k], aa_k, [], code) or aa_k=="*":
										blocks_reassignments[k] = aa_k
										special_chars["J"] = "-"
								else:
									assert k==3
									if j==1 and aa_k=="S":
										# the standard split serine codon block is still there and block k can be merged to it
										#aa_k = "X"
										blocks_reassignments[k] = "X"
										special_chars["J"] = "-"
									else:
										# either j==1 and aa_k!="S", or j==2
										if isAmongNeighbours(blocks[k], aa_k, [], code) or aa_k=="*":
											blocks_reassignments[k] = aa_k
											special_chars["J"] = "-"
											# if j==2, we might be able to merge blocks j and k
											if j==2 and aa_j==aa_k and special_chars["B"]!="-":
												special_chars["B"] = "-"
												blocks_reassignments[j] = aa_j
										else:
											# if j==2, blocks j and k can still be merged; then we don;t need "J" and block k is assigned to "B"
											if j==2 and aa_j==aa_k:
												special_chars["J"] = "-"
												blocks_reassignments[k] = blocks_reassignments[j]

								# compute the robustness of the code
								for codon in blocks[i]:
									code[codon] = aa_i
								for codon in blocks[j]:
									code[codon] = aa_j
								for codon in blocks[k]:
									code[codon] = aa_k;
								rob = robustness(code, physchem_groups)

								# write to output
								f.write("\t".join(["3-"+str(code_num)] + blocks_reassignments + 
									[special_chars["X"], special_chars["Z"], special_chars["B"], special_chars["J"], str(rob)]) + "\n")
								code_num += 1



	########### 4 changes
	code_num = 0

	for aa_0 in aas:
		for aa_1 in aas:
			for aa_2 in aas:
				for aa_3 in aas:
					if (	standard_code[blocks[0][0]]!=aa_0 and 
							standard_code[blocks[1][0]]!=aa_1 and 
							standard_code[blocks[2][0]]!=aa_2 and
							standard_code[blocks[3][0]]!=aa_3):

						code = standard_code.copy()

						blocks_reassignments = ["X", "Z", "B", "J"]

						special_chars = {"X": "-", "Z": "-", "B": "-", "J": "-"}

						# block 0
						if isAmongNeighbours(blocks[0], aa_0, blocks[1], code) or aa_0=="*":
							blocks_reassignments[0] = aa_0
							for codon in blocks[0]:
								code[codon] = aa_0
						else:
							special_chars["X"] = aa_0
							for codon in blocks[0]:
								code[codon] = "X"

						# block 1
						if isAmongNeighbours(blocks[1], aa_1, [], code) or aa_1=="*":
							blocks_reassignments[1] = aa_1
							# check if we now need to re-assign block 0
							if aa_0==aa_1 and special_chars["X"]!="-":
								special_chars["X"] = "-"
								blocks_reassignments[0] = aa_0
						else:
							special_chars["Z"] = aa_1
							# check if blocks 0 and 1 can be merged
							if aa_0==aa_1:
								special_chars["Z"] = "-"
								blocks_reassignments[1] = blocks_reassignments[0]

						# block 2
						if isAmongNeighbours(blocks[2], aa_2, blocks[3], code) or aa_2=="*":
							blocks_reassignments[2] = aa_2
							for codon in blocks[2]:
								code[codon] = aa_2
						else:
							special_chars["B"] = aa_2
							for codon in blocks[2]:
								code[codon] = "B"

						# block 3
						if isAmongNeighbours(blocks[3], aa_3, [], code) or aa_3=="*":
							blocks_reassignments[3] = aa_3
							# check if we now need to re-assign block 2
							if aa_2==aa_3 and special_chars["B"]!="-":
								special_chars["B"] = "-"
								blocks_reassignments[2] = aa_2
						else:
							special_chars["J"] = aa_3
							# check if blocks 2 and 3 can be merged
							if aa_2==aa_3:
								special_chars["J"] = "-"
								blocks_reassignments[3] = blocks_reassignments[2]

						# compute the robustness of the code
						for codon in blocks[0]:
							code[codon] = aa_0
						for codon in blocks[1]:
							code[codon] = aa_1
						for codon in blocks[2]:
							code[codon] = aa_2;
						for codon in blocks[3]:
							code[codon] = aa_3;
						rob = robustness(code, physchem_groups)

						# write to output
						f.write("\t".join(["4-"+str(code_num)] + blocks_reassignments + 
							[special_chars["X"], special_chars["Z"], special_chars["B"], special_chars["J"], str(rob)]) + "\n")
						code_num += 1
