'''
Creates a Ostrov code according to specification and saves it to a file.
Parameters:
	[1]-[4] Meaning of free codon blocks 1-4 (1-letter amino acid code, or "-" if no change relative to the SGC)
		Block 1: UUA, UUG
		Block 2: UAG
		Block 3: AGU, AGC
		Block 4: AGG, AGA
	[5] Name of the output file
	[6] Whether special characters should be used (see below). 0 (no) or 1 (yes).
	[7]-[10] (only if the 6th parameter is 1) Amino acid meaning of the special characters X, Z, B, J (in this order); used if there are split codon blocks.
		(each part of the split codon block is assigned a different character, even though it is the same amino acid in fact)
'''

import sys

# Saves a given genetic code to a file
def saveToFile(code, fileName):
	with open(fileName, 'w') as f:
		f.write("Letter\tCodon\n")
		for codon,aa in code.items():
			f.write(aa+"\t"+codon+"\n")



if __name__ == "__main__":
	# parameters
	aa_blocks = sys.argv[1:5]
	fileName = sys.argv[5]
	withSpecChars = sys.argv[6]
	special_chars = ["-"]*4
	if withSpecChars=="0":
		special_chars = sys.argv[7:11]
	
	# free codon blocks
	blocks = {0: ["UUA", "UUG"],
		1: ["UAG"],
		2: ["AGU", "AGC"],
		3: ["AGG", "AGA"]}

	# the standard code
	standard_code = {}
	file = "input/code_standard.tsv"
	with open(file, 'r') as f:
		lines = f.read().split("\n")
		# remove the header
		lines.pop(0)
		for l in lines:
			if l!="":
				splitLine = l.split("\t")
				standard_code[splitLine[1]] = splitLine[0]

	# if special chars are to be used, just substitute the freed codons with the specified character (standard amino acid code or a special character)
	if withSpecChars=="1":
		for i in range(4):
			aa = aa_blocks[i]
			if aa!="-":
				for codon in blocks[i]:
					standard_code[codon] = aa
	# otherwise, decode the special characters according to parameters 7-10
	else:
		for i in range(4):
			aa = aa_blocks[i]
			if aa!="-":
				if aa=="X":
					aa = special_chars[0]
				elif aa=="Z":
					aa = special_chars[1]
				elif aa=="B":
					aa = special_chars[2]
				elif aa=="J":
					aa = special_chars[3]
				for codon in blocks[i]:
					standard_code[codon] = aa

	# save the genetic code to file
	saveToFile(standard_code, fileName)