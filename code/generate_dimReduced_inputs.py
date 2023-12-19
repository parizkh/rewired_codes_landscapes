### Script to generate inputs with reduced dimensionality, GB1 data.
### Generates all 80 subsets with L=3, randomly samples 100 with L=2.

import random

aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

with open("input/map.tsv", 'r') as f:
	data = f.read().split("\n")[1:]

data_dict = {}
for item in data:
	if item!="":
		split_item = item.split()
		data_dict[split_item[0]] = split_item[1]

# L=3
for i in range(4):
	for background in aas:
		file_name = "input/map_3_" + str(i) + "_" + background + ".tsv"
		with open(file_name, 'w') as f:
			f.write("sequence\tphenotype\n")
			for seq, val in data_dict.items():
				if seq[i]==background:
					new_seq = seq[:i] + seq[(i+1):]
					f.write(new_seq + "\t" + val + "\n")

# L=2
random.seed(0)
for i in range(100):
	# choose 2 random sites and fix the background
	site1 = random.randint(0,3)
	aa1 = aas[random.randint(0,19)]
	site2=site1		# site2 must be different from site1
	while site2==site1:
		site2 = random.randint(0,3)
	aa2 = aas[random.randint(0,19)]
	sites_var = [x for x in [0,1,2,3] if x not in [site1, site2]]

	file_name = "input/map_2_seed_" + str(i) + ".tsv"
	with open(file_name, 'w') as f:
		f.write("sequence\tphenotype\n")
		for seq, val in data_dict.items():
			if seq[site1]==aa1 and seq[site2]==aa2:
				new_seq = seq[sites_var[0]] + seq[sites_var[1]]
				f.write(new_seq + "\t" + val + "\n")