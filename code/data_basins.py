# results, codes that preserve size of the global peak
with open("output/results_constSummitSize.tsv", 'r') as f:
	data = f.read().split("\n")

# save the basin sizes, one basin one code per one line
with open("output/basins.tsv", 'w') as f:
	data.pop(0)
	for line in data:
		if line!="":
			line_split = line.split("\t")
			seed = line_split[0]
			peaks = line_split[5]
			rob = line_split[2]
			peaks_split = peaks.split(",")
			for s_peak in peaks_split:
				if s_peak!="":
					f.write(seed + "\t" + rob + "\t" + s_peak.split(":")[0] + "\t" + s_peak.split(":")[1]+"\n")