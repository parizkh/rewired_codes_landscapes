startSeed=$1
endSeed=$2

# output file
outFile="output/results/results_"$startSeed"-"$endSeed
mkdir -p "output/results"
# temporary files
codeFile="tmp_code_"$startSeed"-"$endSeed
tmpOut="tmp_out_"$startSeed"-"$endSeed

for i in $(seq $startSeed $endSeed); do
	# generate the genetic code
	python3 ../../../code/generate_gen_code.py aa_permutation $i $outFile $codeFile
	# run the ruggedness analysis
	../../../code/./landscape_ruggedness input/map.tsv $codeFile $tmpOut
	# append the results to the output file
	cat $tmpOut >> $outFile
done

# remove the temporary files
rm $codeFile
rm $tmpOut
