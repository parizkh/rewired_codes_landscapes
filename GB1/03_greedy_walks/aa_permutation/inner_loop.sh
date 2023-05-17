startSeed=$1
endSeed=$2

outFile="output/results_"$startSeed"-"$endSeed
mkdir -p "output"
codeFile="tmp_code_"$startSeed"-"$endSeed
tmpOut="tmp_out_"$startSeed"-"$endSeed

for i in $(seq $startSeed $endSeed); do
	python3 ../../../code/generate_gen_code.py aa_permutation $i $outFile $codeFile
	../../../code/./greedy_walk input/map.tsv $codeFile $tmpOut
	cat $tmpOut >> $outFile
done

rm $codeFile
rm $tmpOut