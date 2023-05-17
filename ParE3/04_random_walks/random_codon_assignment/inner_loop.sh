startSeed=$1
endSeed=$2
popSize=$3

outFile="output/results_"$popSize"/results_"$startSeed"-"$endSeed
mkdir -p "output/results_"$popSize
codeFile="tmp_code_"$popSize"_"$startSeed"-"$endSeed
tmpOut="tmp_out_"$popSize"_"$startSeed"-"$endSeed

for i in $(seq $startSeed $endSeed); do
	python3 ../../...code/generate_gen_code.py random $i $outFile $codeFile
	../../../code/./random_walk input/map.tsv $codeFile $tmpOut $popSize
	cat $tmpOut >> $outFile
done

rm $codeFile
rm $tmpOut