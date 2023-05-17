startSeed=$1
endSeed=$2

outFile="output/results/results_"$startSeed"-"$endSeed
mkdir -p "output/results"
codeFile="tmp_code_"$startSeed"-"$endSeed
tmpOut="tmp_out_"$startSeed"-"$endSeed

for i in $(seq $startSeed $endSeed); do
	python3 ../../../code/generate_gen_code_new.py random $i $outFile $codeFile
	../../../code/./landscape_ruggedness input/map.tsv $codeFile $tmpOut
	cat $tmpOut >> $outFile
done

rm $codeFile
rm $tmpOut
