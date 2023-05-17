cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
cd ..

echo -n > output/results

for code in "FS20" "RED20" "OPT" "OPT-NR" "CMC" "CMC2" "REC" "Ostrov"; do
	codeFile="input/codes/code_"$code".tsv"
	../../../code/./greedy_walk input/map.tsv $codeFile tmpOut
	echo -ne "$code\t" >> output/results
	cat tmpOut >> output/results
done

rm tmpOut
