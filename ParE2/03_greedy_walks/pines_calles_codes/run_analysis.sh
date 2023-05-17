cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../../GB1/03_greedy_walks/pines_calles_codes/input/codes .
cd ..

echo -n > output/results

for code in "FS20" "RED20" "OPT" "OPT-NR" "CMC" "CMC2" "REC" "Ostrov"; do
	codeFile="input/codes/code_"$code".tsv"
	../../../code/./greedy_walk_ParD3 input/map.tsv $codeFile tmpOut
	cat tmpOut >> output/results
done

rm tmpOut
