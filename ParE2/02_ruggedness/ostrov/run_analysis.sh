cd input
ln -s ../../aa_permutation/input/code_standard.tsv .
ln -s ../../aa_permutation/input/map.tsv .
ln -s ../../../../GB1/02_ruggedness/ostrov/input/ostrov_codes_summary.tsv .
cd ..

for i in $(seq 0 3); do
	start=$(expr $i \* 50000 + 2)
	end=$(expr $start + 49999)
	sbatch -t 1-00 --wrap="sh inner_loop.sh $start $end"
done
