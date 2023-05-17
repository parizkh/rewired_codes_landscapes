cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../02_ruggedness/aa_permutation/input/code_standard.tsv .
cd ..

for i in $(seq 0 3); do
	start=$(expr $i \* 25000)
	end=$(expr $start + 24999)
	sbatch -t 1-00 --wrap="sh inner_loop.sh $start $end"
done
