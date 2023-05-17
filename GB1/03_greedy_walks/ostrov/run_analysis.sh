cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../02_ruggedness/aa_permutation/input/code_standard.tsv .
ln -s ../../../02_ruggedness/ostrov/input/ostrov_codes_summary.tsv .
cd ..

for i in $(seq 0 97); do
	start=$(expr $i \* 2000 + 2)
	end=$(expr $start + 1999)
	sbatch -t 5-00 --mem-per-cpu=10000 --wrap="sh inner_loop.sh $start $end input/map.tsv"
done