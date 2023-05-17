cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../02_ruggedness/aa_permutation/input/code_standard.tsv .
ln -s ../../../02_ruggedness/ostrov/input/ostrov_codes_summary.tsv .
cd ..

for i in $(seq 0 4); do
	start=$(expr $i \* 40000 + 2)
	end=$(expr $start + 39999)
	sbatch -t 1-00 --wrap="sh inner_loop.sh $start $end input/map.tsv"
done
