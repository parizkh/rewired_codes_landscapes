cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../02_ruggedness/aa_permutation/input/code_standard.tsv .
cd ..

for i in $(seq 0 49); do
	start=$(expr $i \* 2000)
	end=$(expr $start + 1999)
	sbatch -t 5-00 --mem-per-cpu=10000 --wrap="sh inner_loop.sh $start $end"
done