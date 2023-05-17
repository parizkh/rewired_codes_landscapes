popSize=$1	# population size

cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../02_ruggedness/aa_permutation/input/code_standard.tsv .
cd ..

for i in $(seq 0 99); do
	start=$(expr $i \* 1000)
	end=$(expr $start + 999)
	sbatch -t 1-00 --wrap="sh inner_loop.sh $start $end $popSize"
done