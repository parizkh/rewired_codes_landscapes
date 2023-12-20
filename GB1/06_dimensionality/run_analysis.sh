# prepare the input
mkdir input
cd input
ln -s ../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../resources/code_standard.tsv .
cd ..

# generate the inputs
# L=3: all 80 possibilities
# L=2: sample 100
python3 ../../code/generate_dimReduced_inputs.py

# ruggedness
mkdir output_ruggedness
for input in input/map_3*; do
	sbatch -t 5-00 --wrap="sh inner_loop_ruggedness.sh $input 3"
done
for input in input/map_2*; do
	sbatch -t 5-00 --wrap="sh inner_loop_ruggedness.sh $input 2"
done

# greedy walks
mkdir output_greedy_walks
for input in input/map_3*; do
	sbatch -t 5-00 --wrap="sh inner_loop_greedy.sh $input 3"
done
for input in input/map_2*; do
	sbatch -t 5-00 --wrap="sh inner_loop_greedy.sh $input 2"
done

# random walks
mkdir output_random_walks
for input in input/map_*; do
	for popSize in 10 100 10000 1000000; do
		sbatch -t 5-00 --wrap="sh inner_loop_random.sh $input $popSize"
	done
done
