cd input
ln -s ../../aa_permutation/input/code_standard.tsv .
ln -s ../../aa_permutation/input/map.tsv .
cd ..

# generate the summary of the Ostrov codes
python3 ../../../code/generate_codes_ostrov.py input/ostrov_codes_summary.tsv input/code_standard.tsv
sed -i 's/\*/O/g' input/ostrov_codes_summary.tsv		# encode * (stop codon) as O

for i in $(seq 0 129); do
	start=$(expr $i \* 1500 + 2)
	end=$(expr $start + 1499)
	sbatch -t 5-00 --mem-per-cpu=5000 --wrap="sh inner_loop.sh $start $end"
done