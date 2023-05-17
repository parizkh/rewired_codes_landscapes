cd input
ln -s ../../output/results/results
ln -s ../../input/code_standard.tsv .
# sort the results numerically by random seed and only keep the relevant columns
sort -g results | cut -f 4,7-10 > results_sorted.tsv	
cd ..

sbatch -t 1-00 --wrap="python3 ../../../../code/correlations_aaindex.py input/results_sorted.tsv output/results_aaindex 4 W,L,A"
sbatch -t 1-00 --wrap="python3 ../../../../code/correlations_null_aaindex.py input/results_sorted.tsv output/results_aaindex_null 4 W,L,A"
