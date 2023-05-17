cd input
ln -s ../../output/results/results
ln -s ../../input/code_standard.tsv .
ln -s ../../../../../GB1/02_ruggedness/aa_permutation/02_03_aaindex/input/aaindex1 .
ln -s ../../../../../GB1/02_ruggedness/aa_permutation/02_03_aaindex/input/Bartonek_PNAS_2020_S01.csv .
# sort the results numerically by random seed and only keep the relevant columns
sort -g results | cut -f 4,7-10 > results_sorted.tsv	
cd ..

sbatch -t 1-00 --wrap="python3 ../../../../code/correlations_aaindex.py input/results_sorted.tsv output/results_aaindex 4 D,W,E"
sbatch -t 1-00 --wrap="python3 ../../../../code/correlations_null_aaindex.py input/results_sorted.tsv output/results_aaindex_null 4 D,W,E"
