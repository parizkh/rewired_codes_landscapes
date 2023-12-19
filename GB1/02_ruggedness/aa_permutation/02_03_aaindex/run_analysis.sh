mkdir input
cd input
ln -s ../../output/results/results
ln -s ../../../../../resources/code_standard.tsv .
# sort the results numerically by random seed and only keep the relevant columns
sort -g results | cut -f 4,7-10 > results_sorted.tsv	
sort -g results > results_ruggedness
ln -s ../../02_01_global_peak_size/results_globalPeakSize_aaPerm results_globalPeakSize
ln -s ../../../../../resources/aaindex1 .
ln -s ../../input/map.tsv .
cd ..
mkdir output

sbatch -t 1-00 --wrap="python3 ../../../../code/correlations_aaindex.py input/results_sorted.tsv output/results_aaindex 4 W,L,A"
sbatch -t 5-00 --wrap="python3 ../../../../code/correlations_null_aaindex.py input/results_sorted.tsv output/results_aaindex_null 4 W,L,A"
