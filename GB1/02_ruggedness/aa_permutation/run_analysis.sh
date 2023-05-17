# prepare the input
cd input
ln -s ../../../01_vcregression/output/map.txt .
cd ..
cat input/map.txt | sed 's/,/\t/g' > input/map.tsv

# run the ruggedness analysis
for i in $(seq 0 99); do
	start=$(expr $i \* 1000)
	end=$(expr $start + 999)
	sbatch -t 1-00 --mem-per-cpu=5000 --wrap="sh inner_loop.sh $start $end"
done

# compute the global peak sizes for all codes
python3 02_01_global_peak_size/gen_codes_global_peak_size.py aa_permutation 02_01_global_peak_size/results_globalPeakSize_aaPerm
