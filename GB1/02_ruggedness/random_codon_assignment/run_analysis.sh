cd input
ln -s ../../../01_vcregression/output/map.txt .
ln -s ../../aa_permutation/input/code_standard.tsv .
cat map.txt | sed 's/,/\t/g' > map.tsv
cd ..

for i in $(seq 0 49); do
	start=$(expr $i \* 2000)
	end=$(expr $start + 1999)
	sbatch -t 1-00 --mem-per-cpu=5000 --wrap="sh inner_loop.sh $start $end"
done

# compute the global peak sizes for all codes
python3 ../aa_permutation/02_01_global_peak_size/gen_codes_global_peak_size.py random 02_01_global_peak_size/results_globalPeakSize_random