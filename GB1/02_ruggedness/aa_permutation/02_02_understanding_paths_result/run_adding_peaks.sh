for i in 1 5 10 50 100 200 300 400 500 1000 2000 5000 10000; do		# number of peaks added
	for j in $(seq 100); do 		# random seed
		python3 add_random_peaks.py ../input/map.tsv $i $j
		../../../../code/./landscape_ruggedness tmp_input_peaks ../input/code_standard.tsv tmp_out_peaks
		echo -n $i $j " " >> output_peaks
		cat tmp_out_peaks >> output_peaks
	done
done

sed -i 's/ /\t/g' output_peaks