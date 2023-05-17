for i in 1 10 100 1000 2000 3000 5000 7000 10000 50000 100000; do 		# number of reciprocal-sign epistatic squares
	for j in $(seq 100); do 		# random seed
		python3 add_random_rsEpi.py ../input/map.tsv $i $j
		../../../../code/./landscape_ruggedness_new tmp_input_rsEpi ../input/code_standard.tsv tmp_out_rsEpi
		echo -n $i $j " " >> output_rsEpi
		cat tmp_out_rsEpi >> output_rsEpi
	done
done

sed -i 's/ /\t/g' output_rsEpi