cd input
ln -s ../../../01_vcregression/output/map.txt .
ln -s ../../aa_permutation/input/code_standard.tsv .
ln -s ../../aa_permutation/input/map.tsv .
cd ..

outFile="output/results"
codeFile="tmp_code"
tmpOut="tmp_out"

for i in $(seq 0 99999); do
        python3 ../../../code/generate_gen_code.py random $i $outFile $codeFile
        ../../../code/./landscape_ruggedness_ParD3 input/map.tsv $codeFile $tmpOut
        cat $tmpOut >> $outFile
done

# compute the global peak sizes for all codes
python3 ../aa_permutation/02_01_global_peak_size/gen_codes_global_peak_size.py random 02_01_global_peak_size/results_globalPeakSize_random
