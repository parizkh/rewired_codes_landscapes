# prepare the input
cd input
ln -s ../../../01_vcregression/output/map.txt .
ln -s ../../../../GB1/02_ruggedness/aa_permutation/input/code_standard.tsv .
cd ..
cat input/map.txt | sed 's/,/\t/g' > input/map.tsv

# run the ruggedness analysis
mkdir -p "output/results"
outFile="output/results/results"
codeFile="tmp_code"
tmpOut="tmp_out"

for i in $(seq 0 99999); do
        python3 ../../../code/generate_gen_code.py aa_permutation $i $outFile $codeFile
        ../../../code/./landscape_ruggedness_ParD3 input/map.tsv $codeFile $tmpOut
        cat $tmpOut >> $outFile
done

rm $codeFile
rm $tmpOut

# compute the global peak sizes for all codes
python3 02_01_global_peak_size/gen_codes_global_peak_size.py aa_permutation 02_01_global_peak_size/results_globalPeakSize_aaPerm
