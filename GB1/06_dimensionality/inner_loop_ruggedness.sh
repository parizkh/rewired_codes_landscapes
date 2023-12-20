input=$1
L=$2

input_name=$( echo $input | sed 's%input/%%' )
echo $input_name
outFile="output_ruggedness/output_"$input_name
codeFile="tmp_code_ruggedness_"$input_name
tmpOut="tmp_out_ruggedness_"$input_name

for i in $(seq 0 99999); do
        python3 ../../code/generate_gen_code.py aa_permutation_restricted $i $outFile $codeFile
        ../../code/./landscape_ruggedness $input $codeFile $tmpOut $L
        cat $tmpOut >> $outFile
done

rm $codeFile
rm $tmpOut

# here compute the correlation and report to a file (problem with overwriting?)
python3 ../../code/compute_corrs.py $input_name $outFile "output_ruggedness/results"
#rm $outFile