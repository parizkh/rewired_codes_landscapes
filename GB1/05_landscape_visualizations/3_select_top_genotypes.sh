cmd="conda activate gpmap; python filter_landscape.py"
qsub="qsub -cwd -l mem_free=8G"


python select_reduced_protein_sequences.py

codes=`cut -f 1 -d ','  ../genetic_codes/genetic_codes.csv | grep -v code`
for code in $codes
do
	echo "$cmd $code"  | $qsub -N filt.$code
done
