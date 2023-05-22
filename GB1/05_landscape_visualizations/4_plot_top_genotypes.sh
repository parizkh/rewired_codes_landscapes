cmd="conda activate gpmap; python plot_filt_landscape.py"
qsub="qsub -cwd -l mem_free=8G"


codes=`cut -f 1 -d ','  ../genetic_codes/genetic_codes.csv | grep -v code`

for filt in filt1 # filt2
do
	for code in $codes
	do
		echo "$cmd $code $filt" | $qsub -N $filt.$code -o logs/$filt.$code.out
	done
done
