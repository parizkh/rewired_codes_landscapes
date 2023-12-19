data="protein_landscape.csv"
output="output"
n=10
Ns=1

cmd="calc_visualization -C $data -k $n -Ns $Ns"
qsub="qsub -cwd -l mem_free=8G -pe threads 4"

code="Standard"
echo "$cmd -c "Standard" -o $output/$code"  | $qsub -N gb1.$code
exit

codes=`cut -f 1 -d ','  genetic_codes.csv | grep -v code | grep -v Standard | tail -n 4`
for code in $codes
do
	echo "$cmd -c genetic_codes/code_$code.csv -o $output/$code"  | $qsub -N $code
done
