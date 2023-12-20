cmd="conda activate gpmap ; filter_genotypes -e output/edges.npz -n 167730 -l function -nf pq -ef npz"
qsub="qsub -cwd -l mem_free=24G"

codes=`cut -f 1 -d ','  genetic_codes.csv | grep -v code`
for code in $codes
do
        echo "$cmd output/$code.nodes.pq -o output/$code.filtered" | $qsub -N f$code -e logs/f$code.err -o logs/f$code.out
done

