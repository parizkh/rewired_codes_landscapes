qsub="qsub -cwd -l mem_free=8G"
cmd="conda activate gpmap ; filter_genotypes -e output/edges.npz -n 167730 -l function -nf pq -ef npz"

codes=`cut -f 1 -d ','  genetic_codes.csv | grep -v code`
for code in $codes
do
        echo "$cmd output/$code.nodes.pq -o output/$code.filtered" | qsub -N $code.f
done

