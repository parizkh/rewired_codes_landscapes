size=5
edges="output/edges.npz"
cmd="conda activate gpmap; plot_visualization -H $size -W $size -nr 1000 -er 2000 --datashader -e $edges -nc function"
qsub="qsub -cwd -l mem_free=24G"


code="Standard"
nodes="output/$code.nodes.pq"
echo "$cmd $nodes -o plots/$code.2.1 -a 2,1" | $qsub -N p$code.2.1 
echo "$cmd $nodes -o plots/$code.3.1 -a 3,1" | $qsub -N p$code.3.1

code="9002"
nodes="output/$code.nodes.pq"
echo "$cmd $nodes -o plots/$code.1.2 -a 1,2" | $qsub -N p$code.1.2
echo "$cmd $nodes -o plots/$code.3.2 -a 3,2" | $qsub -N p$code.3.2

code="73213"
nodes="output/$code.nodes.pq"
echo "$cmd $nodes -o plots/$code.2.1 -a 2,1" | $qsub -N p$code.2.1
echo "$cmd $nodes -o plots/$code.2.3 -a 2,3" | $qsub -N p$code.2.3

code="5521"
nodes="output/$code.nodes.pq"
echo "$cmd $nodes -o plots/$code.1.2 -a 1,2" | $qsub -N p$code.1.2

code="6037"
nodes="output/$code.nodes.pq"
echo "$cmd $nodes -o plots/$code.1.2 -a 1,2" | $qsub -N p$code.1.2
echo "$cmd $nodes -o plots/$code.3.2 -a 3,2" | $qsub -N p$code.3.2

