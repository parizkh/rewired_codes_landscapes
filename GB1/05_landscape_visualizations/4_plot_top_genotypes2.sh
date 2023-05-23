cmd="conda activate gpmap; python plot_filt_landscape.py"
qsub="qsub -cwd -l mem_free=8G"

code="Standard"
echo "$cmd $code 2 1 3" | $qsub -N p$code
echo "$cmd $code 3 1 2" | $qsub -N p$code

code="9002"
echo "$cmd $code 1 2 3" | $qsub -N p$code
echo "$cmd $code 3 2 1" | $qsub -N p$code

code="73213"
echo "$cmd $code 2 1 3" | $qsub -N p$code
echo "$cmd $code 2 3 1" | $qsub -N p$code

code="5521"
echo "$cmd $code 1 2 3" | $qsub -N p$code

code="6037"
echo "$cmd $code 1 2 3" | $qsub -N p$code
echo "$cmd $code 3 2 1" | $qsub -N p$code

