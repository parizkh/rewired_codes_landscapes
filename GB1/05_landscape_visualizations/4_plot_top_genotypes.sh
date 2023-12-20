cmd="conda activate gpmap; python plot_filt_landscape.py"
qsub="qsub -cwd -l mem_free=8G"

code="Standard"
echo "$cmd $code 2 1" | $qsub -N p$code -e logs/p$code.err -o logs/p$code.out
echo "$cmd $code 3 1" | $qsub -N p$code -e logs/p$code.err -o logs/p$code.out

code="9002"
echo "$cmd $code 1 2" | $qsub -N p$code -e logs/p$code.err -o logs/p$code.out
echo "$cmd $code 3 2" | $qsub -N p$code -e logs/p$code.err -o logs/p$code.out

code="73213"
echo "$cmd $code 2 1" | $qsub -N p$code -e logs/p$code.err -o logs/p$code.out
echo "$cmd $code 2 3" | $qsub -N p$code -e logs/p$code.err -o logs/p$code.out

code="5521"
echo "$cmd $code 1 2" | $qsub -N p$code -e logs/p$code.err -o logs/p$code.out

code="6037"
echo "$cmd $code 1 2" | $qsub -N p$code -e logs/p$code.err -o logs/p$code.out
echo "$cmd $code 3 2" | $qsub -N p$code -e logs/p$code.err -o logs/p$code.out

