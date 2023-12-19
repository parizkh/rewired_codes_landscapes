mkdir input
cd input
ln -s ../../../aa_permutation/input/map.tsv .
ln -s ../../../../../resources/code_standard.tsv .
cd ..
mkdir output

python3 ../../../../../code/compute_emp_robustness.py
