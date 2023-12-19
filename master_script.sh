#!/bin/bash

# parameters
name=$1	#name of the protein
N=$2	#number of loci

cd $name

: '
######## 1: vcregression
echo "Step 1: vcregression"
mkdir 01_vcregression
cd 01_vcregression

# input
mkdir input
cd input
ln -s ../../00_data/output/"$name"_data_processed.csv .
cd ..

# run vcregression
python3 ~/vcregression-master/vcregression/vc_prep.py 20 $N -name vcreg  -data input/"$name"_data_processed.csv 
python3 ~/vcregression-master/vcregression/vc_map_estimate.py 20 $N -name vcreg -data input/"$name"_data_processed.csv -lambdas vcreg/lambdas.txt

# output
mkdir output
cp vcreg/map.txt output/
cd ..

echo "vcregression done"
echo "#####"
'

######## 2: ruggedness
echo "Step 2: ruggedness"
mkdir 02_ruggedness
cd 02_ruggedness

mkdir aa_permutation aa_permutation_restricted ostrov random_codon_assignment

echo "2.1: aa permutation codes"
cd aa_permutation
# prepare input
mkdir input
cd input
ln -s ../../../01_vcregression/output/map.txt .
ln -s ../../../../resources/code_standard.tsv .
cat map.txt | sed 's/,/\t/g' > map.tsv
cd ..

if [[ $N -eq 3 ]]; then
	sbatch -t 1-00 --wrap="sh ../../../scripts/script_ruggedness.sh 0 99999 $N aa_permutation"
else
	for i in $(seq 0 99); do
		start=$(expr $i \* 1000)
		end=$(expr $start + 999)
		sbatch -t 1-00 --mem-per-cpu=5000 --wrap="sh ../../../scripts/script_ruggedness.sh $start $end $N aa_permutation"
	done
fi
cd ..

echo "2.2 restricted aa permutation cdes"
cd aa_permutation_restricted
# prepare input
mkdir input
cd input
ln -s ../../aa_permutation/input/map.tsv .
ln -s ../../../../resources/code_standard.tsv .
cd ..

if [[ $N -eq 3 ]]; then
	sbatch -t 1-00 --wrap="sh ../../../scripts/script_ruggedness.sh 0 99999 $N aa_permutation_restricted"
else
	for i in $(seq 0 99); do
		start=$(expr $i \* 1000)
		end=$(expr $start + 999)
		sbatch -t 1-00 --mem-per-cpu=5000 --wrap="sh ../../../scripts/script_ruggedness.sh $start $end $N aa_permutation_restricted"
	done
fi
cd ..

echo "2.3 random codon assignment cdes"
cd random_codon_assignment
# prepare input
mkdir input
cd input
ln -s ../../aa_permutation/input/map.tsv .
ln -s ../../../../resources/code_standard.tsv .
cd ..

if [[ $N -eq 3 ]]; then
	sbatch -t 1-00 --wrap="sh ../../../scripts/script_ruggedness.sh 0 99999 $N random"
else
	for i in $(seq 0 99); do
		start=$(expr $i \* 1000)
		end=$(expr $start + 999)
		sbatch -t 1-00 --mem-per-cpu=5000 --wrap="sh ../../../scripts/script_ruggedness.sh $start $end $N random"
	done
fi
cd ..

echo "2.4 ostrov"
cd ostrov
# prepare input
mkdir input
cd input
ln -s ../../aa_permutation/input/map.tsv .
ln -s ../../../../resources/code_standard.tsv .
ln -s ../../../../resources/ostrov_codes_summary.tsv .
cd ..

if [[ $N -eq 3 ]]; then
	for i in $(seq 0 3); do
		start=$(expr $i \* 50000 + 2)
		end=$(expr $start + 49999)
		sbatch -t 1-00 --wrap="sh ../../../scripts/script_ruggedness_ostrov.sh $start $end $N"
	done
else
	for i in $(seq 0 129); do
		start=$(expr $i \* 1500 + 2)
		end=$(expr $start + 1499)
		sbatch -t 5-00 --mem-per-cpu=5000 --wrap="sh ../../../scripts/script_ruggedness_ostrov.sh $start $end $N"
	done
fi
cd ..

cd ..
echo 


###### 3: greedy walks
echo "#####"
echo "Step 3: greedy walks"
mkdir 03_greedy_walks
cd 03_greedy_walks

mkdir aa_permutation aa_permutation_restricted ostrov random_codon_assignment pines_calles_codes

echo "3.1: aa permutation codes"
cd aa_permutation
# prepare input
mkdir input
cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../../resources/code_standard.tsv .
cd ..

if [[ $N -eq 3 ]]; then
	for i in $(seq 0 9); do
		start=$(expr $i \* 10000)
		end=$(expr $start + 9999)
		sbatch -t 1-00 --wrap="sh ../../../scripts/script_greedy.sh $start $end $N aa_permutation"
	done
else
	for i in $(seq 0 99); do
		start=$(expr $i \* 1000)
		end=$(expr $start + 999)
		sbatch -t 1-00 --mem-per-cpu=5000 --wrap="sh ../../../scripts/script_greedy.sh $start $end $N aa_permutation"
	done
fi
cd ..

echo "3.2: restricted aa permutation codes"
cd aa_permutation_restricted
# prepare input
mkdir input
cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../../resources/code_standard.tsv .
cd ..

if [[ $N -eq 3 ]]; then
	for i in $(seq 0 9); do
		start=$(expr $i \* 10000)
		end=$(expr $start + 9999)
		sbatch -t 1-00 --wrap="sh ../../../scripts/script_greedy.sh $start $end $N aa_permutation_restricted"
	done
else
	for i in $(seq 0 99); do
		start=$(expr $i \* 1000)
		end=$(expr $start + 999)
		sbatch -t 1-00 --mem-per-cpu=5000 --wrap="sh ../../../scripts/script_greedy.sh $start $end $N aa_permutation_restricted"
	done
fi
cd ..

echo "3.3: random codon assignment codes"
cd random_codon_assignment
# prepare input
mkdir input
cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../../resources/code_standard.tsv .
cd ..

if [[ $N -eq 3 ]]; then
	for i in $(seq 0 9); do
		start=$(expr $i \* 10000)
		end=$(expr $start + 9999)
		sbatch -t 1-00 --wrap="sh ../../../scripts/script_greedy.sh $start $end $N random"
	done
else
	for i in $(seq 0 99); do
		start=$(expr $i \* 1000)
		end=$(expr $start + 999)
		sbatch -t 1-00 --mem-per-cpu=5000 --wrap="sh ../../../scripts/script_greedy.sh $start $end $N random"
	done
fi
cd ..

echo "3.4: ostrov"
cd ostrov
# prepare input
mkdir input
cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../../resources/code_standard.tsv .
ln -s ../../../../resources/ostrov_codes_summary.tsv .
cd ..

if [[ $N -eq 3 ]]; then
	for i in $(seq 0 4); do
		start=$(expr $i \* 40000 + 2)
		end=$(expr $start + 39999)
		sbatch -t 5-00 --mem-per-cpu=10000 --wrap="sh ../../../scripts/script_greedy_ostrov.sh $start $end $N"
	done
else
	for i in $(seq 0 97); do
		start=$(expr $i \* 2000 + 2)
		end=$(expr $start + 1999)
		sbatch -t 5-00 --mem-per-cpu=10000 --wrap="sh ../../../scripts/script_greedy_ostrov.sh $start $end $N"
	done
fi
cd ..

echo "3.5 Pines and Calles codes"
cd pines_calles_codes
# prepare input
mkdir input
cd input
mkdir codes
cd codes
ln -s ../../../../../resources/pines_calles_codes/* .
cd ..
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
cd ..

sbatch -t 1-00 --mem-per-cpu=10000 --wrap="sh ../../../scripts/script_greedy_pines_calles.sh $N"
cd ..

cd ..
echo 


###### 4: random walks
echo "#####"
echo "Step 4: random walks"
mkdir 04_random_walks
cd 04_random_walks

mkdir aa_permutation ostrov random_codon_assignment

echo "4.1: aa permutation codes"
cd aa_permutation
# prepare input
mkdir input
cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../../resources/code_standard.tsv .
cd ..

for i in $(seq 0 99); do
	for popSize in 10 100 10000 1000000; do
		start=$(expr $i \* 1000)
		end=$(expr $start + 999)
		sbatch -t 1-00 --wrap="sh ../../../scripts/script_random.sh $start $end $popSize aa_permutation"
	done
done

cd ..

echo "4.2: random codon assignment codes"
cd random_codon_assignment
# prepare input
mkdir input
cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../../resources/code_standard.tsv .
cd ..

for i in $(seq 0 99); do
	for popSize in 10 100 10000 1000000; do
		start=$(expr $i \* 1000)
		end=$(expr $start + 999)
		sbatch -t 1-00 --wrap="sh ../../../scripts/script_random.sh $start $end $popSize random"
	done
done

cd ..

echo "4.3: Ostrov codes"
cd ostrov
# prepare input
mkdir input
cd input
ln -s ../../../02_ruggedness/aa_permutation/input/map.tsv .
ln -s ../../../../resources/code_standard.tsv .
ln -s ../../../../resources/ostrov_codes_summary.tsv .
cd ..

for i in $(seq 0 97); do
	for popSize in 10 100 10000 1000000; do
		start=$(expr $i \* 2000 + 2)
		end=$(expr $start + 1999)
		sbatch -t 1-00 --wrap="sh ../../../scripts/script_random_ostrov.sh $start $end $popSize"
	done
done

cd ..

echo
echo "#############"
echo "Done"
echo "#############"