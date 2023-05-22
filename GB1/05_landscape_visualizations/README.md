# Visualization of GB1 landscapes under rewired genetic codes

### Dependencies

Visualizations and plots were generated using [gpmap-tools](https://gpmap-tools.readthedocs.io) using commit XXX.
Set up a new python environment with conda and install the library. Use the same exact commit to ensure reproducibility

```
conda create -n gpmap python=3.7
conda activate gpmap
git clone https://github.com/cmarti/gpmap-tools
git checkouot XXX
cd gpmap-tools
python setup.py install
```

### Genetic codes

The file `genetic_codes.csv` contains the information about the genetic codes that were visualized, and the corresponding codes with the
 mapping of codons to aminoacids is found in the folder `genetic codes`

| code     | label        |
|----------|--------------|
| Standard | standard     |
| 9002     | robust_A     |
| 73213    | robust_B     |
| 6037     | non_robust_A |
| 5521     | non_robust_B |

The following codon table represents robust code A representing one of the most robust codes in the set of randomly generated codes, whe
re we can easily appreciate the new connectivity between the different aminoacids compared to the standard genetic code

|   |  T      |  C      |  A      |  G      |   |
|---|---------|---------|---------|---------|---|
| T | TTT G   | TCT T   | TAT R   | TGT P   | T |
| T | TTC G   | TCC T   | TAC R   | TGC P   | C |
| T | TTA Y   | TCA T   | TAA Stop| TGA Stop| A |
| T | TTG Y   | TCG T   | TAG Stop| TGG E   | G |
|---|---------|---------|---------|---------|---|
| C | CTT Y   | CCT Q   | CAT H   | CGT M   | T |
| C | CTC Y   | CCC Q   | CAC H   | CGC M   | C |
| C | CTA Y   | CCA Q   | CAA K   | CGA M   | A |
| C | CTG Y   | CCG Q   | CAG K   | CGG M   | G |
|---|---------|---------|---------|---------|---|
| A | ATT I   | ACT C   | AAT D   | AGT T   | T |
| A | ATC I   | ACC C   | AAC D   | AGC T   | C |
| A | ATA I   | ACA C   | AAA V   | AGA M   | A |
| A | ATG F   | ACG C   | AAG V   | AGG M   | G |
|---|---------|---------|---------|---------|---|
| G | GTT W   | GCT S   | GAT N   | GGT L   | T |
| G | GTC W   | GCC S   | GAC N   | GGC L   | C |
| G | GTA W   | GCA S   | GAA A   | GGA L   | A |
| G | GTG W   | GCG S   | GAG A   | GGG L   | G |
|---|---------|---------|---------|---------|---|

### Calculating coordinates of the visualization under rewired codes

We precalculate these coordinates because with such a big landscape containing over 16 million genotypes it can take up to a couple of hours and lots of memory. We recommend to run these calculations in parallel using a high performance computing cluster requiring at least 20GB per job.

```
#Create link to inferred landscape
ln -s ../01_vcregression/output/map.txt protein_landscape.csv

# Calculate coordinates for standard genetic code
cmd="calc_visualization protein_landscape.csv -k 10 -Ns 1 -C -ef npz -nf pq"
$cmd -o output/standard -e

# Rename edges file to use it for all visualizations
mv output/Standard.edges.npz output/edges.npz

# Calculate coordinates for other codes
codes=`cut -f 1 -d ','  genetic_codes.csv | grep -v code | grep -v Standard`
for code in $codes
do
	$cmd -c genetic_codes/code_$code.csv -o output/$code
done
```

This way, we generate files with the coordinates for 10 Diffusion axis for each possible 12 nucleotide sequence and a single file encoding which genotypes are connected to plot the edges between them

