conda activate gpmap

# Calculate visualziation
calc_visualization protein_landscape.csv -k 20 -Ns 1 -o output/protein

# Plot landscape
python plot_protein_landscape.py
