# prepare the input file
cd input
ln -s ../../00_data/output/GB1_data_processed.csv .
cd ..
cat input/GB1_data_processed.csv | cut -f 1,6,8 -d "," | tail -n +2 > vcregression_input

# run vcregression
python3 ~/vcregression-master/vcregression/vc_prep.py 20 4 -name vcreg  -data vcregression_input 
python3 ~/vcregression-master/vcregression/vc_map_estimate.py 20 4 -name vcreg -data vcregression_input -lambdas vcreg/lambdas.txt

# output
mkdir output
cp vcreg/map.txt output/
rm vcregression_input
