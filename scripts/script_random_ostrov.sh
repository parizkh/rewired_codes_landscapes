from=$1
to=$2
popSize=$3

mkdir -p output/results_"$popSize"

resultsFile="output/results_"$popSize"/results_"$from"-"$to
echo -n > $resultsFile
tmpFile="tmpFile_"$popSize"_"$from"-"$to
tmpCode="tmp_code_"$popSize"_"$from"-"$to
tmpResults="tmp_results_"$popSize"_"$from"-"$to

cat input/ostrov_codes_summary.tsv | head -n $to | tail -n +$from > $tmpFile

while read code block1 block2 block3 block4 X Z B J; do
	# create the code file
	python3 ../../../code/create_ostrov_code.py $block1 $block2 $block3 $block4 $tmpCode 0 $X $Z $B $J
	sed -i 's/O/*/g' $tmpCode
	# run the simulation 
	../../../code/./random_walk input/map.tsv $tmpCode $tmpResults $popSize
	
	if [ $? -eq 0 ]; then
	    cat $tmpResults >> $resultsFile
	else
	    exit 1
	fi

	rm $tmpCode
done < $tmpFile

paste $tmpFile $resultsFile > $tmpResults
mv $tmpResults $resultsFile

rm $tmpFile
