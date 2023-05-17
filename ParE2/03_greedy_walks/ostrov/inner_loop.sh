from=$1
to=$2
inputFile=$3

mkdir -p output

resultsFile="output/results_"$from"-"$to
echo -n > $resultsFile
tmpFile="tmpFile_"$from"-"$to
tmpCode="tmp_code_"$from"-"$to
tmpResults="tmp_results_"$from"-"$to

cat input/ostrov_codes_summary.tsv | head -n $to | tail -n +$from > $tmpFile

while read code block1 block2 block3 block4 X Z B J; do
	# create the code file
	python3 ../../../code/create_ostrov_code.py $block1 $block2 $block3 $block4 $tmpCode 0 $X $Z $B $J
	sed -i 's/O/*/g' $tmpCode
	# run the simulation 
	../../../code/./greedy_walk_ParD3 $inputFile $tmpCode $tmpResults
	
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
