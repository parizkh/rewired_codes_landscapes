from=$1
to=$2

logFile="output/log_"$from"-"$to
echo -n > $logFile
resultsFile="output/results_"$from"-"$to
echo -n > $resultsFile
tmpFile="tmpFile_"$from"-"$to
tmpCode="tmp_code_"$from"-"$to
tmpResults="tmp_results_"$from"-"$to


cat input/ostrov_codes_summary.tsv | head -n $to | tail -n +$from > $tmpFile

while read code block1 block2 block3 block4 X Z B J rob; do
	# create the code file
	python3 ../../../code/create_ostrov_code.py $block1 $block2 $block3 $block4 $tmpCode 0 $X $Z $B $J
	sed -i 's/O/*/g' $tmpCode
	
	# ruggedness
	echo '../../../code/./landscape_ruggedness_ParD3 ' input/map.tsv $tmpCode $tmpResults >> $logFile
	../../../code/./landscape_ruggedness_ParD3 input/map.tsv $tmpCode $tmpResults >> $logFile
	
	if [ $? -eq 0 ]; then
	    cat $tmpResults >> $resultsFile
		echo >> $logFile
	else
	    exit 1
	fi

	rm $tmpCode
	
done < $tmpFile

paste $tmpFile $resultsFile > $tmpResults
mv $tmpResults $resultsFile

rm $tmpFile
