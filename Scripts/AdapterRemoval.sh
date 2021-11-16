#!/bin/bash
#Script to take input fastq files through porechop based on the names they already have. Hence these may need renaming manually after this. 
#.................................................................................................................
#											SET PATHS BASED ON FLAGS
#.................................................................................................................
porechop=$1 # path to porechop-runner.py  script.
input=$2 # path to the demultiplexed fastq files to process
expr=$3 # name of the experiment stated by the user
output=$4 # output directory
THREADS=$5 # number of threads
echo "The path to porechop is: $porechop"
set -e
	#for file in $input/*/ ## use this one for guppy
	for file in $input/*   # use this one for qcat
	do
	echo $file
    echo "Beginning adapter removal for $file..."
	fasa=$(basename "$file")
    echo $fasa
	fasa="${fasa%.*}"
	echo $fasa
	echo $output
	$porechop -i $file -o $output/${fasa}.gz --verbosity 2 --threads $THREADS --format fastq.gz
done
echo "All adapters removed and the files are saved in $output"
echo "You may need to manually change the name of these files to your isolates..."
exit