#!/bin/bash
# Script to carry out formatting of fastq files that are being directed to. 
# Rationale: often fastq files might have some odd encoding, line ends etc. This script seeeks those out and attempts to fix it by calling the tool script convert_fq_to_fa.py.
## has to be run in the Flye or Canu conda environment as it needs bioawk and seqkit! 
#
scpts=$1 # path to directory holding the convert_fq_to_fa script
inp=$2	# path to fastq file to format
#outr="$HOME/All_Isolate_PkMapped_Formatted/"
outr=$3	# path to output directory for formatted fastqs
isol=$4	# isolate name that will be used to name the output file
echo
echo "The file you are working on is $inp"
    # format the file
    python3 $scpts/convert_fq_to_fa.py -i $inp -o $outr/${isol}.fastq
	# remove empty reads
	bioawk -cfastx 'length($seq) > 1 {print "@"$name"\n"$seq"\n+\n"$qual}' $outr/${isol}.fastq > $outr/${isol}_empty0.fastq
	#remove duplicates
	rm $outr/${isol}.fastq
    seqkit rename $outr/${isol}_empty0.fastq > $outr/${isol}.fastq
	rm $outr/${isol}_empty0.fastq
    echo "Formatting done"