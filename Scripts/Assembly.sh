#!/bin/bash
#.................................................................................................................
#											SET PATHS BASED ON FLAGS
#.................................................................................................................
### Entire adapter removed data
input=$1	# raw input reads to assemble on flye. Can be fasta or fastq
output=$2	# output directory
THREADS=$3	# number of threads
nmr=$4	# isolate name for prefix for a specific isolate
geno=$5	# estimated genome-size

    echo "The file you are working on is $input"
    mkdir -p $output/${nmr}
	outr=$output/${nmr}
	flye --nano-raw ${input} --genome-size $geno --out-dir ${outr} --threads $THREADS
exit