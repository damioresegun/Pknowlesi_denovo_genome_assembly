#!/bin/bash
# Script to take racon polished assemblies through medaka. 
rawsPath=$1 	# add the path to the raw fastq files used for de-novo
raconOut=$2	# add path to the racon outputs
outputPath=$3	# add path to output for medaka
QcAssem=$4	# add path to place assembly stats
THREADS=$5 # number of threads
mkdir -p $outputPath
mkdir -p $QcAssem

for i in $raconOut/* ; 
do 
fname=$(basename "$i")
echo "You are working on the ${fname} folder"
mkdir -p $outputPath/${fname}
#cd $outputPath/${fname}
#no more than 8 threads
medaka_consensus -i $rawsPath/${fname}*.fastq -d $raconOut/${fname}/${fname}_iter4.fasta -o $outputPath/${fname} -t $THREADS
#rename the output
mv $outputPath/${fname}/consensus.fasta $outputPath/${fname}/${fname}_consensus.fasta
# take some assembly stats of this
assembly-stats $outputPath/${fname}/${fname}_consensus.fasta >> ${QcAssem}/MaxGuppyAssemStats.txt
echo "Isolate ${fname} done"
done
echo "All done"