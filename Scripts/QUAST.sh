#!/bin/bash
# Script to carry out quality assessment of de novo assembled genomes using QUAST
# Set up with multiple functions for different types of assembled data.
# Input data can be raw assemblies (rawQuas); decontaminated assemblies (cleanQuas); medaka corrected assemblies (meddyQuas);
# pilon corrected assemblies (pilQuas); repeat masked assemblies (maskQuas or maskQuasScript) or complete annotated genomes (genoQuas)
set -e
Assem=$1 # path to folder holding directories of isolate assemblies i.e path/to/all/isolate/assemblies that contains directories of isolate assemblies. Can be raw, medaka, decontaminated 
# input for masked assemblies and annotated genomes differ. These simply need the path to the folder containing isolate fasta files to be assessed
output=$2 # path to output directory
Refre=$3 # path to the reference genome to be compared against
RefGf=$4 # path to the reference genome's gff file
THREADS=$5 # number of threads

cd $Assem
# function to assess the raw de novo genomes.
rawQuas() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	fname="${fna%*.*}"
    echo $fname
	mkdir -p $output/${fname}/Raw
    pwd
    echo
    echo "Running Quast..."
    echo     
	#flye raw
	run="quast -t $THREADS ${Assem}/${fname}/assembly.fasta -o ${output}/${fname}/Raw/ -r ${Refre} -g ${RefGf} --large -f --circos"	
	echo $run
	eval $run
    cd $Assem
	"QUAST for ${fna} done successfully"
done
}
# function to assess the cleaned and decontaminated genomes.
cleanQuas() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	fname="${fna%*.*}"
    echo $fname
	mkdir -p $output/${fname}/Clean
    pwd
    echo
    echo "Running Quast..."
    echo 
	run="quast -t $THREADS ${Assem}/${fname}/clean_assembly.fasta -o ${output}/${fname}/Clean/ -r ${Refre} -g ${RefGf} --large -f --circos"
	echo $run
	eval $run
	cd $Assem
done
}
# function to assess the completed and annotated genomes.
genoQuas() {
for f in ${Assem}/*
do
    echo "The isolate you are working on is $f"
	echo $f
    fna=$(basename "$f")
	echo $fna
	fname="${fna%*.*}"
    echo $fname
    pwd
    echo
    echo "Running Quast..."
    echo 
    mkdir -p $output/${fname}/Annotated    
	quast -t $THREADS ${Assem}/${fname}.fasta -o ${output}/${fname}/Annotated -r ${Refre} -g ${RefGf} --large -f --circos
    cd $Assem
	done
}
# function to assess the medaka corrected genomes.
meddyQuas() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	fname="${fna%*.*}"
    echo $fname
	mkdir -p $output/${fname}/Medaka
    pwd
    echo
    echo "Running Quast..."
    echo 
	run="quast -t $THREADS ${Assem}/${fname}/${fname}_consensus.fasta -o ${output}/${fname}/Medaka/ -r ${Refre} -g ${RefGf} --large -f --circos"
	echo $run
	eval $run
	cd $Assem
done
}
# function to assess the pilon correct genomes
pilQuas() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	fname="${fna%*_*}"
    echo $fname
	mkdir -p $output/${fname}/Pilon
    pwd
    echo
    echo "Running Quast..."
    echo 
	run="quast -t $THREADS ${Assem}/${fname}_iter2/${fname}.fasta -o ${output}/${fname}/Pilon/ -r ${Refre} -g ${RefGf} --large -f --circos"
	echo $run
	eval $run
	cd $Assem
done
}
# function to assess the masked genomes. This function is tuned to be called from a specified quality assesment script
maskQuasScript() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	fname="${fna%*_*}"
    echo $fname
	outy=$output/$fname
	mkdir -p $outy
    pwd
    echo
    echo "Running Quast..."
    echo 
	run="quast -t $THREADS $folder -o $outy -r ${Refre} -g ${RefGf} --large -f --circos"
	echo $run
	eval $run
	cd $Assem
done
}
# function to assess the masked genomes. This function is tuned for this script and to be called from the command line
maskQuas() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	fname="${fna%*_*}"
    echo $fname
	mkdir -p $output/${fname}/Masked
	outy=$output/${fname}/Masked/
    pwd
    echo
    echo "Running Quast..."
    echo 
	run="quast -t $THREADS $folder -o $outy -r ${Refre} -g ${RefGf} --large -f --circos"
	echo $run
	eval $run
	cd $Assem
done
}
