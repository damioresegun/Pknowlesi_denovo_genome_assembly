#!/bin/bash
## Script to take input de novo assembled fasta files through assembly stats. 
## It can be infinitely expanded to take in fasta files of different types. 
## It is currently made for aseembly outputs of Flye. 
## needs the input folder to be the folder holding the folders for the different isolates. 

### input files
INP=$1 # path to folder holding folders of each isolate's assembly. i.e. p'ath/to/AssemblerOutput/' which contains folders like 'isolate_A', isolate_B' etc
OUTP=$2 #path to folder to store output txt files
whichr=$3 # Prefix for single isolate (e.g isolate name) or for all isolates (e.g. experiment name)
#isolates=$4 
mkdir -p $OUTP

# flye assembly stats
if [[ $whichr == "Some" ]]; then
	for i in ${isolates[*]};
	do
	echo "The folder you are working on is $i"
		fname=$i    
		echo $fname
		echo "cool"
		assembly-stats $INP/${fname}/assembly.fasta >> $OUTP/${whichr}_FlyeVsPk_assembly.txt
	done
	#convert flye assem to csv
	cat $OUTP/${whichr}_FlyeVsPk_assembly.txt | tr -s '[:blank:]' ',' > $OUTP/${whichr}_FlyeVsPk_assembly.csv
elif [[ $whichr == "All" ]]; then
	for i in $INP/*
	do
	ls $i
	echo "The folder you are working on is $i"
		fname=$(basename "$i")    
		echo $fname
		assembly-stats $INP/${fname}/assembly.fasta >> $OUTP/${whichr}_FlyeVsPk_assembly.txt
	done
	#convert flye assem to csv
	cat $OUTP/${whichr}_FlyeVsPk_assembly.txt | tr -s '[:blank:]' ',' > $OUTP/${whichr}_FlyeVsPk_assembly.csv
fi

echo "All done"
