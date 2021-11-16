#!/bin/bash
#$ -cwd

# Script to take in raw nanopore files and take them through demultiplexing.
# Set up for guppy demultiplexer and qcat demultiplexer. However qcat is recommended and left as default.
# Uncomment for the ones where necessary. Ensure that this is run in environments where the necessary tools are. 
# Set up for Qcat version 1.1.0

set -e
#Set Paths
Base_IN=$1   # input folder containing basecalled reads
Dem_OUT=$2    # output folder for the demultiplexing
THREADS=$3   # number of threads
#KIT=$4  # sequencing kit used. this version is for guppy e.g. SQK-RBK004
KIT=$4    # sequencing kit used. this version is for qcat e.g. RBK004
ONT=$5	# path to the ONT Guppy e.g. ont-guppy/bin

#mkdir -p $Dem_OUT
#............................................
#				Guppy Barcoder
#............................................
#echo "The files are in $Base_IN"
#gupbar="$ONT/guppy_barcoder -i ${Base_IN}/pass -s ${Dem_OUT}/Guppy --barcode_kits $KIT -q 0 --compress_fastq -x auto"
#echo $gupbar
#eval $gupbar
#pigz --best $Dem_OUT/Guppy/*.fastq
#............................................
#                 QCAT
#............................................
#echo "The files are in $Base_IN"
qdem="cat $Base_IN/pass/* | qcat -b ${Dem_OUT} --detect-middle -t $THREADS --trim -k '$KIT' --guppy"
echo $qdem
eval $qdem
pigz --best $Dem_OUT/*.fastq
echo "Demultiplexing done"