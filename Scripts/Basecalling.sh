#!/bin/bash
#$ -cwd
# Script to take in raw FAST5 nanopore files and take them through basecalling. 
# Currently set up for basecalling using guppy and then qcat demultiplexer using the demultiplexing script.
# Set up based on Guppy version 4.0.15
set -e
# Set Paths
Raw_IN=$1 # input folder for raw FAST5 files
echo "Raw files are in $Raw_IN"
Base_OUT=$2   # output folder for the basecalling
FLOWCELL=$3	# state the flowcell used e.g. FLO-MIN106
KIT=$4  # sequencing kit used. e.g. SQK-RBK004
# path to the ONT Guppy
ONT=$5 # path to the nanopore basecalling package. e.g. ont-guppy/bin
# Begin
mkdir -p $Base_OUT
mkdir -p $Dem_OUT

#guppy basecaller
gupbase=`$ONT/guppy_basecaller -i ${Raw_IN} --save_path ${Base_OUT} --flowcell $FLOWCELL --kit $KIT -r -v -q 0 --qscore_filtering -x auto`

echo $gupbase

echo "All done"