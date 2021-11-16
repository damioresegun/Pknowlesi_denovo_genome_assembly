#!/bin/bash
# Script to take data through part of decontamination after blobtools. 
# Should go through the blobDB text file to keep only entries that have Apicomplexa, no-hit, undef 
# thus theoretically keeping only the phylum for Apicomplexa and other undefined sequences. 
# After this, the script will get the names of the remaining sequences and put them in a list to then use to filter the contigs file. 
# The resulting contig file will then be taken through blobtools again for confirmation of decontamination before being taken forward for polishing

inputFol=$1 # path to the raw blobtools output folder
#path to the assembly fassta file
fastaPath=$2 # path to folder holding the assembly output folders i.e. 'path/to/FlyeOutput' containing assembly folders for 'isolate_A', 'isolate_B' etc
scrpts=$3 # path to scripts folder

for i in $inputFol/* ;
do
ls $i
cd $i
pwd

# take out the lines in the table that are not the header, Apicomplexa, no-hits and undefined
grep -E 'Apicomplexa|#|no-hit|undef' blobDB.table.txt > clean.blobDB.table.txt

#make a list of the contigs we want to keep in our fasta
cut -f 1,1 clean.blobDB.table.txt > nodes.txt

# take off the header in the nodes list
grep -v '^#' nodes.txt > node_names.txt

# remove the first nodes files
rm nodes.txt
cd
done
echo "Nodes cleaned. Moving on to cleaning the FASTA..."

for fas in $inputFol/* ;
do
ls $fas
fname=$(basename "$fas")
echo $fname
# call the fasta tool for the cleaning of the fasta
python3 $scrpts/FastaTool.py ${fastaPath}/${fname}/assembly.fasta ${inputFol}/${fname}/node_names.txt > ${fastaPath}/${fname}/clean_assembly.fasta
done
echo "All cleaning done. Felicitations!"

