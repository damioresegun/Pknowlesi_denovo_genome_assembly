#!/bin/bash
####
# Rationale: script to extract contigs into individual fasta files. Function is also present to delete these contigs from the fasta files
# The script is currently made specifically for contigs to take out API and MIT from fasta files. It can be easily adjusted for other use cases
# Script currently assumes that there are two text files per isolate. Also assumes that multiple fasta files will be manipulated
# These two text files are used individually for the extraction. The extract function works for any fasta and text files
# The delete function assumes more than one text file per isolate, combines these and then deletes the contigs from the file, saving it in the output location
###
scripts=/path/to/scripts/folder # path to directory for the script. NOT the script itself
extract() {
input=$1 #folder containing the text files with the contigs to remove e.g $HOME/API_MIT/ExtractedContigs/contigsToremove
fasta=$2 # the fasta file directory e.g $HOME/PolishedAndCorrected
output=$3 # directory for output e.g. $HOME/API_MIT/ExtractedContigs
for i in $input/*; 
do 
pyt=$(basename "$i"); 
typ=${pyt%.*}; 
echo $typ; 
typa=${pyt%_*.*}; 
python $scripts/extractFasta.py $input/${typ}.txt $fasta/${typa}.fasta $output/${typ}_extractedContigs.fasta;
done
}

delete() {
input=$1 #folder containing the text files with the contigs to remove e.g $HOME/API_MIT/ExtractedContigs/contigsToremove
fasta=$2 # the fasta file directory e.g $HOME/PolishedAndCorrected
output=$3 # directory for output e.g. $HOME/API_MIT_Free_Assemblies
mkdir -p $output
echo "First one ${fasta}"
for i in $input/*; 
do 
pyt=$(basename "$i"); 
typ=${pyt%.*}; 
echo $typ; 
typa=${pyt%_*.*}
cat $input/${typ}.txt >> $input/${typa}_contigs.txt;
done
echo "This is annoying ${fasta}"
for i in $input/*_contigs.txt; 
do pyt=$(basename "$i"); 
typ=${pyt%.*}; 
echo $typ; 
typa=${pyt%_*.*}; 
run="python3 $scripts/deleteFasta.py ${fasta}/${typa}.fasta $i > $output/${typa}.fasta; 
gzip $fasta/${typa}.fasta"; 
echo $run; 
eval $run;
done
}

