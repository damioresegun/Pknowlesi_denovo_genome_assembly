#!/bin/bash
input=$1 
output=$2
mkdir -p $output
for f in $input/*
do
fna=$(basename "$f")
fname=${fna%.*}
mkdir -p $output/${fname}
cd $output/${fname}
csplit -f chr -s -z ${f} '/>/' '{*}'

#rename the new files with the chromosome names and save as fasta files
for i in chr* ; do
  head -n 1 ${i}
  n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ;
  mv ${i} "${n}.fa" ;
 done
 done

