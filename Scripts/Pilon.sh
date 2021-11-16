#!/bin/bash
# Script to carry out rounds of pilon cleaning after medaka IF illumina reads are present
# requires forward and reverse illumina reads, the output of racon or medaka, the number of threads, the output directory and output prefix
set -e
illum=$1 # path to the directory holding folders of isolate Illumina reads
assem=$2 # path to the medaka output folder holding all folders of isolates
output=$3 # path to output pilon results
alignout=$4 # path to output alignments of each iteration against the illumina reads
THREADS=$5 # number of reads
#mkdir -p $alignout

#load bwa for alignment
#module load bwa
for i in $illum/*; 
do
	#iteration one
	fname=$(basename "$i")
	echo $fname
	mkdir -p $alignout/${fname}
	mkdir -p $output/${fname}_iter0
	#carry out alignment of reads against racon/medaka output
	#index first
	bwa index $assem/${fname}/${fname}_consensus.fasta
	#align
	bwa mem $assem/${fname}/${fname}_consensus.fasta $illum/${fname}/*_1.fastq $illum/${fname}/*_2.fastq > $alignout/${fname}/${fname}VsPolishedAssem_0.sam -t $THREADS
	echo "bwa alignment done"
	### samtools
	samtools view -bS -@ $THREADS $alignout/${fname}/${fname}VsPolishedAssem_0.sam > $alignout/${fname}/${fname}VsPolishedAssem_0.bam
	samtools sort -@ $THREADS -o $alignout/${fname}/${fname}VsPolishedAssem_0_sorted.bam $alignout/${fname}/${fname}VsPolishedAssem_0.bam
	rm $alignout/${fname}/${fname}VsPolishedAssem_0.sam
	rm $alignout/${fname}/${fname}VsPolishedAssem_0.bam
	mv $alignout/${fname}/${fname}VsPolishedAssem_0_sorted.bam $alignout/${fname}/${fname}VsPolishedAssem_0.bam
	samtools index $alignout/${fname}/${fname}VsPolishedAssem_0.bam
	#### pilon
	pilon -Xmx120G --genome $assem/${fname}/${fname}_consensus.fasta --bam $alignout/${fname}/${fname}VsPolishedAssem_0.bam --threads $THREADS --outdir $output/${fname}_iter0/ --output ${fname} --tracks --fix all,circles
	echo "iteration one done"
	echo "starting interation two"
###########################################################################
	#index again
	bwa index $output/${fname}_iter0/${fname}.fasta
	#align the round 1 correction with the short reads for round 2
	mkdir -p $output/${fname}_iter1
	bwa mem $output/${fname}_iter0/${fname}.fasta $illum/${fname}/*_1.fastq $illum/${fname}/*_2.fastq > $alignout/${fname}/${fname}VsPolishedAssem_1.sam -t $THREADS
	
	#samtools
	samtools view -@ $THREADS -bS $alignout/${fname}/${fname}VsPolishedAssem_1.sam > $alignout/${fname}/${fname}VsPolishedAssem_1.bam
	samtools sort -@ $THREADS -o $alignout/${fname}/${fname}VsPolishedAssem_1_sorted.bam $alignout/${fname}/${fname}VsPolishedAssem_1.bam
	rm $alignout/${fname}/${fname}VsPolishedAssem_1.sam
	rm $alignout/${fname}/${fname}VsPolishedAssem_1.bam
	mv $alignout/${fname}/${fname}VsPolishedAssem_1_sorted.bam $alignout/${fname}/${fname}VsPolishedAssem_1.bam
	samtools index $alignout/${fname}/${fname}VsPolishedAssem_1.bam
	#### pilon
	pilon -Xmx120G --genome $output/${fname}_iter0/${fname}.fasta --bam $alignout/${fname}/${fname}VsPolishedAssem_1.bam --threads $THREADS --outdir $output/${fname}_iter1/ --output ${fname} --tracks --fix all,circles
	echo "iteration two done"
	echo "iteration three starting"
##################################################################################
	#index again
	bwa index $output/${fname}_iter1/${fname}.fasta
	#align the round 1 correction with the short reads for round 2
	mkdir -p $output/${fname}_iter2
	bwa mem $output/${fname}_iter1/${fname}.fasta $illum/${fname}/*_1.fastq $illum/${fname}/*_2.fastq > $alignout/${fname}/${fname}VsPolishedAssem_2.sam -t $THREADS
	
	#samtools
	samtools view -@ $THREADS -bS $alignout/${fname}/${fname}VsPolishedAssem_2.sam > $alignout/${fname}/${fname}VsPolishedAssem_2.bam
	samtools sort -@ $THREADS -o $alignout/${fname}/${fname}VsPolishedAssem_2_sorted.bam $alignout/${fname}/${fname}VsPolishedAssem_2.bam
	rm $alignout/${fname}/${fname}VsPolishedAssem_2.sam
	rm $alignout/${fname}/${fname}VsPolishedAssem_2.bam
	mv $alignout/${fname}/${fname}VsPolishedAssem_2_sorted.bam $alignout/${fname}/${fname}VsPolishedAssem_2.bam
	samtools index $alignout/${fname}/${fname}VsPolishedAssem_2.bam
	#### pilon
	pilon -Xmx120G --genome $output/${fname}_iter1/${fname}.fasta --bam $alignout/${fname}/${fname}VsPolishedAssem_2.bam --threads $THREADS --outdir $output/${fname}_iter2/ --output ${fname} --tracks --fix all,circles
	echo "iteration three done"
##################################################################################
	echo "Pilon done"
done
