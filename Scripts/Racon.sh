#!/bin/bash
# Script to take through data for racon polishing. Will carry out racon polishing through 4 rounds in order to take this through to Medaka. 

#output to a log file -- only if not using queue system
#exec 1>All_DeN_FlyeRacon.out 2>&1

contigsPath=$1	# path to directory holding folders of raw isolate de novo assembly output
rawsPath=$2	# path to directory holding raw fastq/fasta read files used to generate the de novo assemblies
outputPath=$3 # path to output racon outputs alignments 
alignPath=$4 # 
iterVread=$5 # folder to output intermediate alignments of each iteration versus the raw reads. Ideally a different folder from the final output directory
raconPath="$(which racon)"
THREADS=$6 # number of threads
Stats=$7 # path to output statistical metrics
set -e

#mkdir -p $iterVread
#mkdir -p $outputPath
for i in $contigsPath/* ;
do
fname=$(basename "$i")
echo "You are working on the ${fname} folder"
mkdir -p $outputPath/${fname}
cd $outputPath/${fname}
# Do iteration one for each isolate
$raconPath $rawsPath/${fname}*.fastq ${alignPath}/${fname}_readsVsCleanAssem.sam $contigsPath/${fname}/clean_assembly.fasta -t $THREADS > ${fname}_iter1.fasta
echo "Iteration one done for ${fname}"
pwd
echo "#######################################################################################"
echo "Moving to iteration two"
echo "#######################################################################################"
# aligning output of iteration 1 to the raw reads
minimap2 -ax map-ont ${fname}_iter1.fasta ${rawsPath}/${fname}*.fastq > ${iterVread}/${fname}_readsVsRacon1.sam -t $THREADS
    echo
    echo "Alignment done"
$raconPath $rawsPath/${fname}*.fastq ${iterVread}/${fname}_readsVsRacon1.sam ${fname}_iter1.fasta -t $THREADS > ${fname}_iter2.fasta
# cleaning up. converting and sorting sam file and zipping it
    #convert sam file to bam file
    samtools view --threads $THREADS -bS ${iterVread}/${fname}_readsVsRacon1.sam > ${iterVread}/${fname}_readsVsRacon1.bam
    # sort the bam file
    samtools sort --threads $THREADS -o ${iterVread}/${fname}_readsVsRacon1_sorted.bam ${iterVread}/${fname}_readsVsRacon1.bam
    # clean up
    echo "Deleting intermediate alignment files..."
    rm ${iterVread}/${fname}_readsVsRacon1.sam
    rm ${iterVread}/${fname}_readsVsRacon1.bam
    mv ${iterVread}/${fname}_readsVsRacon1_sorted.bam ${iterVread}/${fname}_readsVsRacon1.bam
    echo
    echo "SAM to BAM conversion done. Sorting of BAM file done."
    echo 
    echo "Iteration two done for ${fname}" 
echo "#######################################################################################"
echo "Moving to iteration three"
echo "#######################################################################################"
# aligning output of iteration 1 to the raw reads
minimap2 -ax map-ont ${fname}_iter2.fasta ${rawsPath}/${fname}*.fastq > ${iterVread}/${fname}_readsVsRacon2.sam -t $THREADS
    echo
    echo "Alignment done"
$raconPath $rawsPath/${fname}*.fastq ${iterVread}/${fname}_readsVsRacon2.sam ${fname}_iter2.fasta -t $THREADS > ${fname}_iter3.fasta
# cleaning up. converting and sorting sam file and zipping it
    #convert sam file to bam file
    samtools view --threads $THREADS -bS ${iterVread}/${fname}_readsVsRacon2.sam > ${iterVread}/${fname}_readsVsRacon2.bam
    # sort the bam file
    samtools sort --threads $THREADS -o ${iterVread}/${fname}_readsVsRacon2_sorted.bam ${iterVread}/${fname}_readsVsRacon2.bam
    # clean up
    echo "Deleting intermediate alignment files..."
    rm ${iterVread}/${fname}_readsVsRacon2.sam
    rm ${iterVread}/${fname}_readsVsRacon2.bam
    mv ${iterVread}/${fname}_readsVsRacon2_sorted.bam ${iterVread}/${fname}_readsVsRacon2.bam
    echo
    echo "SAM to BAM conversion done. Sorting of BAM file done."
    echo 
    echo "Iteration three done for ${fname}" 
echo "#######################################################################################"
echo "Moving to iteration four"
echo "#######################################################################################"
# aligning output of iteration 1 to the raw reads
minimap2 -ax map-ont ${fname}_iter3.fasta ${rawsPath}/${fname}*.fastq > ${iterVread}/${fname}_readsVsRacon3.sam -t $THREADS
    echo
    echo "Alignment done"
$raconPath $rawsPath/${fname}*.fastq ${iterVread}/${fname}_readsVsRacon3.sam ${fname}_iter3.fasta -t $THREADS > ${fname}_iter4.fasta
# cleaning up. converting and sorting sam file and zipping it
    #convert sam file to bam file
    samtools view --threads $THREADS -bS ${iterVread}/${fname}_readsVsRacon3.sam > ${iterVread}/${fname}_readsVsRacon3.bam
    # sort the bam file
    samtools sort --threads $THREADS -o ${iterVread}/${fname}_readsVsRacon3_sorted.bam ${iterVread}/${fname}_readsVsRacon3.bam
    # clean up
    echo "Deleting intermediate alignment files..."
    rm ${iterVread}/${fname}_readsVsRacon3.sam
    rm ${iterVread}/${fname}_readsVsRacon3.bam
    mv ${iterVread}/${fname}_readsVsRacon3_sorted.bam ${iterVread}/${fname}_readsVsRacon3.bam
    echo
    echo "SAM to BAM conversion done. Sorting of BAM file done."
    echo 
    echo "Iteration four done for ${fname}"
# get the stats
assembly-stats ${fname}*.fasta >> ${Stats}/${fname}_stats.txt
done
