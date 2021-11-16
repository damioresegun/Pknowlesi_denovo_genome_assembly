#!/bin/bash
# Script for aligning and sorting SAM files to BAM files for all sequences within the input folder and save to the output folder
# The script aligns an input isolate fastq against the human reference genome, then extracts unmapped reads
# Unmapped reads are then saved and converted to a fastq file and checked for metric statistics
# intermediate files will be deleted i.e. SAM files and unsorted BAM files. Sorted BAM files will be renamed to just .bam rather than _sorted.bam
##NOTE: MIGHT HAVE TO RUN THIS IN A SAMTOOLS CONDA ENVIRONMENT
set -e
#set input and output folders
input=$1	# input isolate fastq to be aligned
outSave=$2	# output directory to hold initial alignment files
fastqSave=$3	# output directory to save generated unmapped fastq file
THREADS=$4	# number of threads
statsFol=$5	# output directory to save statistical assessments
reference=$6	# path to the human reference to check for contimination

mkdir -p $fastqSave
mkdir -p $statsFol/postAlign_Extraction
statsOut=$statsFol/postAlign_Extraction
#make an unmapped and mapped folder
mkdir -p $outSave/Unmapped
unmpd=$outSave/Unmapped
#get the filename
file=$(basename $input)
echo $file
#remove the extension
fname="${file%.*}"
echo $fname
echo
#set the alignment name to the path and the file
alname=${outSave}/${fname}
echo $alname
#align to human reference, convert to bam, and sort
aly="minimap2 -ax map-ont $reference ${input} -t $THREADS | samtools view -@ $THREADS -b - | samtools sort -@ $THREADS -o ${alname}VsHumanRef.bam -"
echo $aly
eval $aly
echo "alignment done"
echo
echo "bam conversion done"
echo ${fname} >> ${statsOut}/${fname}_FlagstatMappedVsHumanRef_stats.txt
samtools flagstat --threads $THREADS ${alname}VsHumanRef.bam >> ${statsOut}/${fname}_FlagstatMappedVsHumanRef_stats.txt
echo "Converted to fastq"
# extract unmapped reads
samtools view --threads $THREADS -f 4 -b ${alname}VsHumanRef.bam > ${unmpd}/${fname}VsHumanRef_unmapped.bam
echo "unmapped extraction done"
#convert unmapped reads to fastq
bedtools bamtofastq -i ${unmpd}/${fname}VsHumanRef_unmapped.bam -fq ${unmpd}/${fname}VsHumanRef_unmapped.fastq
echo
#check the stats
assembly-stats ${unmpd}/${fname}VsHumanRef_unmapped.fastq > ${statsOut}/${fname}_VsHumanRef_FastQ_UnmappedStats.txt
#move fastqs to their own folder
mv ${unmpd}/${fname}VsHumanRef_unmapped.fastq $fastqSave
echo "all done"
