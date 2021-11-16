#!/bin/bash
# script to take your assemblies through for contamination detection using minimap2, samtools, blast and blobtools.
# Requires a text file containing the file paths for the fastq files used to generate assembled isolates
# Will need to change the different commands depending on if Canu or Flye!!!!
set -e
## NOTE: RUN THIS IN THE IN THE BLOBTOOL CONDA ENVIRONMENT!!!!!!!!!!
# $1 = path to blast executable
# $2 = path to blast databases
# $3 = path to txt file containing list of paths for fastq files used to generate assembled isolates
# $4 = path to folder of de novo assemblies
# $5 = path to output blobtools results
# $6 = number of threads

rawBlob() {
# set folder for inputs and outputs
export PATH=${1}:$PATH
export BLASTDB=$2
## inputs
#change these below as appropriate
data=$3
assembly=$4
blobtoole=$5
THREADS=$6
##outputs
alignOutput="${blobtoole}/Raw/Reads_Vs_Assembly"
blastOutput="${blobtoole}/Raw/BLAST_output"
blobOutput="${blobtoole}/Raw/BlobTools_Output"
## make the output folders if not existing
mkdir -p $alignOutput
mkdir -p $blastOutput
mkdir -p $blobOutput

while read -r line;
do
    echo "The isolate you are working on is $line"
    fname=$(basename "$line")    
    echo $fname
    echo
	#get the isolate name
	iso=${fname%.*}
	echo $iso
    echo "Beginning alignment of assembly to reads..."
    echo
    # minimap2 for alignment against the raw reads
	run="minimap2 -ax map-ont ${assembly}/${iso}/assembly.fasta ${line} -t $THREADS | samtools view -@ $THREADS -b - | samtools sort -@ $THREADS -o ${alignOutput}/${iso}_readsVsAssem.bam -"
	echo $run
	eval $run
    echo
    echo "Moving on to Blasting..."
    echo
	blatr="blastn -task megablast -query ${assembly}/${iso}/assembly.fasta -db nt -outfmt '6 qseqid staxids bitscore std scomnames sscinames sblastnames sskingdoms stitle' -evalue 1e-20 -out ${blastOutput}/${iso}_vs_nt.out -num_threads $THREADS"
    echo
	echo $blatr
	eval $blatr
    echo "Blast is finished... Continuing to blobtools"
    echo
    echo "Creating a JSON file of the data"
    mkdir -p ${blobOutput}/${iso}
	echo
	blobcr="blobtools create -i ${assembly}/${iso}/assembly.fasta -b ${alignOutput}/${iso}_readsVsAssem.bam -t ${blastOutput}/${iso}_vs_nt.out -o ${blobOutput}/${iso}/"
    echo $blobcr
	eval $blobcr
	echo
    echo "Creating blobtools view file..."
    blobvi="blobtools view -i ${blobOutput}/${iso}/blobDB.json -o ${blobOutput}/${iso}/"
    echo $blobvi
	eval $blobvi
	echo
    echo "Blobtools table created. Plotting blobplots..."
    echo
    blobpl="blobtools plot -i ${blobOutput}/${iso}/blobDB.json -o ${blobOutput}/${iso}/"
    echo $blobpl
	eval $blobpl
	echo
    echo "Blobplot complete. Drawing covplot..."
    echo
    blobcv="blobtools covplot -i ${blobOutput}/${iso}/blobDB.json -c ${blobOutput}/${iso}/${iso}_readsVsAssem.bam.cov -o ${blobOutput}/${iso}/ --max 1e03"
	echo $blobcv
	eval $blobcv
    echo
    echo "Covplots done"
    done < $data
echo "This iteration is complete. Please look at your results to manually remove contamination. Run this script once contamination has been removed to confirm no more contamination is present"
}

cleanBloby() {
export PATH=${1}:$PATH
export BLASTDB=$2
## inputs
#change these below as appropriate
data=$3
assembly=$4
blobtoole=$5
THREADS=$6
##outputs
alignOutput="${blobtoole}/Clean/Reads_Vs_Assembly"
blastOutput="${blobtoole}/Clean/BLAST_output"
blobOutput="${blobtoole}/Clean/BlobTools_Output"
## make the output folders if not existing
mkdir -p $alignOutput
echo $alignOutput
mkdir -p $blastOutput
mkdir -p $blobOutput
#begin
while read -r line;
do
    echo "The isolate you are working on is $line"
    fname=$(basename "$line")    
    echo $fname
    echo
	#get the isolate name
	iso=${fname%.*}
	echo $iso
    echo "Beginning alignment of assembly to reads..."
    echo
    # minimap2 for alignment against the raw reads
	minimap2 -ax map-ont ${assembly}/${iso}/clean_assembly.fasta ${line} > ${alignOutput}/${iso}_readsVsCleanAssem.sam -t $THREADS
    echo
    echo "Alignment done. Moving on to sam file conversion and sorting..."
    echo
    #convert sam file to bam file
    samtools view --threads $THREADS -bS ${alignOutput}/${iso}_readsVsCleanAssem.sam > ${alignOutput}/${iso}_readsVsCleanAssem.bam
    # sort the bam file
    samtools sort --threads $THREADS -o ${alignOutput}/${iso}_readsVsCleanAssem_sorted.bam ${alignOutput}/${iso}_readsVsCleanAssem.bam
    # clean up
    echo "Deleting intermediate alignment files..."
	# sam files are needed for racon. so they are not deleted here
    rm ${alignOutput}/${iso}_readsVsCleanAssem.bam
    mv ${alignOutput}/${iso}_readsVsCleanAssem_sorted.bam ${alignOutput}/${iso}_readsVsCleanAssem.bam
    echo
    echo
    echo "Moving on to Blasting..."
    echo
	blatr="blastn -task megablast -query ${assembly}/${iso}/clean_assembly.fasta -db nt -outfmt '6 qseqid staxids bitscore std scomnames sscinames sblastnames sskingdoms stitle' -evalue 1e-20 -out ${blastOutput}/${iso}_vs_nt.out -num_threads $THREADS"
    echo
	echo $blatr
	eval $blatr
    echo "Blast is finished... Continuing to blobtools"
    echo
    echo "Creating a JSON file of the data"
    mkdir -p ${blobOutput}/${iso}
	echo
	blobcr="blobtools create -i ${assembly}/${iso}/clean_assembly.fasta -b ${alignOutput}/${iso}_readsVsCleanAssem.bam -t ${blastOutput}/${iso}_vs_nt.out -o ${blobOutput}/${iso}/"
    echo $blobcr
	eval $blobcr
	echo
    echo "Creating blobtools view file..."
    blobvi="blobtools view -i ${blobOutput}/${iso}/blobDB.json -o ${blobOutput}/${iso}/"
    echo $blobvi
	eval $blobvi
	echo
    echo "Blobtools table created. Plotting blobplots..."
    echo
    blobpl="blobtools plot -i ${blobOutput}/${iso}/blobDB.json -o ${blobOutput}/${iso}/"
    echo $blobpl
	eval $blobpl
	echo
    echo "Blobplot complete. Drawing covplot..."
    echo
    blobcv="blobtools covplot -i ${blobOutput}/${iso}/blobDB.json -c ${blobOutput}/${iso}/${iso}_readsVsCleanAssem.bam.cov -o ${blobOutput}/${iso}/ --max 1e03"
	echo $blobcv
	eval $blobcv
    echo
    echo "Covplots done"
    done < $data
echo "Yes This iteration is complete. Please look at your results to manually remove contamination. Run this script once contamination has been removed to confirm no more contamination is present"
}
"$@"