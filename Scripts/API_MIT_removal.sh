#!/bin/bash
# script to separate the apicoplast and mitochondrial DNA sequences from the Pk genome
# This allows the downstream circularisation of the apicoplast and mitochondrial genomes
# The user can also then carry out annotation of the separated sequences

# Enter the path to the scripts folder
bPATH=/path/to/PipelineScripts
# Enter path to the folder containing the FASTA assemblies to remove API and MIT from
APATH=/path/to/assemblies/to/remove/API/MIT/from/folder
# Enter the path to the folder to save the results
SPATH=/path/to/save/the/outputs
# Enter the path to the reference genome's apicoplast
apiPATH=/path/to/reference/API.fasta
# Enter the path to the reference genome's mitochondrial genome
mitPATH=/path/to/reference/MIT.fasta
## API
mkdir -p $SPATH
mkdir -p $SPATH/BLAST_RefVsAssembly
mkdir -p $SPATH/contigsToRemove
mkdir -p $SPATH/ExtractedContigs
for i in $APATH/*;
do
	gf=$(basename "$i")
	gsr=${gf%.*}
	run="blastn -task megablast -query $apiPATH -subject $i -outfmt 6 -out $SPATH/BLAST_RefVsAssembly/RefAPIvs${gsr}.out"
	echo $run
	eval $run
	cat $SPATH/BLAST_RefVsAssembly/RefAPIvs${gsr}.out | awk '{print $2}' | uniq > $SPATH/contigsToRemove/${gsr}_APIcontigs.txt
done

## MIT
mkdir -p $SPATH/BLAST_RefVsAssembly
mkdir -p $SPATH/contigsToRemove
mkdir -p $SPATH/ExtractedContigs
for i in $APATH/*;
do
	gf=$(basename "$i")
	gsr=${gf%.*}
	run2="blastn -task megablast -query $mitPATH -subject $i -outfmt 6 -out $SPATH/BLAST_RefVsAssembly/RefMITvs${gsr}.out"
	echo $run2
	eval $run2
	cat $SPATH/BLAST_RefVsAssembly/RefMITvs${gsr}.out | awk '{print $2}' | uniq > $SPATH/contigsToRemove/${gsr}_MITcontigs.txt
done
## extract the contigs from the assembly file. Ensuring to delete the identified contigs from the file
source $bPATH/fastaExtract_Delete.sh
extract $SPATH/contigsToRemove/ $APATH $SPATH/ExtractedContigs/
extract $SPATH/contigsToRemove/ $APATH $SPATH/ExtractedContigs/
# delete the contigs from the originals
delete $SPATH/contigsToRemove/ $APATH $SPATH/API_MIT_Free
# get assembly metrics of the extracted contigs
assembly-stats $SPATH/ExtractedContigs/* >> $SPATH/ExtractedContigs_Stats.txt
assembly-stats $SPATH/ExtractedContigs/* >> $SPATH/ExtractedContigs_Stats.txt
