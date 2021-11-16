#!/bin/bash
## Script made to take input data through busco version 5.
## In order to run this, it is necessary to have the busco docker image installed. Version 5.0.0 for busco and version 3.2.2 for Augustus. 
## layout of the script: for every isolate in the flye output, run busco, move the run output folder to the output folder, delete the tmp folder and move to the next isolate

set -e
#set the variables needed
ref1=$1 # path to the reference genome being compared against.
inputy=$2 # path to the directory holding the folder for fasta files to be assessed i.e. path/to/directory/with/folder/containing/isolate.fasta
output=$3 # path to output directory to hold busco results for each isolate
busDock=$4 # path to the busco docker image to be used
lineage=$5 # the lineage to be used in the assessment. For plasmodium, recommended is plasmodium_odb10
THREADS=$6 # number of threads

#set the location of the config file
#export AUGUSTUS_CONFIG_PATH=$augConfig
#echo $AUGUSTUS_CONFIG_PATH
#make the output where necessary
mkdir -p $output
refBUS() {
#run busco for the reference files to see how our data match up
	mkdir -p $output/References
	run="singularity exec -H $(pwd) ${busDock} busco -i $ref1 -l ${lineage} --out MyReference --out_path ${output}/References -f -m geno -c $THREADS --long"
	echo $run
	eval $run
	echo "Busco for your reference is complete"
cd
}
# This function takes in the raw assemblies after they have been generated on Flye
rawIso() {
	input=$inputy/SuccessfulAssemblies
	for f in $input/* ; 
	do
	echo "You are in ${f}" 
	fna=$(basename "$f")
	fname=${fna%*.*}
    echo $fname
	echo "The current directory is: $(pwd)"
		#run busco command
		mkdir -p ${output}/Raw
		PREFIX="$fname"
		echo $PREFIX
		echo "current working directory $(pwd)"
		#flye raw
		run="singularity exec -H $(pwd) ${busDock} busco -i ${f}/assembly.fasta -l ${lineage} --out $PREFIX --out_path ${output}/Raw -f -m geno -c $THREADS --long"
		echo $run
		eval $run
		echo "${PREFIX} done successfully"
	done
}
# This function takes in the output of the assemblies after they have been cleaned with blobtools
cleanIso() {
	input=$inputy/SuccessfulAssemblies
	for f in $input/* ; 
	do
	echo "You are in ${f}" 
	fna=$(basename "$f")
	fname=${fna%*.*}
    echo $fname
	echo "The current directory is: $(pwd)"
	#pwd
		#run busco command
		mkdir -p ${output}/Clean
		PREFIX="$fname"
		echo $PREFIX
		echo "current working directory $(pwd)"
		run="singularity exec -H $(pwd) ${busDock} busco -i ${f}/assembly.fasta -l ${lineage} --out $PREFIX --out_path ${output}/Raw -f -m geno -c $THREADS --long"
		echo $run
		eval $run
		echo "${PREFIX} done successfully"
	done
}
# This functions takes in the assemblies after they have gone through gene prediction and annotation
genoBus() {
input=$inputy/CompanionGenomes
for f in $input/* ; 
	do
	echo "You are in ${f}" 
	#cd $f
	fna=$(basename "$f")
	fname=${fna%*.*}
    echo $fname
	echo "power dirc is:"
	pwd
		#run busco command
		mkdir -p ${output}/AnnotatedGenomes
		PREFIX="$fname"
		echo $PREFIX
		echo "current working directory $(pwd)"	
		#####GENOME 
		run="singularity exec -H $(pwd) ${busDock} busco -i $f -l ${lineage} --out $PREFIX --out_path ${output}/AnnotatedGenomes -f -m geno -c $THREADS --long"
		echo $run
		eval $run
		echo "${PREFIX} done successfully"		
done
}
# This functions takes in the assemblies after they have gone through medaka
medyBus() {
input=$inputy/Medaka_Output
for f in $input/* ; 
	do
	echo "You are in ${f}" 
	#cd $f
	fname=$(basename "$f")
    echo $fname
	echo "power dirc is:"
	pwd
		#run busco command
		mkdir -p ${output}/Medaka		
		PREFIX="$fname"
		echo $PREFIX
		echo "current working directory $(pwd)"	
		#####MEDAKA OUTPUT 
		run="singularity exec -H $(pwd) ${busDock} busco -i ${f}/${fname}_consensus.fasta -l ${lineage} --out $PREFIX --out_path ${output}/Medaka -f -m geno -c $THREADS --long"
		echo $run
		eval $run
		echo "${PREFIX} done successfully"		
done
}
# This function is for busco of pilon outputs
pilBus() {
input=$inputy/Pilon_Output
for f in $input/*_iter2 ; 
	do
	echo "You are in ${f}" 
	fna=$(basename "$f")
	fname=${fna%*_*}
    echo $fname
	echo "power dirc is:"
	pwd
		#run busco command
		mkdir -p ${output}/Pilon		
		PREFIX="$fname"
		echo $PREFIX
		echo "current working directory $(pwd)"	
		#####MEDAKA OUTPUT 
		run="singularity exec -H $(pwd) ${busDock} busco -i ${f}/${fname}.fasta -l ${lineage} --out $PREFIX --out_path ${output}/Pilon -f -m geno -c $THREADS --long"
		echo $run
		eval $run
		echo "${PREFIX} done successfully"		
done
}
# This function is for busco of masked assemblies that is command line called
maskBus() {
input=$inputy/MaskedAssemblies
for f in $input/* ; 
	do
	echo "You are in ${f}" 
	fna=$(basename "$f")
	fname=${fna%*_*}
    echo $fname
	echo "power dirc is:"
	pwd
		#run busco command
		mkdir -p ${output}/MaskedGenomes	
		PREFIX="$fname"
		echo $PREFIX
		echo "current working directory $(pwd)"	
		#####MEDAKA OUTPUT 
		run="singularity exec -H $(pwd) ${busDock} busco -i $f -l ${lineage} --out $PREFIX --out_path ${output}/MaskedGenomes -f -m geno -c $THREADS --long"
		echo $run
		eval $run
		echo "${PREFIX} done successfully"		
done
}
# This function is for busco of masked assemblies for masked assemblies quality assessment script
maskBusScript() {
input=$inputy
for f in $input/* ; 
	do
	echo "You are in ${f}" 
	fna=$(basename "$f")
	fname=${fna%*_*}
    echo $fname
	#echo "Current directory is:"
	pwd
		#run busco command
		mkdir -p ${output}/BUSCO	
		PREFIX="$fname"
		echo $PREFIX
		echo "current working directory $(pwd)"	
		run="singularity exec -H $(pwd) ${busDock} busco -i $f -l ${lineage} --out $PREFIX --out_path ${output}/BUSCO -f -m geno -c $THREADS --long"
		echo $run
		eval $run
		echo "${PREFIX} done successfully"		
done
}
