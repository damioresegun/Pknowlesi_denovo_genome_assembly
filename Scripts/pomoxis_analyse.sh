#!/bin/bash
# Script made to carry out different analyses of the assembly after completion. 
# Should hold functions for 1) post raw Assembly 2) post Racon 3) post Medaka 4) post Pilon 5) post repeat masking 6) post RaGOO/RagTag
set -e
input=$1 # path to folder holding directories of isolate assemblies to analyse
output=$2 # output folder to hold the results
refr=$3 # path to the reference file. Note reference file HAS to have a mmi index file generated
thred=$4 # number of threads

# pomoxis for the assembly after cleaning on blobtools
rawPomo() {
	for i in $input/*;
	do
		ls $i
		cpt=$(basename "$i")
		echo $cpt
		mkdir -p $output/postAssembly/${cpt}
		cd $output/postAssembly/${cpt}
		# copy the assembly into the current working directory
		cp $i/assembly.fasta $output/postAssembly/${cpt}/${cpt}.fasta
		assess_assembly -r $refr -i $output/postAssembly/${cpt}/${cpt}.fasta -H -t $thred
		cd
	done
}

# pomoxis for the racon polished genome assembly
raconPomo() {
	for i in $input/*;
	do
		ls $i
		cpt=$(basename "$i")
		echo $cpt
		# check if the correct file is present
		check2r=$i/${cpt}_iter4.fasta
		if [[ -f "$check2r" ]]; then
		mkdir -p $output/postRacon/${cpt}
		cd $output/postRacon/${cpt}
		# copy the assembly into the current working directory
		cp $i/${cpt}_iter4.fasta $output/postRacon/${cpt}/${cpt}.fasta
		assess_assembly -r $refr -i $output/postRacon/${cpt}/${cpt}.fasta -H -t $thred
		cd
		fi
	done
}

#pomoxis for medaka polished genome assembly
medakaPomo() {
	for i in $input/*;
	do
		ls $i
		cpt=$(basename "$i")
		echo $cpt
		mkdir -p $output/postMedaka/${cpt}
		cd $output/postMedaka/${cpt}
		# copy the assembly into the current working directory
		cp $i/${cpt}_consensus.fasta $output/postMedaka/${cpt}/${cpt}.fasta
		assess_assembly -r $refr -i $output/postMedaka/${cpt}/${cpt}.fasta -H -t $thred
		cd
	done
}

#pomoxis for pilon polished genome assembly
pilonPomo() {
	for i in $input/*_iter2/*.fasta;
	do
		ls $i
		cpt=$(basename "$i")
		echo $cpt
		cptr=${cpt%.*}
		echo $cptr
		mkdir -p $output/postPilon/${cptr}
		cd $output/postPilon/${cptr}
		# copy the assembly into the current working directory
		cp $i $output/postPilon/${cptr}
		assess_assembly -r $refr -i $output/postPilon/${cptr}/${cptr}.fasta -H -t $thred
		cd
	done
}

# pomoxis for masked assemblies to be called from this script or command line
maskedPomo() {
	for i in $input/*;
	do
		ls $i
		cpt=$(basename "$i")
		echo $cpt
		cptr=${cpt%.*}
		echo $cptr
		mkdir -p $output/postMasking/${cptr}
		cd $output/postMasking/${cptr}
		# copy the assembly into the current working directory
		cp $i $output/postMasking/${cptr}
		assess_assembly -r $refr -i $output/postMasking/${cptr}/${cptr}.fasta -H -t $thred
		cd
	done
} 

# pomoxis for de-chimerised assemblies
ragtagPomo() {
	for i in $input/*;
	do
		ls $i
		cpt=$(basename "$i")
		echo $cpt
		cptr=${cpt%.*}
		echo $cptr
		mkdir -p $output/postRagTag/${cptr}
		cd $output/postRagTag/${cptr}
		# copy the assembly into the current working directory
		cp $i $output/postRagTag/${cptr}
		assess_assembly -r $refr -i $output/postRagTag/${cptr}/${cptr}.fasta -H -t $thred
		cd
	done
}

# pomoxis for annotated assemblies
annotPomo() {
	for i in $input/*;
	do
		ls $i
		cpt=$(basename "$i")
		echo $cpt
		cptr=${cpt%.*}
		echo $cptr
		mkdir -p $output/postAnnotation/${cptr}
		cd $output/postAnnotation/${cptr}
		# copy the assembly into the current working directory
		cp $i $output/postAnnotation/${cptr}
		assess_assembly -r $refr -i $output/postAnnotation/${cptr}/${cptr}.fasta -H -t $thred
		cd
	done
}
# pomoxis for masked assemblies to be called from a designated script
maskedPomoScript() {
	for i in $input/*;
	do
		ls $i
		cpt=$(basename "$i")
		echo $cpt
		cptr=${cpt%*_*}
		echo $cptr
		mkdir -p $output/${cptr}
		cd $output/${cptr}
		# copy the assembly into the current working directory
		cp $i $output/${cptr}
		assess_assembly -r $refr -i $output/${cptr}/${cptr}*.fasta -H -t $thred
		rm $output/${cptr}/${cptr}*.fasta
		cd
	done
}