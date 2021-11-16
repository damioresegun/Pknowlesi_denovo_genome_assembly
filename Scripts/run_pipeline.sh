#!/bin/bash
set -e
#############################################################################################################################################################################
#													Automation script for nanopore sequencing data analysis																	#
#																		Version 0.1.0																						#
# 								Designed to take in raw nanopore data for basecalling, demultiplexing, adapter removal, alignment and de novo assembly.						#
# 											This script works as a wrapper script that calls multiple smaller scripts for different steps. 									#
#											Hence this script is quite modular and can be modified further to take in other scripts.										#
#										It is assumed that there is some human contamination in the raw reads being processed												#
#		 Requirements: 																																						#
#			 1: The called scripts need to be in the $SCPTS variable.																										#
#			 2: A GPU is required for certain aspects of this pipeline. A GPU is needed for the basecalling and demultiplexing steps. 										#
#			 3: pigz needs to be in the $PATH																																#
# 			 4: assembly-stats is required in the $PATH																														#
#			 5: samtools is required in the the $PATH																														#
#			 6: path to busco docker is required to be provided. Note: This script was developed using busco v3.0.0															#
#			 7: path to augustus version 3.2.2 need to be provided																											#
#			 8: QUAST needs to in the $PATH																																	#
#			 9: RACON is needed in the $PATH																																#
#			 10: Porechop was downloaded from the associated link: https://github.com/rrwick/Porechop																		#
#																																											#
# 		NOTE: It is important to note that if a step fails or needs to be re-run, already complete steps will need to be commented out. 									#
# 			  This wrapper script is currently not dynamic enough to determine if a step has been completed.																#
#																																											#
#																																											#
#															PIPELINE OUTLINE																								#
#		Basecalling (Guppy) -> Demultiplexing (Qcat) -> Adapter Removal (Porechop) -> Align Vs Ref (minimap2) -> Extract FASTQs -> 	De Novo assembly (FLYE) ->				#
#																																											#
#			Assemblies -> Decontimatination (Blobtools) -> Quality Assessment (BUSCO/QUAST) -> Pilon -> Medaka -> Racon (x4) -> Quality Assessment (BUSCO/QUAST/FRCbam)	    #
#																																											#
#				-> Remove MIT/API -> Circularise MIT/API -> Mask repeats on MIT/API-free Assemblies (RepeatModeler/Masker/CD-HIT/Transposon PSI/OCTFTA)  					#
#																																											#
#				-> Masked Assemblies -> Remove chimeric contigs -> Statistics -> Ready for Manual/Automated Annotation														#
#	Output Folders:																																							#
#		1. Basecalled: Holds the output of the basecalling step. Raw reads still with adapters and not demultiplexed														#
#		2. Demultiplexed: Output of the demultiplexing step. Reads have been separated into their barcoded bins. Still possessing adapters									#
#		3. AdapterRemoved: Output of the adapter removal step based on the experiment name provided by user. 																#
#		4. IsolatesToAlign: Containing adapter removed fastqs that have been renamed based on the isolate IDs provided by the user. Reads are now ready for alignment		#
#		5. ReadsVsRefAlignments: Contains bam alignment files of the reads against the provided reference genome															#
#			5a. Mapped: extracted alignments of only the mapped regions of the reads vs reference																			#
#			5b. Unmapped: extracted alignments and fastq of the unmapped regions of the reads vs reference																	#
#		6. Mapped/UnmappedIsolateFastQs: Extracted fastqs of the mapped/unmapped reads from 5a depending on user's entry													#
#		7. IsolatesToAssemble: Formatted fastqs of the extracted reads from (6). Ready to assemble 																			#
#		8. DeNoVoAssembly: containing assemblies based on the isolate ID given by user																						#
#		9. PreviousRunData: Containing previously run data with the same isolate ID as a current run. Data has been chosen to be kept by user								#
#	   10. Jellyfish: holding output of jellyfish for the isolates provided with Illumina sequencing. Outputs here to be uploaded onto GenomeScope							#
#	   11. InputFiles: Holding textfiles for user chosen isolates to take forward for downstream analysis after de novo assembly											#
#	   12. BlobTools: Containing the process of decontaminating the user chosen assemblies.																					#
#	   		12a. 																																							#
#																																											#
#																																											#
#																																											#
#	   13. Stats: containing various statistical results for different steps of the pipeline																				#
#	   14. LogFiles: contains logfiles created as the pipeline proceeds																										#
#	   15. Duplicate_Isolates_AdapterRemoved: folder to hold adapter removed fasta files of duplicate isolates. Useful if you wish to investigate the yield of each barcode	#
#	   16. PreviousRunData: folder holding any previously run data detected and chosen to be saved by the user																#
#	   17. SuccessfulAssemblies: Folder to hold FASTA files of de novo assemblies the user chose to carry forward															#
#############################################################################################################################################################################
#.......................................................
#					SET PATHS
#  These are paths that you will be required to change to suit your particular paths
#.......................................................
# Path to pipeline scripts
SCPTS=/path/to/where/you/saved/package/Scripts
# Path to raw sequencing reads
RPATH=	# leave blank if you have already basecalled your data
# Is the data already basecalled? Yes/Y/No/N
BSECL="yes"
# Path to the MinION Basecalled Data to be analysed
DPATH= # leave blank if you started with raw reads
# If you have access to the Illumina Sequencing data for the same isolates. Provide the path to those here.
ILLUMPATH= # leave blank if unneccesary
# Path to save output. There will be multiple outputs per modular script so this should be a path to hold the outputs of each of the isolates
SPATH= # the folder to hold outputs
# Where is your reference genome. If you have indexed it with minimap; provide the .mmi path
REFPATH=installationPath/Index/GCA_000001405.28_GRCh38.p13_genomic.fna # be aware that you do need the fasta file to carry out busco
REFPATHR=installationPath/Index/GCA_000001405.28_GRCh38.p13_genomic.fna.mmi
PKREF=installationPath/Index/Pknowlesi_contigs.fasta # pknowlesi reference
# Where is the gff file for your reference genome
REFGFF=installationPath/Index/GCA_000001405.28_GRCh38.p13_genomic.gff # HUMAN GFF
PKGFF=installationPath/Index/Pknowlesi.gff # PKNOWLESI GFF
# What is the size of your genome to the nearest megabase/gigabase e.g. 24m or 24g
GENOME="25m" # can be changed to your liking
# What is the name of the experiment
EXPMT="Exp07_Cultured_ClinicalIsolates"
# Which barcodes were used
BCODES=("barcode01" "barcode02" "barcode03" "barcode04" "barcode05")
# Set the isolate names -- correspond this with the barcodes they were used with!!! i.e. barcode01 was Cultured
ISOLATES=("Cultured" "Cultured" "Patient1" "Patient2" "Patient3" "Patient3") # can be named anything you wish
# Is an isolate duplicated? This is necessary if a single isolate was used for multiple barcodes
DUPL="Y" # or N/No/no
# Which isolate(s) is duplicated
DUPPY=("Cultured" "Patient3")	# separate isolates with space
# Enter isolates you have Illumina sequences for
ILLUMISO=("Cultured" "Patient1" "Patient2" "Patient3")
# Which protocol was used for the sequencing run - enter both
KIT="SQK-RBK004" 	# for guppy
KIT2="RBK004"		# for qcat. Note: If you are using PBK004 for your kit, please enter PBK004/LWB001. Qcat documentation for further information on convention
# Which flowcell was used for sequencing
FLOWCELL="FLO-MIN106" # may have used FLO-MIN110 or FLO-MIN111
# Number of threads
THREADS="12" # no maximum thread count
# Enter a kmer length to use for processing your files. Must be an odd number
KMER="27" # set a kmer number
# Enter the species to use for your quality assessments in BUSCO, QUAST and AUGUSTUS. Be aware, that your organism may not have a model
# and you will have to enter the closest organism that has been modelled.
SPEC="pfalciparum" # pfalciparum is best used for pknowlesi though pipeline may also be suitable for other organisms


#........................................................
#				SET PATHS FOR PROGRAMS
#			Change these to suit you
#........................................................
# path to guppy bin folder
ONT="/path/to/guppy/bin"
#path to porechop.py 
PORCP="/path/to/Porechop/porechop-runner.py"
# path to blast bin
blastPATH="/path/to/ncbi-blast-2.7.1+/bin/"
#path to blast databases
blastDATA="/path/to/blastntnr/blastDatabases"
# path to the folder that holds the augustus config
augConfig="installationPath/Tools/augustusConfig/config/"
# path to the busco docker image
busDock="installationPath/Tools/BUSCOv5/busco_v5.0.0_cv1.sif"
# path to the lineage library you wish to use
LINEAGE="plasmodium_ensembl"

##############################################################################################################################################################
# DO NOT CHANGE BELOW THIS. SCRIPT IS READY TO RUN
##############################################################################################################################################################
#........................................................
#			CHECKS
# do not change these
#........................................................

  # Loop forever until the user enters a valid response (Y/N or Yes/No).
while true; do
	read -r -p 'There will more prompts for you to confirm at different stages of this script. Please be aware of this. Proceed? Y/N: ' starting
    case "$starting" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo "Beginning the script..."
		break
		;;
      [Nn][Oo]|[Nn])  # No or N.
        #return 1
		echo "The script will end here. Thank you."
		exit
        ;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
  done

# confirming the paths
echo "Your raw files are in ${RPATH}"
echo "Your basecalled files are/will be saved in ${DPATH}"
echo "Your outputs will be saved in ${SPATH}/${EXPMT}"
echo "Your barcodes are ${BCODES[*]}"
echo "Your isolates are ${ISOLATES[*]}"
if [[ $BSECL = "Yes" ]] || [[ $BSECL = "y" ]] || [[ $BSECL = "YES" ]] || [[ $BSECL = "yes" ]] || [[ $BSECL = "Y" ]]; then
	echo "Your data is basecalled, so pipeline will not basecall. Next step will be demultiplexing"
else
	echo "Your data is not basecalled, so pipeline will begin basecalling"
fi
if [[ $DUPL = "Yes" ]] || [[ $DUPL = "y" ]] || [[ $DUPL = "YES" ]] || [[ $DUPL = "yes" ]] || [[ $DUPL = "Y" ]]; then
	echo "Your data is duplicated. The duplicated isolates are ${DUPPY[*]}"
	
else
	echo "Your data is not duplicated"
fi

echo "The flowcell you used is $FLOWCELL"

# Prior checks for save destination
if [ -d "$SPATH" ]; then
	echo "${SPATH} exists. Subsequent folders will be made here"
else
	echo "{SPATH} does not exist. The folder and parent folders in the path will be made to hold subsequent folders and data"
mkdir -p $SPATH
echo "The save folder and parent folders [${SPATH}] have been created"
fi
echo


#......................................................
#				Starting
#......................................................
WPATH=$SPATH
cd $WPATH
echo "Your working directory is currently ${WPATH}"
echo
#......................................................
#				BASECALLING
#......................................................
# Make directory for basecalling outputs

# Carry out basecalling if not previously done. Call the basecalling script
if [[ $BSECL = "No" ]] || [[ $BSECL = "n" ]] || [[ $BSECL = "NO" ]] || [[ $BSECL = "no" ]] || [[ $BSECL = "N" ]]; then
	echo "You said your files are not basecalled. So basecalling will be carried out"
	echo "Beginning basecalling..."
	mkdir -p $WPATH/Basecalled
	DPATH=$WPATH/Basecalled
	${SCPTS}/Basecalling.sh $RPATH $DPATH $FLOWCELL $KIT $ONT
	echo "Base calling done. Proceeding to Demultiplexing..."
else 
		echo "Basecalling has been previously done. Proceeding to Demultiplexing..."
fi
#................................................................
#
#................................................................
#						DEMULTIPLEXING
#................................................................
# make the demultiplex output directory
mkdir -p $WPATH/Demultiplexed/${EXPMT}
DEMULP=$WPATH/Demultiplexed/${EXPMT}
LOGFF=$WPATH/LogFiles
mkdir -p $LOGFF
mkdir -p $LOGFF/Demultiplexing
# call demultiplex script
${SCPTS}/Demultiplexing.sh $DPATH $DEMULP $THREADS $KIT2 $ONT 2>&1 | tee $LOGFF/Demultiplexing/${EXPMT}.txt
echo "Demultiplexing complete. Moving on to adapter removal..."
#
#................................................................
#						ADAPTER REMOVAL
#................................................................
echo "Beginning Adapter Removal..."
APTRM=$WPATH/AdapterRemoved/${EXPMT}
mkdir -p $LOGFF/AdapterRemoval
mkdir -p $APTRM
${SCPTS}/AdapterRemoval.sh $PORCP $DEMULP $EXPMT $APTRM $THREADS 2>&1 | tee $LOGFF/AdapterRemoval/${EXPMT}.txt
#................................................................
#			RENAME BARCODES TO ISOLATES and GET SOME STATS
#................................................................
#make output folder
REISOL=$WPATH/IsolatesToAlign
mkdir -p $REISOL
# some isolates may be repeated and you may want to know the coverage for each iteration of the repeated isolates. So these will be copied into a different folder
CPYF=$WPATH/Duplicate_Isolates_AdapterRemoved
mkdir -p $CPYF
# generate a random number to append to file name to make it unique
nnmm=$RANDOM
# generate the date the adapter removal is done, in case an isolate is repeated again
ttdday=$(date +"%d_%m_%Y")
# set the iteration number to 0
n=0
# for each barcode stated, take the corresponding isolate in the ISOLATES variable and rename
echo
for i in ${BCODES[*]}; 
do
echo "The barcoded file to be renamed is $i...";
echo "It will be renamed ${ISOLATES[$n]}.fastq";
FILE="$REISOL/${ISOLATES[$n]}.fastq"
# check if the isolate already has a previous file
	if [[ -f "$FILE" ]]; then
		echo "The file for isolate ${ISOLATES[$n]}.fastq exists. Will now concatenate to this file"
		cp $APTRM/$i.* $CPYF/${ISOLATES[$n]}_${nnmm}_${ttdday}.fastq.gz
		zcat $APTRM/$i.* >> $FILE;
		echo "The data has been concatenated and saved in ${FILE}..."
		echo "A copy of the new addition has been saved in ${CPYF}/${ISOLATES[$n]}_${nnmm}_${ttdday}.fastq.gz"
	else
	# if file doesn't exist then write to file.
		zcat $APTRM/$i.* > $FILE
		echo "${ISOLATES[$n]} has been written to ${FILE}"
		echo
	fi
# add one to the iteration
n=`expr $n + 1`;

# do some stats
mkdir -p $WPATH/Stats
STTS=$WPATH/Stats
mkdir -p $WPATH/Stats/AdapterRemovedData
TOALI=$WPATH/Stats/AdapterRemovedData
assembly-stats $FILE >> $TOALI/${EXPMT}_stats_${ttdday}.txt
echo "The stats have been saved in $TOALI/${EXPMT}_stats_${ttdday}.txt"
done
#................................................................
#					ALIGNMENT AGAINST THE REFERENCE
#................................................................
echo
echo "Starting alignment against the reference provided"
echo "PLEASE NOTE: Alignment will be done based on the isolates you have provided" 
echo "If you have already sequenced the same isolate previously, the new sequence data will be added to the previous one"
echo "Then the concatenated sequence will be re-aligned against the reference"
echo "Furthermore, any files previously generated from the isolate data will also be deleted"
echo "If this is a problem, stop the script here. If not, answer the prompt as necessary..."
echo
echo "Do you want the alignment to proceed? Be aware that previous isolates with the same name can be overwritten! or you can choose to move them to another folder Proceed? Y/N"
read -p 'Proceed with alignment: ' contProc
if [[ $contProc = "No" ]] || [[ $contProc = "n" ]] || [[ $contProc = "NO" ]] || [[ $contProc = "no" ]] || [[ $contProc = "N" ]]; then
	echo "The script will end here. It may be necessary to rename your new isolate and then re-run this script."
	echo
	echo "To save time for the re-run, consider commenting the demultiplexing and adapter removal stages."
	echo
	echo "Script will end. Please come again"
	exit
else
	echo "The script will continue. Any previous duplicate isolate analysis will be overwritten"
fi

#make folder to hold the alignments of the reads vs the human reference
ALIPATH=$WPATH/ReadsVsHumanRef
mkdir -p $ALIPATH
echo
echo "The alignment of the isolates of your experiment will be saved in ${ALIPATH}"
# make folder to hold the extracted fastqs of unmapped reads to be taken forward for de novo assembly
ALIFAST=$WPATH/UnMappedIsolateFastQs
mkdir -p $ALIFAST
echo "Once alignments are complete, the fastq of unmapped reads will be saved in ${ALIFAST}"
ALIFORM=$WPATH/IsolatesToAssemble
mkdir -p $ALIFORM
echo "The unmapped FASTQs will be formatted and then saved in ${ALIFORM}..."
mkdir -p $STTS/UnMappedExtractedFastQs
mkdir -p $LOGFF/Alignment
# in the case there are duplicate isolates indicated, unique isolates names are taken forward as the data will already be concatenated
ISOLATTE="$( echo "${ISOLATES[@]}" | xargs -n1 | sort -u | xargs )"
ISOLATES=${ISOLATTE}
echo "The unique isolates are ${ISOLATES}"
##### can manually change the next few lines
# if you want to align all isolates in the IsolatesToAlign folder then uncomment the next 2 lines (and comment the third and fourth line)
ISOLATES=$REISOL	# uncomment this if you are doing the entire IsolatesToAlign folder 
for k in $ISOLATES/*;	# uncomment this if you are doing the entire IsolatesToAlign folder
#for i in ${ISOLATES[*]};
#for i in ${ISOLATES}; # uncomment if you are aligning just a particular experiment/set of isolates
do
vr=$(basename "$k") # uncomment this if you are doing the entire IsolatesToAlign folder
i=${vr%.*}	# uncomment this if you are doing the entire IsolatesToAlign folder
echo $i
check=$ALIPATH/${i}VsHumanRef.bam
echo $check
	if [[ -f "$check" ]]; then
	echo
	lasmood=$(date -r ${check} +"%d_%m_%Y")
		echo "The isolate ${i} already exists. This will be overwritten!!!"
		echo "Alternatively, you can choose to save these files in another directory"
		 read -p 'Do you want to delete? Yes = delete; No = save in new location; Skip = skip this isolate, move to the next: ' newSave
		if [[ $newSave == "No" ]] || [[ $newSave = "n" ]] || [[ $newSave = "NO" ]] || [[ $newSave = "no" ]] || [[ $newSave = "N" ]]; then
			mkdir -p $WPATH/PreviousRunData/${i}_${lasmood}
			echo "The previously generated files will be moved to ${WPATH}/PreviousRunData/${i}..."
			mv $check $WPATH/PreviousRunData/${i}_${lasmood}
			mv $STTS/postAlign_Extraction/${i}_FlagstatMappedVsHumanRef_stats.txt $WPATH/PreviousRunData/${i}_${lasmood}
			mv $ALIPATH/Unmapped/${i}VsHumanRef_unmapped.bam $WPATH/PreviousRunData/${i}_${lasmood}
			mv $STTS/postAlign_Extraction/${i}_VsHumanRef_FastQ_UnmappedStats.txt $WPATH/PreviousRunData/${i}_${lasmood}
			mv $ALIFAST/${i}VsHumanRef_unmapped.fastq $WPATH/PreviousRunData/${i}_${lasmood}
			mv $ALIFORM/${i}* $WPATH/PreviousRunData/${i}_${lasmood}
			mv $STTS/UnMappedExtractedFastQs/${i}_UnMappedReads.txt $WPATH/PreviousRunData/${i}_${lasmood}
			mv $LOGFF/Alignment/${i}* $WPATH/PreviousRunData/${i}_${lasmood}
			echo
			echo "Your files have been moved. The new run will now take place"
		elif [[ $newSave == "Yes" ]] || [[ $newSave = "y" ]] || [[ $newSave = "YES" ]] || [[ $newSave = "yes" ]] || [[ $newSave = "Y" ]]; then
			rm $check
			rm $STTS/postAlign_Extraction/${i}*
			rm $ALIPATH/Unmapped/${i}VsHumanRef_unmapped.bam
			rm $ALIFAST/${i}VsHumanRef_unmapped.fastq
			rm $ALIFORM/${i}*
			rm $STTS/UnMappedExtractedFastQs/${i}_UnMappedReads.txt
			rm $LOGFF/Alignment/${i}*
			echo "Files have been removed"
			echo "New run will now take place"
		elif [ $newSave == "Skip" ]] || [[ $newSave = "s" ]] || [[ $newSave = "SKIP" ]] || [[ $newSave = "skip" ]] || [[ $newSave = "S" ]]; then
		echo "The isolate ${i} will be skipped"
		continue
		fi
	fi
echo
echo "Running alignment.."
align="${SCPTS}/Alignment.sh ${REISOL}/${i}.fastq $ALIPATH $ALIFAST $THREADS $STTS $REFPATHR 2>&1 | tee $LOGFF/Alignment/${i}_alignments.txt" 
echo $align
eval $align
echo "Alignment for ${i}.fastq done"
# format the fastq file
formatr="${SCPTS}/fastqformatter.sh ${SCPTS} ${WPATH}/UnMappedIsolateFastQs/${i}* ${ALIFORM} ${i} 2>&1 | tee $LOGFF/Alignment/${i}_formatting.txt" 
eval $formatr
echo "Formatting for ${i}.fastq complete..."
pigz --best ${WPATH}/MappedIsolateFastQs/${i}*.fastq
pigz --best $ALIPATH/Mapped/${i}*.fastq
formsta="assembly-stats $ALIFORM/${i}* >> $STTS/UnMappedExtractedFastQs/${i}_UnMappedReads.txt"
echo $formsta
eval $formsta
done
# Alignment stage complete
echo
echo "Alignments, conversion and formatting complete"
echo
echo "Alignments can be found in $ALIPATH. Mapped and Unmapped alignment files are held within"
echo 
echo "Converted and formatted fastq files are in $ALIFORM. These are ready to be assembled de novo"
echo
echo "Next step is de novo assembly with Flye"
#................................................................
#					DE NOVO ASSEMBLY USING FLYE
#................................................................
echo "Do you want to use the barcodes of this ${EXPMT} experiment [1] or use all the isolates in the folder [2]?"
  # Loop forever until the user enters a valid response (Y/N or Yes/No).
  Assems=$WPATH/DeNoVoAssembly
  mkdir -p $Assems
  mkdir -p $LOGFF/Assembly
while true; do
	read -r -p 'Enter 1 or 2: ' assemStart
    case "$assemStart" in
      [1]) # Yes or Y (case-insensitive).
        #return 0
		echo "Beginning the script..."
		echo
		whichr="Some"
		echo "You will be assembling ${ISOLATES} on Flye"
		break
		;;
      [2])  # No or N.
        #return 1
		whichr="All"
		echo "Assembling all isolates"
		ALIFORM=$WPATH/IsolatesToAssemble
		ALICONT="$(ls ${ALIFORM})"
		ISOLATES=$ALIFORM
		echo $whichr
		break
        ;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
  done
  echo
  echo
### if all assemblies are chosen then do this:
if [[ $whichr == "Some" ]]; then
### if only experiment isolates are chosen then run this:
	echo "Only isolates in the ${EXPMT} dataset will be assembled"
	#for i in ${ISOLATES[*]};
	for i in ${ISOLATES};
		do
		echo $i
		check=$Assems/${i}/assembly.fasta
		echo $check
		if [[ -f "$check" ]]; then
			echo
			echo "The isolate ${i} already exists. This will be overwritten!!!"
			echo "Alternatively, you can choose to save these files in another directory"
			read -p 'Proceed to delete? Yes = delete; No = save in new location; Skip = skip this isolate, move to the next: ' newSave
			if [[ $newSave == "No" ]] || [[ $newSave = "n" ]] || [[ $newSave = "NO" ]] || [[ $newSave = "no" ]] || [[ $newSave = "N" ]]; then
				mkdir -p $WPATH/PreviousRunData/${i}/assembly
				echo "The previously generated files will be moved to ${WPATH}/PreviousRunData/${i}..."
				mv $Assems/${i}/* $WPATH/PreviousRunData/${i}/assembly
				echo
				echo "Your files have been moved. The new run will now take place"
			elif [[ $newSave == "Yes" ]] || [[ $newSave = "y" ]] || [[ $newSave = "YES" ]] || [[ $newSave = "yes" ]] || [[ $newSave = "Y" ]]; then
				rm -rf $Assems/${i}/*
				echo "Files have been removed"
				echo "New run will now take place"
			elif [ $newSave == "Skip" ]] || [[ $newSave = "s" ]] || [[ $newSave = "SKIP" ]] || [[ $newSave = "skip" ]] || [[ $newSave = "S" ]]; then
				echo "The isolate ${i} will be skipped"
				continue
			fi
		echo "${SCPTS}/Assembly.sh $ALIFORM/${i}* $Assems $THREADS $i $GENOME"
		fi
		${SCPTS}/Assembly.sh $ALIFORM/${i}* $Assems $THREADS $i $GENOME  2>&1 | tee $LOGFF/Assembly/${i}_assembly.txt
	done
	## carry out stats on the assembly
	${SCPTS}/batchAssemStats.sh $Assems $STTS/Assembly $whichr $ISOLATES
	
elif [[ $whichr == "All" ]]; then
	echo "All isolates in the $ALIFORM folder will be assembled..."
	for i in $ISOLATES/*;
		do
		echo $i
		isola=$(basename "$i")
		isol=${isola%.*}
		check=$Assems/${isol}/assembly.fasta
		echo $check
			if [[ -f "$check" ]]; then
			echo
				echo "The isolate ${isol} already exists. This will be overwritten!!!"
				echo "Alternatively, you can choose to save these files in another directory"
				read -p 'Proceed to delete? Yes = delete; No = save in new location; Skip = skip this isolate, move to the next: ' newSave
				if [[ $newSave == "No" ]] || [[ $newSave = "n" ]] || [[ $newSave = "NO" ]] || [[ $newSave = "no" ]] || [[ $newSave = "N" ]]; then
					mkdir -p $WPATH/PreviousRunData/${isol}/assembly
					echo "The previously generated files will be moved to ${WPATH}/PreviousRunData/${isol}..."
					mv $Assems/${isol}/* $WPATH/PreviousRunData/${isol}/assembly
					echo
					echo "Your files have been moved. The new run will now take place"
				elif [[ $newSave == "Yes" ]] || [[ $newSave = "y" ]] || [[ $newSave = "YES" ]] || [[ $newSave = "yes" ]] || [[ $newSave = "Y" ]]; then
					rm -rf $Assems/${isol}/*
					echo "Files have been removed"
					echo "New run will now take place"
				elif [[ $newSave == "Skip" ]] || [[ $newSave = "s" ]] || [[ $newSave = "SKIP" ]] || [[ $newSave = "skip" ]] || [[ $newSave = "S" ]]; then
					echo "The isolate ${isol} will be skipped"
					continue
				fi
			fi
	echo $isol
	${SCPTS}/Assembly.sh $i $Assems $THREADS $isol $GENOME 2>&1 | tee $LOGFF/Assembly/${isol}_assembly.txt
	done
	${SCPTS}/batchAssemStats.sh $Assems $STTS/Assembly $whichr $ISOLATES
fi
echo "De Novo Assembly completed."
echo "However not all isolates may have assembled. This may be due to errors during the assembly"
echo "More likely, the isolates may have shorter than expected lengths which may be due to the input coverage of your data"
echo "To troubleshoot errors in the process, have a look at the log file produced when you ran the script"
echo
echo
echo
echo "For assemblies that did not proceed as intended, you can choose to remove these from the analysis"
echo "To determine which isolates are unsuitable, if you have access to some Illumina sequences for these isolates, please carry out a jellyfish analysis after this"
echo "Following the jellyfish, upload the outputs to genomescope to determine the estimated expected length"
echo "DO NOT STOP THIS SCRIPT OR YOU WILL HAVE TO START AGAIN"
echo
echo
echo
#................................................................
#					JELLYFISH FOR ESTIMATED GENOME SIZE
#................................................................
# user input required 
echo "CHECK THAT YOU HAVE SAVED CORRECT PATHS IN THE INPUTFILES FILE!!!"
echo "Would you like to carry out the jellyfish run? Y/N." 
echo "Enter 'No' if you do not have access to Illumina sequences. This step will be skipped"
while true; do
	read -r -p 'Continue with Jellyfish? [Y/N]: ' jello
    case "$jello" in
      [yY][eE][sS]|[yY]) # Yes or Y (case-insensitive).
        #return 0
		echo "Beginning jellyfish..."
		jelly=$WPATH/Jellyfish
		mkdir $jelly
		mkdir -p $LOGFF/Jellyfish
		for i in ${ILLUMISO[*]};
		do
		mkdir -p ${jelly}/${i}
		jellyfish count -C -m $KMER -s 10000000000 -t $THREADS $ILLUMPATH/${i}/*.fastq -o $jelly/${i}/reads.jf 2>&1 | tee $LOGFF/Jellyfish/${i}_count.txt
		jellyfish histo -t $THREADS $jelly/${i}/reads.jf > $jelly/${i}/reads.histo 2>&1 | tee $LOGFF/Jellyfish/${i}_histo.txt
		done
		echo "Please upload ${jelly}/${i}/reads.histo the results to genomescope: http://qb.cshl.edu/genomescope/"
		echo "Jellyfish done"
		break
		;;
      [nN][oO]|[nN])  # No or N.
        #return 1
		echo "Skipping the jellyfish stage..."
		break
        ;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
  done
 
jells=$WPATH/InputFiles
mkdir -p $jells
touch $jells/AssembliesToProcess.txt

echo "Once you have completed the jellyfish and genomescope, choose which isolates meet your criteria"
echo "Then enter the FULL PATHS of the ASSEMBLY INPUT FASTQ FILE and save these in $jells/AssembliesToProcess.txt"
echo "The FASTQ files can be found in ${ALIFORM}"
echo "If no Illumina sequences are available then just enter the paths to the isolates you wish to keep in the text file above"

while true; do
	read -r -p 'Have you entered the paths to the fastq files of successful assemblies you want to process in ${jells}/AssembliesToProcess.txt? [Y/N]: ' ctnno
    case "$ctnno" in
      [yY][eE][sS]|[yY]) # Yes or Y (case-insensitive).
	  echo "Continuing to next step..."
	  break
	  ;;
	  [nN][oO]|[nN])  # No or N.
	  echo "You have said you have not entered the paths"
	  echo "Please enter the path of the fastq assemblies you want to take forward in the $jells/AssembliesToProcess.txt file"
	  echo "The paths to these fastq files can be found in ${ALIFORM}"
	  echo "To know which assemblies were successful, the stats file can be found in ${STTS}/Assembly"
	  read -n 1 -s -r -p "Once you have done this, press any key to continue: "
	  break
	  ;;
	  *)
	  ;;
	  esac
  done
#................................................................
#					BLOBTOOLS DECONTAMINATION
#................................................................
#make the output folders
blobOut=$WPATH/BlobTools
mkdir -p $blobOut
mkdir -p $LOGFF/BlobTools
# call blobtools
echo "Begining blobtools decontamination..."
. ${SCPTS}/Blobtools.sh 
runr="rawBlob $blastPATH $blastDATA ${jells}/AssembliesToProcess.txt $Assems $blobOut $THREADS 2>&1 | tee $LOGFF/BlobTools/iteration1.txt"
echo $runr
eval $runr
echo
echo "BlobTools iteration one done"
# Decontaminating the reads
echo
echo "Decontaminating your assembly"
${SCPTS}/Decontamination.sh $blobOut/Raw/BlobTools_Output $Assems $SCPTS 2>&1 | tee $LOGFF/BlobTools/iteration1_decontamination.txt
echo
fnAssem=${WPATH}/SuccessfulAssemblies
echo "Carrying out BlobTools on the clean assembly to ensure all unwanted sequences have been removed"
. ${SCPTS}/Blobtools.sh 
runng="cleanBloby $blastPATH $blastDATA $jells/AssembliesToProcess.txt $Assems $blobOut $THREADS 2>&1 | tee $LOGFF/BlobTools/iteration2_clean.txt"
echo $runng
eval $runng
echo
echo "Blobtools for the clean assemblies completed. You may want to check the results to ensure that these sequences have been removed before proceeding"
# check if user wants to proceed
while true; do
	read -r -p ' Have you checked your Blobtools? Do you want to proceed to the next step (Y)? Or alternatively, you can manually carry out another BlobTools iteration ensuring that you name the outputs the same as the previous automated run (N): ' blobCheck
    case "$blobCheck" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo "You have chosen that your BlobTools results are satisfactory. Proceeding to running BUSCO..."
		break
		;;
      [Nn][Oo]|[Nn])  # No or N.
        #return 1
		echo "You have chosen to carry out further Blobtools iterations."
		echo "Once you are ready, answer Yes to the next question to proceed"
		read -p "Ready to proceed?: Y " blobRead
		if [[ $blobRead == "Yes" ]] || [[ $blobRead = "y" ]] || [[ $blobRead = "YES" ]] || [[ $blobRead = "yes" ]] || [[ $blobRead = "Y" ]]; then
			echo "Proceeding..."
			break
		else 
			echo "Your answer has caused a crash.. Please try again"
			exit
		fi
        ;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
  done
fnAssem=${WPATH}/SuccessfulAssemblies
mkdir -p $fnAssem
while read -r line;
do
    fname=$(basename "$line")    
    echo $fname
    echo
	#get the isolate name
	iso=${fname%.*}
	echo $iso
	mv $Assems/${iso} $fnAssem
	done < $jells/AssembliesToProcess.txt
echo "Your chosen assemblies have been moved to ${fnAssem}"
echo
echo "Decontamination completed. Moving to carry out assembly quality analysis..."
#................................................................
#					QUALITY ASSESSMENT OF ASSEMBLIES
#................................................................
echo
echo "First step is BUSCO..."
buscOut=${WPATH}/BUSCO_Output
quastOut=${WPATH}/Quast_Output
pomoOut=${WPATH}/Pomoxis
mkdir -p $LOGFF/BUSCO
mkdir -p $LOGFF/QUAST
mkdir -p $LOGFF/Pomoxis
mkdir -p $pomoOut
# Loop forever until the user enters a valid response (Y/N or Yes/No).
echo "To provide a point of comparison, your reference genome will also need to be taken through BUSCO..."
echo "If you have already carried this out or you do not wish to carry this out, you can skip this by answering No"
while true; do
	read -r -p ' Carry out BUSCO on your reference genome? Proceed? Y/N: ' refBUSstarting
    case "$refBUSstarting" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo "Beginning busco for your reference genome..."
		source ${SCPTS}/BUSCO.sh $PKREF $WPATH $buscOut $busDock $LINEAGE $THREADS
		refBUS 2>&1 | tee $LOGFF/BUSCO/referenceGenome.txt
		break
		;;
      [Nn][Oo]|[Nn])  # No or N.
        #return 1
		echo "Busco will not be run on your reference genome. Skipping to BUSCO your assemblies"
		break
        ;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
  done
echo
echo "Carrying out BUSCO for the assemblies you have chosen and cleaned using blobtools"
# call the busco script with the necessary variables
. ${SCPTS}/BUSCO.sh $PKREF $WPATH $buscOut $busDock $LINEAGE $THREADS
# call the functions needed
	cleanIso 2>&1 | tee $LOGFF/BUSCO/postBlobtools.txt
echo "BUSCO is done. The results of this are saved in ${buscOut}"
echo
echo "QUAST will now be carried out"
. ${SCPTS}/QUAST.sh $fnAssem $quastOut $PKREF $PKGFF $THREADS
	cleanQuas 2>&1 | tee $LOGFF/QUAST/postBlobtools.txt
echo "Quast done and the results are saved in ${quastOut}"
echo
echo "Some assembly analysis will be carried out on the cleaned assembly file"
. ${SCPTS}/pomoxis_analyse.sh $fnAssem $pomoOut $PKREF $THREADS
# call the raw assembly assessment functions
	rawPomo 2>&1 | tee $LOGFF/Pomoxis/postBlobtoolsPomoxis.txt
echo "Pomoxis analysis of the cleaned assembly is complete. Results are in ${pomoOut}/postAssembly"

#................................................................
#					POLISHING WITH RACON
#................................................................
#output
racOut=$WPATH/Racon_Polishing
racIter=${racOut}/IterationVsReads
mkdir -p $racIter
mkdir -p $STTS/RaconStats
mkdir -p $LOGFF/Racon
echo "Racon polishing starting on your chosen assemblies..."
${SCPTS}/Racon.sh $fnAssem $ALIFORM $racOut ${blobOut}/Clean/Reads_Vs_Assembly $racIter $THREADS $STTS/RaconStats 2>&1 | tee $LOGFF/Racon/AllRacon.txt
echo
echo "Racon polishing is complete"
echo "Racon outputs are saved in ${racOut}"
echo
echo "Some assembly analysis will be carried out on the cleaned assembly file"
. ${SCPTS}/pomoxis_analyse.sh $racOut $pomoOut $PKREF $THREADS
# call the racon polished assembly assessment functions
	raconPomo 2>&1 | tee $LOGFF/Pomoxis/postRaconPomoxis.txt
echo "Pomoxis analysis of the cleaned assembly is complete. Results are in ${pomoOut}/postRacon"
echo
echo "Next step is to carry out correction with medaka"
#................................................................
#					CORRECTION WITH MEDAKA
#................................................................
medOut=$WPATH/Medaka_Output
mkdir -p $medOut
mkdir -p $LOGFF/Medaka
medRun="${SCPTS}/Medaka.sh $ALIFORM $racOut $medOut ${STTS}/Medaka $THREADS 2>&1 | tee $LOGFF/Medaka/AllMedaka.txt"
echo $medRun
eval $medRun
echo
echo "Medaka complete"
echo
# remove redundant files:
check2r=$medOut/IterationVsReads
	if [[ -d "$check2r" ]]; then
		rm -r $check2r
	fi
check2y=$buscOut/IterationVsReads
	if [[ -d "$check2y" ]]; then
		rm -r $check2y
	fi
echo "Quality Assessment will be carried out on the now corrected assemblies"
#call the busco script with the necessary variables
. ${SCPTS}/BUSCO.sh $PKREF $WPATH $buscOut $busDock $LINEAGE $THREADS
# call the functions needed
	medyBus 2>&1 | tee $LOGFF/BUSCO/postMedaka.txt
echo "BUSCO is done. The results of this are saved in ${buscOut/Medaka}"
echo
echo "QUAST will now be carried out"
. ${SCPTS}/QUAST.sh $medOut $quastOut $PKREF $PKGFF $THREADS
	meddyQuas 2>&1 | tee $LOGFF/QUAST/postMedaka.txt
echo "Quast done and the results are saved in ${quastOut}"
echo
echo "Some assembly analysis will be carried out on the cleaned assembly file"
# . ${SCPTS}/pomoxis_analyse.sh $medOut $pomoOut $PKREF $THREADS
# call the medaka assembly assessment functions
	medakaPomo 2>&1 | tee $LOGFF/Pomoxis/postMedakaPomoxis.txt
echo "Pomoxis analysis of the cleaned assembly is complete. Results are in ${pomoOut}/postMedaka"
echo
echo "QA for Medaka is now complete. It will be advised for you to take this time to see your results"
echo
echo "If you have access to Illumina reads for the same isolates, you should have provided the paths to these in the script"
echo
echo "If you have provided these, you can proceed with the Pilon step after this"
echo
echo "If you do not have access to this, then this is as far as this script will take you"
echo
#................................................................
#					CORRECTION WITH PILON
#................................................................
while true; do
	read -r -p ' Proceed with Pilon correction step? Y/N: ' pilonStart
    case "$pilonStart" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo "Beginning Pilon for the isolates in the indicated path..."
		pilOut=$WPATH/Pilon_Output
		illumy=$WPATH/IllumVsAssembly
		mkdir -p $pilOut
		mkdir -p $illumy
		mkdir -p $LOGFF/Pilon
		pil="${SCPTS}/Pilon.sh $ILLUMPATH $medOut $pilOut $illumy $THREADS 2>&1 | tee $LOGFF/Pilon/AllPilon.txt"
		echo $pil
		eval $pil
		echo "Pilon for your isolates complete"
		echo "Carrying out Quast and BUSCO on these generated isolates"
		# BUSCO on Pilon output
		#call the busco script with the necessary variables
		. ${SCPTS}/BUSCO.sh $PKREF $WPATH $buscOut $busDock $LINEAGE $THREADS
		# call the functions needed
		pilBus 2>&1 | tee $LOGFF/BUSCO/postPilon.txt
		echo "BUSCO is done. The results of this are saved in ${buscOut/Pilon}"
		echo
		echo "QUAST will now be carried out"
		. ${SCPTS}/QUAST.sh $pilOut $quastOut $PKREF $PKGFF $THREADS
			pilQuas 2>&1 | tee $LOGFF/QUAST/postPilon.txt
		echo "Quast done and the results are saved in ${quastOut}"
		echo
		echo
		echo "Some assembly analysis will be carried out on the cleaned assembly file"
		. ${SCPTS}/pomoxis_analyse.sh $pilOut $pomoOut $PKREF $THREADS
		# call the pilon assembly assessment functions
		pilonPomo 2>&1 | tee $LOGFF/Pomoxis/postPilonPomoxis.txt
		echo "Pomoxis analysis of the cleaned assembly is complete. Results are in ${pomoOut}/postPilon"
		echo
		break
		;;
      [Nn][Oo]|[Nn])  # No or N.
        #return 1
		echo "Pilon will not run. This is the end of the script."
		break
        ;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
  done
#................................................................
#					ADVICE ON NEXT STEPS
#................................................................
echo
echo "If you have chosen to do Pilon, it is now complete" 
echo "Next step involves you ending this script here and carrying out some analyses outside the script"
echo "These analyses are personal to you. It is advised that you carry out the following steps to prepare your data"
echo "1. Remove apicoplast and mitochondrial genome for circularisation. This can be done manually."
echo "a. Use BLAST or align reference sequences to your assembly and remove what maps. Then use Circlator for circularisation followed by Prokka"
echo "2. Carry out repeatmasking with the repeatmasking script and carry out quality assessment check with indicated script"
echo "a. Use the RepeatMasking.sh script for this step. The script is set up to run automagically once set up although there might be some issues"
echo "b. You can use the MaskedQualityAssessment.sh script to check the quality of the masked genomes after repeatmasking"
echo "3. Check for chimeric contigs and scaffold contigs with RagTag"
echo "a. This is completely optional although there is a chance that you get chimeric contigs. Please check this with a tool of your choice"
echo "4. Upload onto Companion for genome annotation and prediction" 
echo "5. Alternatively carry out manual annotation using your own pipeline"
echo
echo
echo "This pipeline is now complete at this point. There are further scripts that can aid in carrying out further downstream processing. Please have a look in the README to carry on"
exit
