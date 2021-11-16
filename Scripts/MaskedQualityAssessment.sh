#!/bin/bash
set -e
#############################################################################################################################################################################
#													Automation script for carrying out quality assessment of masked, de novo assemblies										#	#																																											#
#																						Version 0.1.0																		#
# 								Designed to take in masked fasta files from a designated directory, carry out basic statistics, QUAST, BUSCO and pomoxis					#
# 											This script works as a primary script which calls sesourcery scripts to carry out the analyses specified. 						#
#											The user can indicate which processes they want to carry out in the configuration step of the script							#
#############################################################################################################################################################################
#.........................................................................................................
#													CONFIGURATION STEPS
#.........................................................................................................
# Path to pipeline scripts
SCPTS=/path/to/PipelineScripts
# Path to your directory holding your masked assemblies
maskedAssem=/path/to/MaskedAssemblies
# Path to directory to save analyses
outputdir=/path/to/outputs/MaskedAssemblies_QA
# Number of threads
THREADS="18" # no maximum thread count
# Where is your reference genome and the gff associated
PKREF=/path/to/your/reference.fasta # pknowlesi reference
PKGFF=/path/to/your/reference_annotations.gff # PKNOWLESI GFF
# Where is the busco docker image
busDock=/path/to/BUSCO/singularity/image/busco_v5.0.0_cv1.sif
# What is the lineage you wish to use for the busco assessment?
LINEAGE="plasmodium_odb10"

echo "This script must be run in a conda environment with quast and assembly-stats!"

while true; do
	read -r -p ' Continue with the script? Is your assembly masked? Y/N: ' maskStart
    case "$maskStart" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
		LOGFF=$outputdir/LogFiles
		mkdir -p $LOGFF
		# BUSCO on masked assembly
		mkdir -p $outputdir
		#call the busco script with the necessary variables
		. ${SCPTS}/BUSCO.sh $PKREF $maskedAssem $outputdir $busDock $LINEAGE $THREADS
		# call the functions needed
		maskBusScript 2>&1 | tee $LOGFF/BUSCO.txt
		echo "BUSCO is done. The results of this are saved in ${buscOut/MaskedAssemblies}"
		echo
		echo "QUAST will now be carried out"
		quastOut=$outputdir/QUAST
		mkdir -p $quastOut
		. ${SCPTS}/QUAST.sh $maskedAssem $quastOut $PKREF $PKGFF $THREADS
		maskQuasScript 2>&1 | tee $LOGFF/QUAST.txt
		echo "Quast done and the results are saved in ${quastOut}"
		echo
		echo
		echo "Some assembly analysis will be carried out on the masked assembly file"
		pomoOut=$outputdir/Pomoxis
		mkdir -p $pomoOut
		. ${SCPTS}/pomoxis_analyse.sh $maskedAssem $pomoOut $PKREF $THREADS
		maskedPomoScript 2>&1 | tee $LOGFF/Pomoxis.txt
		echo
		echo "Pomoxis complete. Assembly stats to be done"
		assembly-stats $maskedAssem/* >> $outputdir/All_MaskedAssemblies_assemblyStats.txt
		cat $outputdir/All_MaskedAssemblies_assemblyStats.txt | tr -s '[:blank:]' ',' > $outputdir/All_MaskedAssemblies_assemblyStats.csv
		break
		;;
      [Nn][Oo]|[Nn])  # No or N.
        #return 1
		echo "You are not ready to start. This is the end of the script."
		break
		exit
        ;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
  done
