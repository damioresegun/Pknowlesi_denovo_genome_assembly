#!/bin/bash
# Script to carry out repeatmasking of identified assemblies
# Multiple steps are being taken with different python scripts
# You will need to have installed OneCodeToFindThemAll locally
# Note: you can easily comment out each steps if you want to skip them for steps later on in the script
### set paths and variables
# set path to assemblies to mask
aPATH="/path/to/assembles/you/want/to/mask"
# set path to save directory
sPATH="/path/to/output/directory"
# set path to scripts location
scpts="path/to/directory/for/scripts/RepeatMasking"
# set path to one code to find them all directory holding the scripts
ocfta="path/to/where/you/installed/Tools/Onecodetofindthemall"
# ONLY fill if transposonPSI in your conda environment fails. Enter path to transposonPSI perl script
trapsi=""	# path/to/transposonPSI.pl
# set number of threads
## note that the default thread number is 32
THREADS=32
############################################################################
# step 1 - rename contigs incrementally
#### there are downstream steps which may break if the contig ID are too long
while true; do
	read -r -p 'Step 1 of repeatmasking protocol: The first step reformats your assemblies for use. It does not affect your original files. Proceed?: ' step1
    case "$step1" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		LOGFF=$sPATH/Logs
		mkdir -p $LOGFF
		echo "Beginning repeatmasking process"
		echo 
		echo "First step is to rename the contigs to incremental numbers"
		for i in $aPATH/*; 
			do cl=$(basename "$i"); 
			clr=${cl%.*};
			outr=$sPATH/ReWrittenAssemblies;
			mkdir -p $outr;
			python $scpts/rename_fasta.py -i $i	-o $outr/${clr}.fasta --pre contig_ 2>&1 | tee $LOGFF/Renaming_Contigs.txt; 
		done
		echo "Rewritten assemblies are now saved in ${sPATH}/ReWrittenAssemblies."
		echo
		echo "Renaming complete. Moving on to repeatmodeler step"
		break
		;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
done
###############################################################################
# Step 2 - repeatmodeling
while true; do
	read -r -p 'Step 2 of repeatmasking protocol: This carries out repeat element prediction with repeatmodeler Proceed?: ' step2
    case "$step2" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo "Begining repeatmodeler step..."
		echo
		for i in $aPATH/*; 
			do foi=$(basename "$i");
			foiu=${foi%.*}; 
			echo $foiu; 
			yup=$sPATH/RepeatModeler/${foiu}; 
			mkdir -p $yup; 
			cd $yup; 
			$scpts/repeatModeler_step.py $foiu $i 2>&1 | tee $LOGFF/RepeatModeller_Contigs.txt; 
			cd; 
		done
		echo
		echo "RepeatModeler complete"
		echo
		echo "RepeatModeler outputs are saved in ${sPATH}/RepeatModeler"
		echo
		break
		;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
done
###############################################################################
#step 3 - classify censor results
# ask if they have uploaded to censor and downloaded results
# it must be done on FIREFOX due to line wrapping issues with Chrome/Chromium brwosers
# must name the results "Censor_results" and palce in RM* folder in the repeatmodeler directory
# Then run censor classification step with:
echo "To proceed to next step, you will need to complete a step OUTSIDE of this script."
echo "Please leave the script running..."
echo "Before the next step of the script, you have to upload the results of repeatmodeler to CENSOR"
echo "To do this, find the consensi.fa.classfied e.g ${sPATH}/RepeatModeler/isolateID/RM_*/consensi.fa.classfied"
echo "Using the FIRFOX browser, Upload the consensi.fa.classfied to CENSOR at 'https://www.girinst.org/censor/'"
echo "Select the most appropriate sequence source e.g. Eukaryota for Plasmodium"
echo "Select 'Report simple repeats' and submit"
echo "Copy the resulting webpage into a blank txt file saved as ${sPATH}/RepeatModeler/isolateID/RM_*/Censor_results for each isolate"
echo "It is important that CENSOR is run in a FIREFOX browser"
echo "When complete, agree to the question below to proceed"
while true; do
	read -r -p 'Step 3 of repeatmasking protocol: You need to have carried out the Censor analysis steps above BEFORE proceeding. Have you done that?: ' censoree
    case "$censoree" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo "Beginning the censor classification..."
		for i in $aPATH/*; 
			do foi=$(basename "$i"); 
			foiu=${foi%.*}; 
			echo $foiu; 
			cd $sPATH/RepeatModeler/${foiu}/RM*; 
			dos2unix Censor_results; 
			$scpts/classify_Censor.py $foiu $i 2>&1 | tee $LOGFF/censor_classification.txt; 
			cd;
		done
		break
		;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
  done
echo
echo "Censor classification complete..."
echo
echo "Censor classifications are saved as ${sPATH}/RepeatModeler/isolateID/RM_*/consensi.fa.censor"
echo "proceeding to CD-HIT"
#################################################################################################################
# step 4 - CD-HIT
# next step is to do CD-HIT analysis
# if multiple isolates OF the same species are being masked then the classfied censor outputs must be conc'd
# make CDHIT output directory
while true; do
	read -r -p 'Step 4 of repeatmasking protocol: Removes redundancies to make a single repeat library. Proceed? ' step4
    case "$step4" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo "To begin, the output of the censor classification will be concatenated"
		echo
		mkdir -p $sPATH/CDHIT
		# concatenate all isolate censor outputs into one file to run CDHIT to remove redundancy
		for i in $sPATH/RepeatModeler/*/RM*/*.fa.censor; 
		do 
		cat $i >> $sPATH/RepeatModeler/CombinedIsolates.censor.fa; 
		done
		echo
		echo "The combined consensi.fa.censor file is now saved as $sPATH/RepeatModeler/CombinedIsolates.censor.fa"
		echo
		echo "beginning CD-HIT..."
		 # run CDHIT
		cd-hit-est -i $sPATH/RepeatModeler/CombinedIsolates.censor.fa -o $sPATH/CDHIT/Pk_RepeatLib.censor.fa -c 1.0 -n 10 -d 0 -g 1 -M 2000 -T $THREADS 2>&1 | tee $LOGFF/CDHIT.txt
		echo
		echo "CDHIT complete"
		echo
		echo "CDHIT output is saved as ${sPATH}/CDHIT/Pk_RepeatLib.fa. Output is the repeat library for all isolates to be masked"
		echo
		echo "proceeding to repeatmasker"
		break
		;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
 done
####################################################################################################################
# step 5 - repeatmasker step
while true; do
	read -r -p 'Step 5 of repeatmasking protocol: The repeat library is used to call repeats using repeatmasker. Proceed? ' step5
    case "$step5" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo "begining repeatmasker..."
		for i in $aPATH/*; 
		do foi=$(basename "$i"); 
		foiu=${foi%.*}; 
		mkdir -p $sPATH/RepeatMasking/${foiu}; 
		$scpts/repeatmasker_step.py $i $sPATH/CDHIT/Pk_RepeatLib.censor.fa $sPATH/RepeatMasking/${foiu} 2>&1 | tee $LOGFF/Repeatmasking_step.txt; 
		done
		echo
		echo "repeatmasking complete... Results are saved in ${sPATH}/RepeatMasking"
		echo
		echo "proceeding to one code to find them all"
		break
		;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
 done
##############################################################################
#step 6 - one code to find them all
while true; do
	read -r -p 'Step 6 of repeatmasking protocol: OneCodeToFindThemAll does its thing. Proceed? ' step6
    case "$step6" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo "begining one code to find them all"
		for i in $sPATH/RepeatMasking/*; 
		do ls $i; 
		rfc=$(basename "$i"); 
		cd $i; 
		$ocfta/build_dictionary.pl --rm $i/${rfc}.fasta.out --fuzzy > ${rfc}_fuzzy.txt 2>&1 | tee $LOGFF/OCTFTA_dictionary.txt; 
		$ocfta/one_code_to_find_them_all.pl --rm $i/${rfc}.fasta.out --ltr ${rfc}_fuzzy.txt --fasta 2>&1 | tee $LOGFF/OCTFTA_library.txt; 
		done
		echo
		echo "ocfta complete. Output are in the isolate directories in ${sPATH}/RepeatMasking/"
		echo
		echo "proceeding to carry out ltrharvest"
		echo
		break
		;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
 done
########################################################################################
# step 7 - ltr harvest
while true; do
	read -r -p 'Step 7 of repeatmasking protocol: LTRHarvest will begin. Proceed? ' step7
    case "$step7" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo "beginning ltrharvest..."
		for i in $aPATH/*; 
		do puf=$(basename "$i"); 
		pufr=${puf%.*}; 
		echo $pufr; 
		outy=$sPATH/LTRHarvest/${pufr}; 
		mkdir -p $outy; 
		cd $outy; 
		gt suffixerator -db $i -indexname $pufr -tis -suf -lcp -des -ssp -sds -dna 2>&1 | tee $LOGFF/ltrharvest_suffix.txt; 
		gt ltrharvest -index $pufr -mintsd 5 -maxtsd 100 > ${pufr}_ltr.out 2>&1 | tee $LOGFF/ltrharvest.txt; 
		cd; 
		done
		echo
		echo "ltrharvest complete... Outputs are saved in ${sPATH}/LTRHarvest/"
		echo
		echo "proceeding to transposonPSI"
		break
		;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
 done
####################################################################################
# step 8 - transposonPSI
# ask for the specific path to transposonPSI 
while true; do
	read -r -p 'Step 8 of repeatmasking protocol: TransposonPSI searches. Proceed? ' step8
    case "$step8" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo "beginning transposonPSI.."
		echo "Please be aware, the script may or may not break here due to attempts to find transposonPSI paths"
		echo "If this does break and you see an error, then best way to fix it is to manually download transposonPSI"
		echo "If downloaded then please go to line 249 in the script for further instructions on how to proceed"
		echo "If the script works then, GREAT"
		# in case the script breaks, follow the instructions in the comments below.
		# download transposonPSI locally and enter path into variable at the top of this script
		# comment lines 252-255 (inclusive), go to line 263 and then line 265
		ftty=`which transposonPSI.pl`; 
		ftt=$(dirname "$ftty");
		fti=$(dirname "$ftt")
		ftyy=${fti}/share/transposonPSI/transposonPSI.pl
		for i in $aPATH/*; 
		do puf=$(basename "$i"); 
		pufr=${puf%.*}; 
		echo $pufr; 
		outy=$sPATH/TransposonPSI/${pufr}; 
		mkdir -p $outy; 
		cd $outy; 
		# comment this line if conda transposonPSI did not work and you are using transposonPSI locally
		$ftyy $i nuc 2>&1 | tee $LOGFF/transposonPSI.txt; 
		# uncomment this line if conda transposonPSI did not work and you are using transposonPSI locally
		#$trapsi $i nuc 2>&1 | tee $LOGFF/transposonPSI.txt;
		cd; 
		done
		echo
		echo "TransposonPSI complete... Results are saved in ${sPATH}/TransposonPSI/"
		echo
		echo "proceeding to parsing and reducing redundancy"
		break
		;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
 done
###################################################################################################################################################
#step 9 - removing redundancy
while true; do
	read -r -p 'Step 9 of repeatmasking protocol: All independent repeat elements are collapsed to be non-redundant. Proceed? ' step9
    case "$step9" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo
		echo "beginning to remove redundancy..."
		for i in $aPATH/*; 
		do foit=$(basename "$i"); 
		foitr=${foit%.*}; 
		echo $foitr; 
		outr=$sPATH/LociToMask/${foitr}; 
		mkdir -p $outr; 
		cd $outr; 
		$scpts/generate_mask_loci.py $foitr $i $sPATH/RepeatMasking/${foitr} $sPATH/LTRHarvest/${foitr}/${foitr}_ltr.out	$sPATH/TransposonPSI/${foitr}/${foitr}.fasta.TPSI.allHits.chains.bestPerLocus 2>&1 | tee $LOGFF/redundancyRemoval.txt; 
		cd; 
		done
		echo "redundancy removal complete. Results are in ${sPATH}/LociToMask"
		echo
		echo "proceeding to final masking"
		break
		;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
 done
###################################################################################################################################################
# step 10 - mask the assemblies based on the loci gffs
while true; do
	read -r -p 'Step 10 of repeatmasking protocol: Assemblies will be masked now. Proceed? ' step10
    case "$step10" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        #return 0
		echo
		echo "beginning masking of identified repeats..."
		for i in $aPATH/*; 
		do bmk=$(basename "$i"); 
		echo $bmk; 
		bmkr=${bmk%.*}; 
		echo $bmkr; 
		outr=$sPATH/MaskedAssemblies; 
		mkdir -p $outr; 
		bedtools maskfasta -fi $i -fo ${outr}/${bmkr}_masked.fasta -bed $sPATH/LociToMask/${bmkr}/${bmkr}.transposable_elements.gff3 -soft 2>&1 | tee $LOGFF/finalMasking.txt;
		done
		echo
		echo "Masking complete. Soft masked assemblies are saved in ${sPATH}/MaskedAssemblies"
		echo 
		echo "Masked assemblies are ready for protein prediction and annotation"
		echo
		echo "However, I would suggest checking for chimeric contigs using ragtag" 
		echo "Thank you"
		break
		;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
 done 