#!/usr/bin/env python2
# script to carry ninth step of the whole repeatmasking process. 
# This should take the isolate ID, input assembly, repeatmasker outputs,
# ltrharvest outputs and transposonPSI outputs 
# Should then generate a directory containing sub-directories for each isolate.
# Sub-directories contain GFF files of the loci to mask 
# Written and adapted by Peter Thorpe and Damilola Oresegun
# Original commands come from TE discovery cookbook written by Amir Szitzenberg
# Original post: https://github.com/HullUni-bioinformatics/TE-search-tools
import sys
from TE_cookbook import *

"""This is the step to generate the non-redundant repeat loci.
GFF files are saved which can then be used to mask the assembly
downstream.
"""
isolate = sys.argv[1] # the isolate name
input = sys.argv[2] # fasta file
Outputy = sys.argv[3]   # path to output repeatmasker
ltrHavy = sys.argv[4]   # path to ltrharvest file for the isolate
TPSIr = sys.argv[5]     # path to the transposonPSI file of the isolate

print(input)
#########################################################
# puting OneCodeToFindThemAll parsed RepeatMasker results in the data structure
print ("parsing OCTFTA csv files")
TEs_from_OCFA, serial = parse_ocfa_elem_stored(Outputy)
#
print ("DONE parsing OCTFTA csv files")
# adding LTRharvest results to the data structure
print ("reduce LTRharvest into the data structure")
TEs, serial = integrate_ltrharvest_to_RM_TEs(ltrHavy,      
                                              input,  
                                              TEs_from_OCFA,
                                             serial)
print ("reduce transposonPSI into the data structure")
TEs = integrate_TransposonPSI_to_RM_TEs(TPSIr, 
                                         input, 
                                         TEs, 
                                         serial)
print ("write non-redundant loci to GFF")
write_gff(TEs, isolate+'.transposable_elements.gff3', max_RM_OC_score=False)
write_gff(TEs, isolate+'.transposable_elements_RM_Score.gff3', max_RM_OC_score=True)
exit()

############################################################################################################