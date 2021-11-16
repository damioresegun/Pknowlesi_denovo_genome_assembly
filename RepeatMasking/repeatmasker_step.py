#!/usr/bin/env python2
# script to carry fifth step of the whole repeatmasking process. 
# This should take the input assembly, combined consensi.classified.censor 
# file and the output directory, through the repeatmasker tool.
# Should then generate a directory containing sub-directories for each isolate.
# Sub-directories contain masked outputs which are not to be used at this point 
# Written and adapted by Peter Thorpe and Damilola Oresegun
# Original commands come from TE discovery cookbook written by Amir Szitzenberg
# Original post: https://github.com/HullUni-bioinformatics/TE-search-tools
import sys
from TE_cookbook import *

print("""This is the repeatmasker step
Repeatmasker will carry out first masking of the genome
based on the classified loci determined by the repeatmodeling
model and the censor outputs
""")
input = sys.argv[1] # fasta file
CeNsor = sys.argv[2]    # path to combined censor file
Outputy = sys.argv[3]   # path to output repeatmasker
print(input)
#########################################################

print("...beginning repeatmasker")
run_repeat_masker(input,
                  lib = CeNsor,
                  dir = Outputy)  
print ("repeatmasker done")
exit()



