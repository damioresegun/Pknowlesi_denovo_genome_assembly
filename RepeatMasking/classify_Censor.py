#!/usr/bin/env python2
# script to carry second step of the whole repeatmasking process. 
# This should take the output from CENSOR (Censor_results) and classify the repeats reported
# The output taken forward from this step is the /path/to/isolate/RM_*/consensi.fa.censor
# Written and adapted by Peter Thorpe and Damilola Oresegun
# Original commands come from TE discovery cookbook written by Amir Szitzenberg
# Original post: https://github.com/HullUni-bioinformatics/TE-search-tools
import sys
from TE_cookbook import *

print("""This is the 'classify_Censor' step

This should classify the types of repeats/TEs
that have been reported by CENSOR.
""")
# isolate = sys.argv[1] # the isolate name
# input = sys.argv[2] # fasta file
# print(input)
#########################################################
# classifying censore results

censor_classifications = parse_online_censor('Censor_results')

print_online_censor(censor_classifications,
                   'censor_classifications_file')

put_censor_classification_in_repeatmodeler_lib('consensi.fa.classified',
                                                censor_classifications,
                                                'consensi.fa.censor')
# if carrying out repeatmasking of multiple isolates of same species
print ("""
 If you do this with multiple isolates of the same species,
 You will need to combine all the censor files to make a single
 classified file and then run CD-HIT at 100% to reduce redundancy
 among the elements reported for each isolate in their
 consensi.fa.censor files.""")
exit()