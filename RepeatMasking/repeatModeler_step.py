#!/usr/bin/env python2
# script to carry first step of the whole repeatmasking process. 
# This should take the input assemblies through the repeatmodelling and repeatscout software
# Should then generate a directory of TE models for each assembly. 
# The output taken forward from this step is the /path/to/isolate/RM_*/consensi.fa.classified
# Note: it is important that the assembly contigs have been renamed prior to starting this script
# Written and adapted by Peter Thorpe and Damilola Oresegun
# Original commands come from TE discovery cookbook written by Amir Szitzenberg
# Original post: https://github.com/HullUni-bioinformatics/TE-search-tools

import sys
from TE_cookbook import *

print("""This is the 'Run RepeatModeler step'

RepeatModeler will produce consensus sequences representing
clusters of de novo repeat sequences
""")
isolate = sys.argv[1] # the isolate name
input = sys.argv[2] # fasta assembly file
print(input)
print ("Have you renamed your contigs/scaffolds using the 'rename_fasta.py'?")
print ("If not, please stop this script NOW and do that")
#########################################################
# run
make_repeatmodeler_database(name=isolate+'_lib.DB',
                            input_filename=input)
print ("make_repeatmodeler_database done")

print ("...run_repeatmodeler")
run_repeatmodeler(isolate+'_lib.DB') 

print ("""RepeatModeler complete. Using the firefox brower,
upload the result (RM_*/consensi.fa.classified) from within the
repeatmodeler directory onto CENSOR online.
Select the closest kingdom/family option, 
select 'report simple repeats'. 
Copy the resulting webpage into a txt file named 'Censor_results'.
Note: If the txt file is created on Windows, it is possible 
the next step might break, so it is advised to open the txt file
on a Linux based system.
CENSOR is found at this link:http:/www.girinst.org/censor/
 """)
exit()
#########################################################

