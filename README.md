<!-- Auto generated table of contents (TOC) -->

[TOC]

<!-- tocstop -->
# P.knowlesi de novo genome assembly
## Abstract:
A pipeline to take input data from raw nanopore reads to being ready for manual or automated annotation. 
It is important to note that while this pipeline automates alot of commands, it cannot automate everything. 
User input is needed at various points to confirm which branch of the pipeline to follow as well as when to continue at certain 'break-points'
I'm hoping that this has been fairly prepped in a fairly straight-forward manner that is somewhat approachable
However, if you want to tweak the script, the vast majority is written in Bash with some secondary scripts in python.
Please be aware that the first portion of this script is tuned for barcoded ONT sequencing reads however, it is possible to use non-barcoded ONT sequenced read. 

