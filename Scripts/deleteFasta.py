#!/usr/bin/env python3
# Script to delete identified reads or contigs from a fasta file based on the read name/ID in a txt file
# The script is designed to take in a list of read names and then search for these in the fasta file, delete them then save remainder of the fasta file in a given path
#usage : deleteFasta.py input.fasta list_ofReads.txt > filtered.fasta
from Bio import SeqIO
import sys

ffile = SeqIO.parse(sys.argv[1], "fasta")
header_set = set(line.strip() for line in open(sys.argv[2]))

for seq_record in ffile:
    try:
        header_set.remove(seq_record.name)
    except KeyError:
        print(seq_record.format("fasta"))
        continue
if len(header_set) != 0:
    print(len(header_set),'of the headers from list were not identified in the input fasta file ',header_set, file=sys.stderr)
