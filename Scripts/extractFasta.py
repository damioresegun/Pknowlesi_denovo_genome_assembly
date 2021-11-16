# Script to extract identified reads or contigs from a fasta file based on the read name/ID in a txt file
# The script is designed to take in a list of read names and then search for these in the fasta file, extract them and save in a new fasta file
# Note that this does not delete the sequences from the original fasta file

# Usage: python extractFasta.py readsListTxtFile file.fasta output.fasta

from Bio import SeqIO
import sys

readsFile = open(sys.argv[1], 'r')
#readsFile = open(sys.argv[1], 'r')
fastaFile = sys.argv[2]
outputFile = open(sys.argv[3], 'w')

wanted = set()
with readsFile as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)

fasta_sequences = SeqIO.parse(open(fastaFile),'fasta')

with outputFile as i:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], i, "fasta")

readsFile.close()
outputFile.close()

# from Bio import SeqIO
# import sys

# readsList = open(sys.argv[1], 'rU')
# fastafile = sys.argv[2]
# outputfile = open(sys.argv[3], 'w')

# wanted = set()
# with readsList as f:
    # for line in f:
        # line = line.strip()
        # if line != "":
            # wanted.add(line)

# fasta_sequences = SeqIO.parse(open(fastafile),'fasta')

# with outputfile as i:
    # for seq in fasta_sequences:
        # if seq.id in wanted:
            # SeqIO.write([seq], i, "fasta")

# readsList.close()
# outputfile.close()
