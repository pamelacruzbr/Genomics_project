""" This script allows parsing through a FASTA file with multi-entries using BioPython.
        It will retrieve information from Gene IDs, sequence and sequence lengths and store in lists.
        It will return the number of entries in the FASTA file """

# Import Biopython library
from typing import List

from Bio import SeqIO

ids = []
sequences = []
lengths = []

# Parse FASTA file

for seq_record in SeqIO.parse("dna.example.fasta", "fasta"):
        ids.append(seq_record.id) # List all IDs
        sequences.append(seq_record.seq) # List all sequences
        lengths.append(len(seq_record)) # List all lengths

#count number of entries in the FASTA

entries = print("Total entries in this FASTA file: {}".format(len(ids)))

# Find the biggest and smallest sequences

max_seq = max(lengths)
min_seq = min(lengths)

# Are there more than one sequences in the max and min ? Give their ID.

def count_seqmaxmin():
    global lengths
    print("The maximum length of sequence found in this FASTA is : %d bp" %max(lengths)) #Print MAX seq
    print("The minimum length of sequence found in this FASTA is : %d bp" %min(lengths)) #Print MIN seq
    if lengths.count(max(lengths))> 1: #Count if there are more than 1 seq with the MAX value
        print("There are %d with %d basepairs, the ID are: {}".format((ids[lengths.index(max_seq)])) %(lengths.count(max_seq), max_seq))
    else:
        print("There is only one sequence with %d basepairs, its ID is: {}".format((ids[lengths.index(max_seq)])) %max_seq)
    if lengths.count(min(lengths))> 1: #Count if there are more than 1 seq with the MIN value
        return print("There are %d with %d basepairs, the ID are: {}".format((ids[lengths.index(min_seq)])) %(lengths.count(min_seq), max_seq))
    else:
        print("There is only one sequence with %d basepairs, its ID is: {}".format((ids[lengths.index(min_seq)])) %min_seq)

count_seqmaxmin()

# Functionalities 1 e 2 of project: Done
