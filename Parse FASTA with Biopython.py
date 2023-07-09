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

entries = print(len(ids))
