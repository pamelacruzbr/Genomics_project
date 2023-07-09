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
