# First step, parse in the sequences
# Import Biopython library
from typing import List

from Bio import SeqIO
import re

ids = []
sequences = []
lengths = []

# Parse FASTA file

for seq_record in SeqIO.parse("dna.example.fasta", "fasta"):
        ids.append(seq_record.id) # List all IDs
        sequences.append(seq_record.seq) # List all sequences
        lengths.append(len(seq_record)) # List all lengths


seq_fasta = dict(map(lambda gene,seq : (gene,seq) , ids, sequences))



pattern = re.compile(r'((?:...)*?)(ATG(?:...)*?(?:TAG|TGA|TAA))')


def ORF_finder2(frame=3):
    lengthorf = 1
    orf = {
        'length': 1,
        'start': 0,
        'end': 0,
        'name': '',
        'strand': ''
    }
    if frame > 3 or frame <= 0:
        print("Error: frame has to be between 0 and 2")
    else:
        for name, seq in seq_fasta.items():
            for match in re.finditer(pattern, str(seq[frame-1:])):  #search ORF in positive strand
                group1 = match.group(1)
                group2 = match.group(2)
                if group1 is not None:
                    start = match.start() + len(group1)
                    end = 0 + match.end()
                    length = len(group2)
                    if length >= lengthorf:
                        lengthorf = length
                        orf = {
                            'length': length,
                            'start': start + frame,
                            'end': end,
                            'name': name,
                            'seq': group2,
                            'strand': 'positive'
                        }
            for match in re.finditer(pattern, str(seq.reverse_complement()[frame-1:])):  #search ORF in negative strand
                group1 = match.group(1)
                group2 = match.group(2)
                if group1 is not None:
                    start = match.start() + len(group1)
                    end = 0 + match.end()
                    length = len(group2)
                    if length >= lengthorf:
                        lengthorf = length
                        orf = {
                            'length': length,
                            'start': start + frame,
                            'end': end,
                            'name': name,
                            'seq': group2,
                            'strand': 'negative'
                        }
        return orf
