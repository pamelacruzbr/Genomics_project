""" Regex solution 2"""

# First step, parse in the sequences
# Import Biopython library
from typing import List
from Bio import SeqIO
import re

# Parse FASTA file
def import_fasta(file="dna.example.fasta"):
    ids = []
    sequences = []
    lengths = []
    for seq_record in SeqIO.parse(file, "fasta"):
            ids.append(seq_record.id) # List all IDs
            sequences.append(seq_record.seq) # List all sequences
            lengths.append(len(seq_record)) # List all lengths
    seq_fasta = dict(map(lambda gene,seq : (gene,seq) , ids, sequences))
    return seq_fasta


def ORF_finder2(frame=3):
    seq_fasta = user_choice()['data']
    pattern = re.compile(r'((?:...)*?)(ATG(?:...)*?(?:TAG|TGA|TAA))')
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

#ROUTER !!
def user_choice():
    file = ''
    typeseq = input("Single seq or FASTA ? (seq or FASTA)")
    if typeseq == 'FASTA':
        file = input('Give the name of the FASTA file')
        data = import_fasta(file)
    elif typeseq == 'Single seq':
        data = {}
        identif = input('Give the name of the gene')
        seq_fasta = import_fasta()
        data[identif] = seq_fasta.get(identif, "identif is not a key")
    else: print('error')
    return data

print(user_choice())
