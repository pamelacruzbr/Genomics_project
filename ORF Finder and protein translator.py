""" Finding ORFS in our FASTA sequences (Open Reading Frames) and translate them to proteins """

# First step, parse in the sequences
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

# Start ORF function

from Bio.Data import CodonTable
table = CodonTable.unambiguous_dna_by_name["Standard"]
min_pro_len = 100
proteins = []

def parse_sequences(choice):
    global sequences
    if choice > 2 or choice < 0:
        print("Error: A reading frame has to be between 1 and 3")
    else:
        for sequence in sequences:
            for strand, base in [(1, sequence), (-1, sequence.reverse_complement())]:
                for frame in range(3):
                    length = 3 * ((len(sequence) - frame) // 3) # Multiple of three
                    for protein in base[frame : frame + length].translate(table).split("*"):
                        if len(protein) >=50 and protein[0] == 'M':
                            prot_base = sequences.index(sequence)
                            print(
                                "%s..%s - ID {}, length %i, strand %i, frame %i".format(ids[prot_base])
                                %(protein[:3], protein[-3:], len(protein), strand, frame)
                            )
parse_sequences(1)

# Identify ALL ORFS (ATG.. TAA,TAG,TGA
# Length of ORF, longest and its ID
# A given ID, find longest ORF, find its start position (frame)

""" Regex solution"""

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

sequencestr = str(sequences)

pattern = re.compile(r'(ATG(?:...)*?(?:TAG|TGA|TAA))')

ORF_total = []

def ORF_finder(frame=0):
    if frame > 2 or frame < 0:
        print("Error: frame has to be between 0 and 2")
    else:
        for seq in sequences:
            strand_positve = pattern.findall(str(seq[frame:]))
            revcomp = seq.reverse_complement()
            strand_negative = pattern.findall(str(revcomp[frame:]))
            ORF_total.append(strand_positve + strand_negative)
    return ORF_total

ORF_finder(2)
print(len(ORF_total))

for ORF in ORF_total:
                print
                ident = print(
                ids[sequences.index(seq)])
                for ORF in ORF_total:
                    print('\nTotal ORFS = ',cnt,'\n', ORF,'\nlenght = ', len(ORF),'\n---------------------------\n\n')

def ORF_informations():
    ORF_list = ORF_finder(2)
    max_ORF = max(ORF_list)
    print(len(max_ORF))

ORF_informations()
