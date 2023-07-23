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
    return {
        'seq_fasta': seq_fasta,
        'ids': ids,
        'lengths': lengths,
        'sequences': sequences
    }

def count_seqmaxmin():
    lengths = import_fasta()['lengths']
    ids = import_fasta()['ids']
    max_seq = max(lengths)
    min_seq = min(lengths)
    print("Your FASTA file has {} entries \n ----------------- \n".format(len(ids)))
    print("The maximum length of sequence found in this FASTA is : %d bp \n ----------------- \n" %max(lengths)) #Print MAX seq
    print("The minimum length of sequence found in this FASTA is : %d bp" %min(lengths)) #Print MIN seq
    if lengths.count(max(lengths))> 1: #Count if there are more than 1 seq with the MAX value
        print("There are %d with %d basepairs, the ID are: {}".format((ids[lengths.index(max_seq)])) %(lengths.count(max_seq), max_seq))
    else:
        print("There is only one sequence with %d basepairs, its ID is: {}".format((ids[lengths.index(max_seq)])) %max_seq)
    if lengths.count(min(lengths))> 1: #Count if there are more than 1 seq with the MIN value
        return print("There are %d with %d basepairs, the ID are: {}".format((ids[lengths.index(min_seq)])) %(lengths.count(min_seq), max_seq))
    else:
        print("There is only one sequence with %d basepairs, its ID is: {} \n ----------------- \n".format((ids[lengths.index(min_seq)])) %min_seq)


def ORF_finder2(data, frame=3):
    seq_fasta = data
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
                    if length >= lengthorf and length < 2000:
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
                    if length >= lengthorf and length < 2000:
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
        data = import_fasta(file)['seq_fasta']
    elif typeseq == 'Single seq':
        data = {}
        identif = ''
        identif = input('Give the name of the gene')
        seq_fasta = import_fasta()['seq_fasta']
        data[identif] = seq_fasta.get(identif)
    else: print('error')
    return data


def ORF_informations():
    data = user_choice()
    sequence = list(data.values())[0]
    if len(data) == 1:
        print(" Your sequence has a length of: {}".format(len(sequence)))
    else:
        count_seqmaxmin()
    print("Longest ORF in reading frame 1: {}, its start position is: {}, id of the gene: {} \n ".format(ORF_finder2(data, frame=1)['length'], ORF_finder2(data, frame=1)['start'], ORF_finder2(data, frame=1)['name']))
    print("Longest ORF in reading frame 2: {}, its start position is: {}, id of the gene: {} \n ".format(ORF_finder2(data, frame=2)['length'], ORF_finder2(data, frame=2)['start'], ORF_finder2(data, frame=2)['name']))
    print("Longest ORF in reading frame 3: {}, its start position is: {}, id of the gene: {} \n ".format(ORF_finder2(data, frame=3)['length'], ORF_finder2(data, frame=3)['start'], ORF_finder2(data, frame=3)['name']))


ORF_informations()
