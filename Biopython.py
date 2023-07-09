""" This script allows BLAST alignments using Biopython.
    It uses a Fasta file containing a DNA sequence of interest for input.
    All BLAST parameters are standard """


import Bio
from Bio.Blast import NCBIWWW

fasta_string = open("../Studies/Randomseq.fa", 'r').read() #read sequence

#blast command form NCBI
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)

#fomartting BLAST results XML
from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle)

#Parsing BLAST output
len(blast_record.alignments)

E_VALUE_TRESH = 0.01
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_TRESH:
            print('****Alignment****')
            print('sequence:', alignment.title)
            print('length:', alignment.title)
            print('e value:', hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)

from Bio.Seq import Seq

print('reverse complement is %s' % my_seq.reverse_complement())

