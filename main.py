#This is a program that reads fastas and organizes sequences in a dictionary

try:
    f = open(r'C:\Users\pamel\Desktop\Data_science\Data\protfasta.fasta', 'r')
except IOError:
    print("File protfasta does not exist!!")

seqs = {}
for line in f:
    #discard the newline at the end (if any)
    line=line.rstrip()
    #distinguish header from sequence
    if line[0] =='>': #line.startwith('>')
        words=line.split()
        name=words[0][1:]
        seqs[name]=''
    else: #sequence, not header
        seqs[name] = seqs[name] + line
