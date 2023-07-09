#This is a program that reads fastas and organizes sequences in a dictionary

try:
    filename = 'dna.example.fasta'
    f = open(filename, 'r') #open and read fasta file
except IOError:
    print("File %s does not exist!!" % filename) #stderr

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

f.close


print((seqs))



