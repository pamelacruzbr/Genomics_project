
#This function counts the GC% in a DNA sequence

def gc(dna):
    nbases=dna.count('n')+dna.count('N') #Here, N or n characters are taken into account to be removed after
    gcpercent=float(dna.count('c')+dna.count('C')+dna.count('g')+dna.count('G'))//(len(dna)-nbases)
    return gcpercent
