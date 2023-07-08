import random
def create_dna(n, alphabet='actg'):
    return ''.join([random.choice(alphabet) for i in range(n)])
