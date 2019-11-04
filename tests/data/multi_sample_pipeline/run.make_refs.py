#!/usr/bin/env python3

import random
import pyfastaq

random.seed(1)

nucleotides = ['A', 'C', 'G', 'T']
seq0 = [random.choice(nucleotides) for _ in range(1000)]
seq1 = seq0[:800] + ['A'] + seq0[800:]
seq2 = seq0[:600] + seq0[630:]
seq0[149] = 'G'
seq1[149] = 'A'
seq2[149] = 'T'
seq0[449] = 'T'
seq1[449] = 'C'
seq2[449] = 'T'
seq0[609] = 'A'
seq1[609] = 'G'
seq0[74] = 'A'
seq1[74] = 'A'
seq2[74] = 'G'

seqs = [seq0, seq1, seq2]
for i, seq in enumerate(seqs):
    name = 'ref.' + str(i)
    fasta = pyfastaq.sequences.Fasta(name, ''.join(seq))
    with open('run.' + name + '.fa', 'w') as f:
        print(fasta, file=f)

