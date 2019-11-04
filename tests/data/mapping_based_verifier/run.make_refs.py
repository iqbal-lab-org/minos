#!/usr/bin/env python3

import copy
import random
import pyfastaq

random.seed(42)

def mutated_nucleotide(nucleotide):
    other_nucleotides = ['A', 'C', 'G', 'T']
    other_nucleotides.remove(nucleotide)
    return random.choice(other_nucleotides)


nucleotides = ['A', 'C', 'G', 'T']
seq1 = [random.choice(nucleotides) for _ in range(1000)]
seq2 = [random.choice(nucleotides) for _ in range(1000)]
seq3 = seq1[:300] + [random.choice(nucleotides) for _ in range(700)]

ref1 = pyfastaq.sequences.Fasta('ref.1', ''.join(seq1))
ref2 = pyfastaq.sequences.Fasta('ref.2', ''.join(seq2))
ref3 = pyfastaq.sequences.Fasta('ref.3', ''.join(seq3))


mutated_seq1 = copy.copy(seq1)
mutated_seq1[99] = mutated_nucleotide(mutated_seq1[99])
mutated_seq1.pop(297)
mutated_seq1.pop(297)
mutated_seq1.insert(499, 'AGTC')
to_delete = mutated_seq1[699:701]
insertion = [mutated_nucleotide(x) for x in to_delete] + ['C']
mutated_seq1 = mutated_seq1[:699] + insertion + mutated_seq1[701:]

mutated_ref1 = pyfastaq.sequences.Fasta('mutated.ref.1', ''.join(mutated_seq1))

with open('run.ref.fa', 'w') as f:
    print(ref1, ref2, ref3, sep='\n', file=f)


ref2.id = 'mutated.ref.2'
ref3.id = 'mutated.ref.3'

with open('run.ref.mutated.fa', 'w') as f:
    print(mutated_ref1, ref2, ref3, sep='\n', file=f)

with open('run.false_call_vcf.vcf', 'w') as f:
    # put in an incorrect variant
    ref_at_900 = seq1[899]
    wrong_nucl = mutated_nucleotide(ref_at_900)
    line = 'ref.1\t900\t.\t' + ref_at_900 + '\t' + wrong_nucl + '\t42\t.\tDP=11;VDB=0.40105;SGB=-0.616816;RPB=0.279932;MQB=0.503877;BQB=1.00775;MQ0F=0.636364;AC=2;AN=2;DP4=4,0,6,0;MQ=12\tGT:PL\t1/1:85,6,0'
    print(line, file=f)

    # add a het call
    alts = {'A': 'C,G', 'C': 'A,T', 'G':'A,C', 'T':'A,C'}[seq1[949]]
    line = 'ref.1\t950\t.\t' + seq1[949] + '\t' + alts + '\t42\t.\tDP=11\tGT:PL\t1/2:85,6,0'
    print(line, file=f)

    # add a line with no genotype call
    line = 'ref.1\t975\t.\t' + seq1[974] + '\t' + mutated_nucleotide(seq1[974]) + '\t42\t.\tDP=11\tFOO\tBAR'
    print(line, file=f)


