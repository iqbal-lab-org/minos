#!/usr/bin/env python3

import copy
import random
import pyfastaq

random.seed(1)

def mutated_nucleotide(nucleotide):
    other_nucleotides = ['A', 'C', 'G', 'T']
    other_nucleotides.remove(nucleotide)
    return random.choice(other_nucleotides)


nucleotides = ['A', 'C', 'G', 'T']
seq1 = [random.choice(nucleotides) for _ in range(1000)]
seq2 = [random.choice(nucleotides) for _ in range(1000)]
seq3 = seq1[:300] + [random.choice(nucleotides) for _ in range(700)]
seq1[99] = 'G'
seq1[50] = 'n'

ref1 = pyfastaq.sequences.Fasta('ref.1', ''.join(seq1))
ref2 = pyfastaq.sequences.Fasta('ref.2', ''.join(seq2))
ref3 = pyfastaq.sequences.Fasta('ref.3', ''.join(seq3))


mutated_seq1 = copy.copy(seq1)
mutated_seq1[99] = 'T'
mutated_seq1.pop(299)
mutated_seq1.pop(299)
mutated_seq1.insert(499, 'AGTC')
to_delete = mutated_seq1[699:703]
insertion = [mutated_nucleotide(x) for x in to_delete]
mutated_seq1 = mutated_seq1[:699] + insertion + mutated_seq1[703:]

mutated_ref1 = pyfastaq.sequences.Fasta('mutated.ref.1', ''.join(mutated_seq1))

with open('run.ref.fa', 'w') as f:
    print(ref1, ref2, ref3, sep='\n', file=f)


ref2.id = 'mutated.ref.2'
ref3.id = 'mutated.ref.3'

with open('run.ref.mutated.fa', 'w') as f:
    print(mutated_ref1, ref2, ref3, sep='\n', file=f)

with open('run.false_call_vcf.vcf', 'w') as f:
    ref_at_900 = seq1[899]
    wrong_nucl = mutated_nucleotide(ref_at_900)
    line = 'ref.1\t900\t.\t' + ref_at_900 + '\t' + wrong_nucl + '\t42\t.\tDP=11;VDB=0.40105;SGB=-0.616816;RPB=0.279932;MQB=0.503877;BQB=1.00775;MQ0F=0.636364;AC=2;AN=2;DP4=4,0,6,0;MQ=12\tGT:PL\t1/1:85,6,0'
    print(line, file=f)

