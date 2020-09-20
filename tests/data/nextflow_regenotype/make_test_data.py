#!/usr/bin/env python3

import copy
import os
import random
import subprocess

import pyfastaq

def make_ref_seqs():
    random.seed(42)
    ref_seq = [random.choice(["A", "C", "G", "T"]) for _ in range(2000)]

    ref_seq[49] = "A"
    ref_seq[99] = "G"
    ref_seq[499] = "A"
    ref_seq[999] = "C"
    ref_seq[1499] = "G"
    ref_seq[1504] = "G"
    ref_seq[1509] = "G"

    # this is to make sample 3 have too many snps, so it gets removed
    for i in range(1600, 1800, 20):
        ref_seq[i] = "T"

    sample_seqs = [
        copy.copy(ref_seq),
        copy.copy(ref_seq),
        copy.copy(ref_seq),
        copy.copy(ref_seq),
    ]

    # this is to make sample 3 have too many snps, so it gets removed
    for i in range(1600, 1800, 20):
        sample_seqs[3][i] = "A"

    sample_seqs[0][49] = "G"
    sample_seqs[0][99] = "C"
    sample_seqs[1][99] = "T"
    sample_seqs[2][99] = "GG"

    sample_seqs[0][499] = "T"
    sample_seqs[1][499] = "T"
    sample_seqs[2][499] = "T"

    sample_seqs[0][999] = "A"
    sample_seqs[1][999] = "G"
    sample_seqs[2][999] = "T"

    sample_seqs[0][1499] = "GCGCAGACGGCA"
    sample_seqs[2][1499] = "A"
    sample_seqs[1][1499] = "T"
    sample_seqs[1][1504] = "T"
    sample_seqs[1][1509] = "T"

    ref_seq.pop(1251)
    ref_seq = pyfastaq.sequences.Fasta("ref", "".join(ref_seq))
    sample_seqs = [pyfastaq.sequences.Fasta(f"sample.{i}", "".join(x)) for (i, x) in enumerate(sample_seqs)]
    return ref_seq, sample_seqs


def process_one_sample(ref_fa, sample_seq, outprefix):
    sample_fa = f"{outprefix}.fa"
    with open(sample_fa, "w") as f:
        print(sample_seq, file=f)

    reads_fq = f"{outprefix}.reads.fq"
    subprocess.check_output(f"fastaq to_perfect_reads {sample_fa} {reads_fq} 150 1 25 70", shell=True)
    bam = f"{outprefix}.bam"
    subprocess.check_output(f"bwa mem {ref_fa} {reads_fq} | samtools sort -O BAM -o {bam}", shell=True)
    subprocess.check_output(f"samtools index {bam}", shell=True)
    vcf = f"{outprefix}.vcf"
    subprocess.check_output(f"bcftools mpileup -f {ref_fa} {bam} | bcftools call -cv > {vcf}", shell=True)
    os.unlink(reads_fq)


ref_seq, sample_seqs = make_ref_seqs()
ref_fa = "data.ref.fa"
with open(ref_fa, "w") as f:
    print(ref_seq, file=f)
subprocess.check_output(f"bwa index {ref_fa}", shell=True)
subprocess.check_output(f"samtools faidx {ref_fa}", shell=True)

for seq in sample_seqs:
    print(seq.id)
    process_one_sample(ref_fa, seq, f"data.{seq.id}")


with open("data.sample.2.vcf") as f_in, open("data.sample.4.vcf", "w") as f_out:
    for line in f_in:
        if line.startswith("#CHROM"):
            line = line.replace("sample.2", "sample.4")
        print(line, end="", file=f_out)
