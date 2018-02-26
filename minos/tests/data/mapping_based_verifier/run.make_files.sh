#!/usr/bin/env bash
set -e

./run.make_refs.py
fastaq to_perfect_reads --seed 42 run.ref.mutated.fa - 230 1 20 100 | fastaq deinterleave - run.reads_1.fq run.reads_2.fq
bwa index run.ref.fa
bwa index run.ref.mutated.fa
bwa mem run.ref.fa run.reads_1.fq run.reads_2.fq | samtools sort -o run.bwa.bam
samtools mpileup -ugf run.ref.fa run.bwa.bam | bcftools call -vm -O v -o run.calls.vcf
cat run.false_call_vcf.vcf >> run.calls.vcf

rm run.reads* run.bwa.bam run.false_call_vcf.vcf

# samtools gets one of the variants wrong. Looks like this:
# ref ACGCTA-CCGCAACG
# alt ACGCTAGTCGCAACG
# at position 700. Add in correct allele, but
# leave the old one so hits has >1 ALT, which
# we also want to test.
sed -i 's/700\t.\tA\tAG/700\t.\tAC\tAGC,AGT/' run.calls.vcf

# Add GT_CONF
awk -F"\t" '/^#/ {print} !/^#/ {OFS="\t"; $9=$9":GT_CONF"; $10=$10":"NR; print $0}' run.calls.vcf > run.calls.vcf.tmp.$$
mv run.calls.vcf.tmp.$$ run.calls.vcf

