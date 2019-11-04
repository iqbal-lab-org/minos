#!/usr/bin/env bash
set -e

./run_gramtools.make_refs.py
pre=run_gramtools
fastaq to_perfect_reads --seed 42 $pre.ref.mutated.fa $pre.reads.fq  230 1 20 100
fastaq deinterleave $pre.reads.fq $pre.reads_1.fq $pre.reads_2.fq
bwa index $pre.ref.fa
bwa index $pre.ref.mutated.fa
bwa mem $pre.ref.fa $pre.reads_1.fq $pre.reads_2.fq | samtools sort -o $pre.bwa.bam
samtools mpileup -ugf $pre.ref.fa $pre.bwa.bam | bcftools call -vm -O v -o $pre.calls.vcf
cat $pre.false_call_vcf.vcf >> $pre.calls.vcf

#rm $pre.bwa.bam $pre.false_call_vcf.vcf

# samtools gets one of the variants wrong. Looks like this:
# ref ACGCTA-CCGCAACG
# alt ACGCTAGTCGCAACG
# at position 700. Fix it:
sed -i 's/700\t.\tA\tAG/700\t.\tAC\tAGT/' $pre.calls.vcf

# gramtools ignores lines that don't PASS the filter.
# Fix the VCF so every line passes
awk -F"\t" '$1!~/^#/ {OFS="\t"; $7="PASS"} 1' $pre.calls.vcf > $pre.calls.vcf.$$
mv $pre.calls.vcf.$$ $pre.calls.vcf
