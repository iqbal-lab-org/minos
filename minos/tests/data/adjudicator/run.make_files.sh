#!/usr/bin/env bash
set -e

./run.make_refs.py
pre=run
fastaq to_perfect_reads --seed 42 $pre.ref.mutated.fa $pre.reads.fq  230 1 20 100
fastaq deinterleave $pre.reads.fq $pre.reads_1.fq $pre.reads_2.fq
bwa index $pre.ref.fa
bwa index $pre.ref.mutated.fa
bwa mem $pre.ref.fa $pre.reads_1.fq $pre.reads_2.fq | samtools sort -o $pre.bwa.bam
samtools mpileup -ugf $pre.ref.fa $pre.bwa.bam | bcftools call -vm -O v -o $pre.calls.1.vcf
cat $pre.false_call_vcf.vcf >> $pre.calls.1.vcf

#rm $pre.bwa.bam $pre.false_call_vcf.vcf

# make a second vcf with a couple of differences
awk -F"\t" '$2!=100' $pre.calls.1.vcf | awk -F"\t" 'BEGIN{OFS="\t"} {if($2==299) {$5="GA,CC"}; $1=$1; print $0}' > $pre.calls.2.vcf
