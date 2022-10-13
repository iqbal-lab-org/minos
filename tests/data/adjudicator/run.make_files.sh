#!/usr/bin/env bash
set -e

./run.make_refs.py
pre=run
fastaq to_perfect_reads --seed 42 $pre.ref.mutated.fa $pre.reads.fq  230 1 20 100
fastaq deinterleave $pre.reads.fq $pre.reads_1.fq $pre.reads_2.fq

# add in another ref sequence here, so it has no reads mapped.
# This is to test when the genotyper gets no info
echo ">ref.4
ACGTACTGCATGCATGCATCGTCTGACATATTCCGTACGTACGTACGTACGTGTCTCGAC
CATCGACTTGACTGACTGACTGCATCTGTCGATCGCTATGCATGACTGCATTGATATGCT
CGTACTGACGTAACTCTACATCATCATCTGACTGACTGATGCATCTGACTCTACTGACTG" >> $pre.ref.fa

bwa index $pre.ref.fa
bwa index $pre.ref.mutated.fa
bwa mem $pre.ref.fa $pre.reads_1.fq $pre.reads_2.fq | samtools sort -o $pre.bwa.bam
samtools index $pre.bwa.bam
samtools mpileup -ugf $pre.ref.fa $pre.bwa.bam | bcftools call -vm -O v -o $pre.calls.1.vcf
cat $pre.false_call_vcf.vcf >> $pre.calls.1.vcf

#rm $pre.bwa.bam $pre.false_call_vcf.vcf

# make a second vcf with a couple of differences
awk -F"\t" '$2!=100' $pre.calls.1.vcf | awk -F"\t" 'BEGIN{OFS="\t"} {if($2==299) {$5="GA,CC"}; $1=$1; print $0}' > $pre.calls.2.vcf

# add a variant to the reference that has no reads
echo -e "ref.4\t61\t.\tC\tA,T\t.\t.\tDP=42\tGT\t1/2" >> run.calls.1.vcf

# add negative coord to vcf 1
grep '^#' $pre.calls.1.vcf > $pre.calls.1.vcf.$$
echo -e "ref.1\t-1\t.\tA\tG\t61\t.\tDP=7;VDB=0.40105;SGB=-0.616816;MQ0F=0.428571;AC=2;AN=2;DP4=0,0,6,0;MQ=20\tGT:PL\t1/1:88,18,0"  >> $pre.calls.1.vcf.$$
grep -v '^#' $pre.calls.1.vcf >> $pre.calls.1.vcf.$$

# change some calls to het instead of hom and add a non-ACGT allele
awk -F"\t" 'BEGIN{OFS="\t"} $2==703 {$5="T,G"; $10="0/1"substr($10,4)} $2==704 {$5="T,<FOO>"} 1' $pre.calls.1.vcf.$$ > $pre.calls.1.vcf
rm $pre.calls.1.vcf.$$

# add a variant at the final coord of ref1 - should get removed before gramtools
# because it makes gramtools build crash
echo -e "ref.1\t1000\t.\tC\tA\t.\t.\tDP=100\tGT\t1/1" >> run.calls.1.vcf
