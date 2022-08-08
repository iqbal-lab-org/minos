# If this is run as a script, then it makes a directory `test_data_files`
# with all the test data files.
# This is not indended to be run by the user.
# It assumes samtools and minimap2 are installed.
#
# There's a separate commnd line task make_test_data. It copies the
# Files directory, but then also edits the file `manifest.tsv` (used by the
# nextflow joint genotype piepline), to put absolute paths in there.
# Has to be done at run time because that's when we know the full path.
import copy
import os
import random
import subprocess
import pyfastaq

from minos import utils

def _call_variants(ref, reads, out, mates=None):
    bam = f"{out}.bam"
    if mates is None:
        command = f"minimap2 -a {ref} {reads}"
    else:
        command = f"minimap2 -x sr -a {ref} {reads} {mates}"
    command += f" | samtools sort -O BAM -o {bam}"
    subprocess.check_output(command, shell=True)
    command = f"samtools mpileup -ugf {ref} {bam} | bcftools call -vm -O v -o {out}.vcf"
    subprocess.check_output(command, shell=True)
    subprocess.check_output(f"samtools index {bam}", shell=True)


def _make_test_data():
    random.seed(42)

    outdir = "test_data_files"
    os.mkdir(outdir)

    ref_fasta = os.path.join(outdir, "ref.fa")
    sample1_pre = os.path.join(outdir, "sample1")
    sample2_pre = os.path.join(outdir, "sample2")
    reads1_1 = f"{sample1_pre}.reads_1.fastq.gz"
    reads1_2 = f"{sample1_pre}.reads_2.fastq.gz"
    reads2 = f"{sample2_pre}.reads.fastq.gz"

    ref_seq = random.choices(["A", "C", "G", "T"], k=1000)
    mut_seq1 = copy.copy(ref_seq)
    mut_seq2 = copy.copy(ref_seq)

    ref_seq[199] = "A"
    mut_seq1[199] = "C"
    mut_seq2[199] = "C"

    ref_seq[399] = "G"
    mut_seq1[399] = "A"
    mut_seq2[399] = "G"

    ref_seq[599] = "T"
    mut_seq1[599] = "T"
    mut_seq2[599] = "C"

    with open(ref_fasta, "w") as f:
        print(">ref", "".join(ref_seq), sep="\n", file=f)

    mut_fasta1 = os.path.join(outdir, "tmp.mut1.fa")
    with open(mut_fasta1, "w") as f:
        print(">m1", "".join(mut_seq1), sep="\n", file=f)

    mut_fasta2 = os.path.join(outdir, "tmp.mut2.fa")
    with open(mut_fasta2, "w") as f:
        print(">m2", "".join(mut_seq2), sep="\n", file=f)

    subprocess.check_output(f"fastaq to_perfect_reads {mut_fasta1} - 250 1 20 100 | fastaq deinterleave - {reads1_1} {reads1_2}", shell=True)

    subprocess.check_output(f"fastaq to_perfect_reads {mut_fasta2} {reads2} 250 1 20 100", shell=True)

    call_variants(ref_fasta, reads1_1, sample1_pre, mates=reads1_2)
    call_variants(ref_fasta, reads2, sample2_pre)

    os.unlink(mut_fasta1)
    os.unlink(mut_fasta2)

    os.rename(f"{sample1_pre}.vcf", os.path.join(outdir, "in.1.vcf"))
    os.rename(f"{sample2_pre}.vcf", os.path.join(outdir, "in.2.vcf"))

    with open(os.path.join(outdir, "manifest.tsv"), "w") as f:
        print("name\treads\tvcf", file=f)
        print("sample1\tsample1.bam\tin.1.vcf", file=f)
        print("sample2\tsample2.bam\tin.2.vcf", file=f)


def make_test_data_dir(outdir):
    if os.path.exists(outdir):
        raise Exception(f"Output directory {outdir} already exists. Cannot continue")
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    dir_to_copy = os.path.join(this_file_dir, "test_data_files")
    if not os.path.exists(dir_to_copy):
        raise Exception(f"Error! Did not find test data files. Expected to find them here: {dir_to_copy}")

    utils.syscall(f"cp -rp {dir_to_copy} {outdir}")

    # Need to fix the manifest file so that the filenames are absolute paths,
    # required by the joint genotype pipeline
    outdir = os.path.abspath(outdir)
    manifest = os.path.join(outdir, "manifest.tsv")
    lines_out = []
    with open(manifest) as f:
        for line in f:
            if len(lines_out) == 0:
                lines_out.append(line.rstrip().split("\t"))
            else:
                fields = line.rstrip().split("\t")
                for i in [1,2]:
                    fields[i] = os.path.join(outdir, fields[i])
                lines_out.append(fields)

    with open(manifest, "w") as f:
        for t in lines_out:
            print(*t, sep="\t", file=f)


if __name__ == "__main__":
    _make_test_data()

