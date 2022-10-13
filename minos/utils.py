import subprocess
import sys

import pyfastaq
import pysam


def fasta_to_upper_and_ACGT_only(infile, outfile):
    f = pyfastaq.sequences.file_reader(infile)
    with open(outfile, "w") as f_out:
        nucs = {"A", "C", "G", "T"}
        for seq in f:
            new_seq = list(seq.seq)
            for i, c in enumerate(new_seq):
                new_seq[i] = c.upper()
                if new_seq[i] not in nucs:
                    new_seq[i] = "C"
            seq.seq = "".join(new_seq)
            print(seq, file=f_out)


def rm_rf(*paths):
    for path in paths:
        subprocess.check_output(f"rm -rf {path}", shell=True)


def syscall(command, allow_fail=False, cwd=None):
    completed_process = subprocess.run(
        command,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        cwd=cwd,
    )
    if (not allow_fail) and completed_process.returncode != 0:
        print("Error running this command:", command, file=sys.stderr)
        print("Return code:", completed_process.returncode, file=sys.stderr)
        print(
            "\nOutput from stdout:", completed_process.stdout, sep="\n", file=sys.stderr
        )
        print(
            "\nOutput from stderr:", completed_process.stderr, sep="\n", file=sys.stderr
        )
        raise Exception("Error in system call. Cannot continue")

    return completed_process


def estimate_max_read_length_and_read_error_rate_from_qual_scores(
    infile, number_of_reads=10000
):
    """Estimates the maximum read length, and error rate from a file of reads, using the
    quality scores. Calculated by converting the mean phred quality score
    into the probability (formula is P = 10 ^ (- mean qual / 10)).
    Returns a tuple max_read_length, error_rate.
    Error rate is None if no quality
    scores found. Uses the first number_of_reads found in the file.
    File type can be BAM, SAM, FASTQ and is auto detected"""
    try:
        f = pysam.AlignmentFile(infile)
        filetype = "pysam"
    except ValueError:
        f = pyfastaq.sequences.file_reader(infile)
        filetype = "fastaq"
    except:
        raise Exception("Error opening file " + infile)

    total = 0
    base_count = 0
    read_count = 0
    max_read_length = 0

    for sequence in f:
        read_count += 1
        if read_count > number_of_reads:
            break

        if filetype == "pysam":
            max_read_length = max(max_read_length, sequence.query_length)
            if sequence.query_qualities is not None:
                for qual_score in sequence.query_qualities:
                    total += qual_score
                base_count += len(sequence.query_qualities)
        else:
            max_read_length = max(max_read_length, len(sequence))
            if type(sequence) != pyfastaq.sequences.Fastq:
                continue

            for qual_score in [ord(x) - 33 for x in sequence.qual]:
                total += qual_score
            base_count += len(sequence)

    if base_count == 0:
        return max_read_length, None

    mean_qual_score = total / base_count
    return max_read_length, pow(10, -mean_qual_score / 10)


def remove_vars_from_vcf_at_contig_ends(
    vcf_in, vcf_out, ref_lengths=None, ref_fasta=None
):
    if ref_lengths is None:
        assert ref_fasta is not None
        ref_lengths = {
            x.id.split()[0]: len(x) for x in pyfastaq.sequences.file_reader(ref_fasta)
        }

    lines = []
    with open(vcf_in) as f:
        for line in map(str.rstrip, f):
            if not line.startswith("#"):
                fields = line.split("\t")
                ref_length = ref_lengths[fields[0]]
                ref_end = int(fields[1]) + len(fields[3]) - 1
                if ref_end >= ref_length:
                    continue
            lines.append(line)

    with open(vcf_out, "w") as f:
        print(*lines, sep="\n", file=f)
