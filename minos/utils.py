import subprocess
import sys

import pyfastaq
import pysam

class Error (Exception): pass

def syscall(command, allow_fail=False):
    completed_process = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    if (not allow_fail) and completed_process.returncode != 0:
        print('Error running this command:', command, file=sys.stderr)
        print('Return code:', completed_process.returncode, file=sys.stderr)
        print('\nOutput from stdout:', completed_process.stdout, sep='\n', file=sys.stderr)
        print('\nOutput from stderr:', completed_process.stderr, sep='\n', file=sys.stderr)
        raise Error('Error in system call. Cannot continue')

    return completed_process


def estimate_max_read_length_and_read_error_rate_from_qual_scores(infile, number_of_reads=10000):
    '''Estimates the maximum read length, and error rate from a file of reads, using the
    quality scores. Calculated by converting the mean phred quality score
    into the probability (formula is P = 10 ^ (- mean qual / 10)).
    Returns a tuple max_read_length, error_rate.
    Error rate is None if no quality
    scores found. Uses the first number_of_reads found in the file.
    File type can be BAM, SAM, FASTQ and is auto detected'''
    try:
        f = pysam.AlignmentFile(infile)
        filetype = 'pysam'
    except ValueError:
        f = pyfastaq.sequences.file_reader(infile)
        filetype = 'fastaq'
    except:
        raise Error('Error opening file ' + infile)

    total = 0
    base_count = 0
    read_count = 0
    max_read_length = 0

    for sequence in f:
        read_count += 1
        if read_count > number_of_reads:
            break

        if filetype == 'pysam':
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

