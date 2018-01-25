import os
import unittest

from minos import utils

modules_dir = os.path.dirname(os.path.abspath(utils.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'utils')

class TestUtils(unittest.TestCase):
    def test_estimate_read_error_rate_from_qual_scores_fastq_file(self):
        '''test estimate_read_error_rate_from_qual_scores fastq file'''
        tmp_file = 'tmp.estimate_read_error_rate_from_qual_scores.fq'
        with open(tmp_file, 'w') as f:
            print('@1', 'ACGT', '+', 'IIHH', sep='\n', file=f)
            print('@2', 'ACGT', '+', 'GGFF', sep='\n', file=f)

        got = utils.estimate_read_error_rate_from_qual_scores(tmp_file, number_of_reads=1)
        expect = pow(10, -39.5 / 10)
        self.assertAlmostEqual(expect, got)
        got = utils.estimate_read_error_rate_from_qual_scores(tmp_file, number_of_reads=2)
        expect = pow(10, -38.5 / 10)
        self.assertAlmostEqual(expect, got)
        os.unlink(tmp_file)


    def test_estimate_read_error_rate_from_qual_scores_sam_file(self):
        '''test estimate_read_error_rate_from_qual_scores sam file'''
        tmp_file = 'tmp.estimate_read_error_rate_from_qual_scores.sam'
        with open(tmp_file, 'w') as f:
            print('@SQ\tSN:ref\tLN:1000', file=f)
            print(1, 0, 'ref', 42, 43, '4M', '*', 0, 0, 'ACGT', 'IIHH', sep='\t', file=f)
            print(2, 0, 'ref', 42, 43, '4M', '*', 0, 0, 'ACGT', 'GGFF', sep='\t', file=f)

        got = utils.estimate_read_error_rate_from_qual_scores(tmp_file, number_of_reads=1)
        expect = pow(10, -39.5 / 10)
        self.assertAlmostEqual(expect, got)
        got = utils.estimate_read_error_rate_from_qual_scores(tmp_file, number_of_reads=2)
        expect = pow(10, -38.5 / 10)
        self.assertAlmostEqual(expect, got)
        os.unlink(tmp_file)


    def test_estimate_read_error_rate_from_qual_scores_sam_file_no_quals(self):
        '''test estimate_read_error_rate_from_qual_scores sam file with no quals'''
        tmp_file = 'tmp.estimate_read_error_rate_from_qual_scores.sam'
        with open(tmp_file, 'w') as f:
            print('@SQ\tSN:ref\tLN:1000', file=f)
            print(1, 0, 'ref', 42, 43, '4M', '*', 0, 0, 'ACGT', '*', sep='\t', file=f)
            print(2, 0, 'ref', 42, 43, '4M', '*', 0, 0, 'ACGT', '*', sep='\t', file=f)

        got = utils.estimate_read_error_rate_from_qual_scores(tmp_file, number_of_reads=1)
        self.assertEqual(None, got)
        got = utils.estimate_read_error_rate_from_qual_scores(tmp_file, number_of_reads=2)
        self.assertEqual(None, got)
        os.unlink(tmp_file)


    def test_estimate_read_error_rate_from_qual_scores_fasta_file(self):
        '''test estimate_read_error_rate_from_qual_scores fasta file'''
        tmp_file = 'tmp.estimate_read_error_rate_from_qual_scores.fa'
        with open(tmp_file, 'w') as f:
            print('>1', 'ACGT', sep='\n', file=f)
            print('>2', 'ACGT', sep='\n', file=f)

        got = utils.estimate_read_error_rate_from_qual_scores(tmp_file)
        self.assertEqual(None, got)
        os.unlink(tmp_file)

