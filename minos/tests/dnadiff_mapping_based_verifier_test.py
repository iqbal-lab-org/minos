import filecmp
import os
import unittest

from cluster_vcf_records import vcf_record
import pysam
import pyfastaq

from minos import dnadiff_mapping_based_verifier

modules_dir = os.path.dirname(os.path.abspath(dnadiff_mapping_based_verifier.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'dnadiff_mapping_based_verifier')


class TestDnadiffMappingBasedVerifier(unittest.TestCase):
    def test_write_dnadiff_plus_flanks_to_fastas(self):
        dnadiff_file_in = os.path.join(data_dir, 'test.snps')
        sample1_file_in = os.path.join(data_dir, 'sample1.fa')
        sample2_file_in = os.path.join(data_dir, 'sample2.fa')
        tmp_out1 = 'tmp.write_dnadiff_plus_flanks_to_fastas.out.1.fa'
        tmp_out2 = 'tmp.write_dnadiff_plus_flanks_to_fastas.out.2.fa'
        expected_out1 = os.path.join(data_dir, 'sample1.plusflanks.fa')
        expected_out2 = os.path.join(data_dir, 'sample2.plusflanks.fa')

        flank = 5
        dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier._write_dnadiff_plus_flanks_to_fastas(dnadiff_file_in,sample1_file_in, sample2_file_in, tmp_out1, tmp_out2, flank)

        self.assertTrue(filecmp.cmp(expected_out1, tmp_out1, shallow=False))
        self.assertTrue(filecmp.cmp(expected_out2, tmp_out2, shallow=False))
        os.unlink(tmp_out1)
        os.unlink(tmp_out2)
