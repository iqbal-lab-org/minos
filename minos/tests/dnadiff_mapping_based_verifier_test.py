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

    def test_parse_sam_file_and_vcf(self):
        samfile = os.path.join(data_dir, 'test.1.sam')
        vcffile = os.path.join(data_dir, 'sample1a.vcf')

        flank = 5
        allow_mismatches = False
        found, gt_conf = dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier._parse_sam_file_and_vcf(samfile, vcffile, flank, allow_mismatches)

        exp_found = ['1','1','0','0','1','1','1']
        exp_gt_conf = [42.42, 42.42, None, None, 42.42, 32.32, None]
        self.assertEqual(expect_found, found)
        self.assertEqual(expect_gt_conf, gt_conf)


