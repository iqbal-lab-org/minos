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

    def test_parse_sam_file_and_vcf1(self):
        samfile = os.path.join(data_dir, 'sample1.sam')
        vcffile = os.path.join(data_dir, 'sample1a.vcf')

        flank = 5
        allow_mismatches = False
        found, gt_conf = dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier._parse_sam_file_and_vcf(samfile, vcffile, flank, allow_mismatches)

        exp_found = ['1','1','0','0','1','1','1']
        exp_gt_conf = [42, 42, None, None, 42, 32, None]
        self.assertEqual(exp_found, found)
        self.assertEqual(exp_gt_conf, gt_conf)

    def test_parse_sam_file_and_vcf2(self):
        samfile = os.path.join(data_dir, 'sample2.sam')
        vcffile = os.path.join(data_dir, 'sample2a.vcf')

        flank = 5
        allow_mismatches = False
        found, gt_conf = dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier._parse_sam_file_and_vcf(samfile, vcffile, flank, allow_mismatches)

        exp_found = ['1','1','1','0','1','0','1']
        exp_gt_conf = [42, 42, 52, None, 42, None, None]
        self.assertEqual(exp_found, found)
        self.assertEqual(exp_gt_conf, gt_conf)

    def test_run_with_filter_cluster_include_ref_alleles(self):
        '''test run with filtering and clustering'''
        dnadiff_file_in = os.path.join(data_dir, 'test.snps')
        sample1_file_in = os.path.join(data_dir, 'sample1.fa')
        sample2_file_in = os.path.join(data_dir, 'sample2.fa')
        vcf_file1_in = os.path.join(data_dir, 'sample1a.vcf')
        vcf_file2_in = os.path.join(data_dir, 'sample2a.vcf')
        vcf_reference_file = os.path.join(data_dir, 'vcfref.fa')

        tmp_out = 'tmp.dnadiff_mapping_based_verifier.out'
        verifier = dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier(dnadiff_file_in, sample1_file_in, sample2_file_in, vcf_file1_in, vcf_file2_in, vcf_reference_file, tmp_out, flank_length=5, discard_ref_calls=False)
        verifier.run()



