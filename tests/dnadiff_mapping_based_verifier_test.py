import filecmp
import os
import unittest
import glob

import pyfastaq

from cluster_vcf_records import vcf_file_read

from minos import dnadiff_mapping_based_verifier

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "dnadiff_mapping_based_verifier")


class TestDnadiffMappingBasedVerifier(unittest.TestCase):
    def test_write_dnadiff_plus_flanks_to_fastas(self):
        dnadiff_file_in = os.path.join(data_dir, "test.snps")
        sample1_file_in = os.path.join(data_dir, "sample1.fa")
        sample2_file_in = os.path.join(data_dir, "sample2.fa")
        tmp_out1 = "tmp.write_dnadiff_plus_flanks_to_fastas.out.1.fa"
        tmp_out2 = "tmp.write_dnadiff_plus_flanks_to_fastas.out.2.fa"
        expected_out1 = os.path.join(data_dir, "sample1.plusflanks.fa")
        expected_out2 = os.path.join(data_dir, "sample2.plusflanks.fa")

        flank = 5
        dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier._write_dnadiff_plus_flanks_to_fastas(
            dnadiff_file_in, sample1_file_in, sample2_file_in, tmp_out1, tmp_out2, flank
        )

        self.assertTrue(filecmp.cmp(expected_out1, tmp_out1, shallow=False))
        self.assertTrue(filecmp.cmp(expected_out2, tmp_out2, shallow=False))
        os.unlink(tmp_out1)
        os.unlink(tmp_out2)

    def test_write_vars_plus_flanks_to_fasta1(self):
        vcfref_file_in = os.path.join(data_dir, "vcfref.fa")
        sample_file_in = os.path.join(data_dir, "sample1a.vcf")
        tmp_out = "tmp.write_vars_plus_flanks_to_fasta.out.1.fa"
        expected_out = os.path.join(data_dir, "sample1a.plusflanks.fa")

        vcf_header, vcf_records = vcf_file_read.vcf_file_to_dict(
            sample_file_in, sort=True, remove_useless_start_nucleotides=True
        )
        vcf_ref_seqs = {}
        pyfastaq.tasks.file_to_dict(vcfref_file_in, vcf_ref_seqs)
        flank = 5
        ns = 1
        dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier._write_vars_plus_flanks_to_fasta(
            tmp_out, vcf_records, vcf_ref_seqs, flank, ns
        )

        self.assertTrue(filecmp.cmp(expected_out, tmp_out, shallow=False))
        os.unlink(tmp_out)

    def test_write_vars_plus_flanks_to_fasta2(self):
        vcfref_file_in = os.path.join(data_dir, "vcfref.fa")
        sample_file_in = os.path.join(data_dir, "sample2a.vcf")
        tmp_out = "tmp.write_vars_plus_flanks_to_fasta.out.2.fa"
        expected_out = os.path.join(data_dir, "sample2a.plusflanks.fa")

        vcf_header, vcf_records = vcf_file_read.vcf_file_to_dict(
            sample_file_in, sort=True, remove_useless_start_nucleotides=True
        )
        vcf_ref_seqs = {}
        pyfastaq.tasks.file_to_dict(vcfref_file_in, vcf_ref_seqs)
        flank = 5
        ns = 1
        dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier._write_vars_plus_flanks_to_fasta(
            tmp_out, vcf_records, vcf_ref_seqs, flank, ns
        )

        self.assertTrue(filecmp.cmp(expected_out, tmp_out, shallow=False))
        os.unlink(tmp_out)

    def test_parse_sam_file_and_vcf1(self):
        samfile = os.path.join(data_dir, "sample1.sam")
        vcffile = os.path.join(data_dir, "sample1a.vcf")
        reffile = os.path.join(data_dir, "sample1.plusflanks.fa")

        flank = 5
        allow_mismatches = False
        dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier._index_vcf(vcffile)
        found, gt_conf, allele, match_flag, allele_flag = dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier._parse_sam_file_and_vcf(
            samfile, vcffile + ".gz", reffile, flank, allow_mismatches
        )

        exp_found = [
            "1",
            "1",
            "0",
            "0",
            "1",
            "0",
            "1",
        ]  # nb doesn't currently handle '.' alleles
        exp_gt_conf = [42, 42, 0, 0, 42, 0, 0]
        exp_allele = ["0", "1", "0", "0", "0", "0", "0"]
        self.assertEqual(exp_found, found)
        self.assertEqual(exp_gt_conf, gt_conf)
        self.assertEqual(exp_allele, allele)
        os.unlink(vcffile + ".gz")
        os.unlink(vcffile + ".gz.tbi")

    def test_parse_sam_file_and_vcf2(self):
        samfile = os.path.join(data_dir, "sample2.sam")
        vcffile = os.path.join(data_dir, "sample2a.vcf")
        reffile = os.path.join(data_dir, "sample2.plusflanks.fa")

        flank = 5
        allow_mismatches = False
        dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier._index_vcf(vcffile)
        found, gt_conf, allele, match_flag, allele_flag = dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier._parse_sam_file_and_vcf(
            samfile, vcffile + ".gz", reffile, flank, allow_mismatches
        )

        exp_found = ["1", "1", "1", "0", "1", "0", "1"]
        exp_gt_conf = [42, 42, 52, 0, 42, 0, 0]
        exp_allele = ["1", "0", "1", "0", "1", "0", "1"]
        self.assertEqual(exp_found, found)
        self.assertEqual(exp_gt_conf, gt_conf)
        self.assertEqual(exp_allele, allele)
        os.unlink(vcffile + ".gz")
        os.unlink(vcffile + ".gz.tbi")

    def test_run_with_filter_cluster_include_ref_alleles(self):
        """test run with filtering and clustering"""
        dnadiff_file_in = os.path.join(data_dir, "test.snps")
        sample1_file_in = os.path.join(data_dir, "sample1.fa")
        sample2_file_in = os.path.join(data_dir, "sample2.fa")
        vcf_file1_in = os.path.join(data_dir, "sample1a.vcf")
        vcf_file2_in = os.path.join(data_dir, "sample2a.vcf")
        vcf_reference_file = os.path.join(data_dir, "vcfref.fa")
        exp_out = os.path.join(data_dir, "exp_stats.tsv")
        exp_gt_conf = os.path.join(data_dir, "exp_gt_conf.tsv")
        exp_summary = os.path.join(data_dir, "exp_summary.tsv")

        tmp_out = "tmp.dnadiff_mapping_based_verifier.out"
        verifier = dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier(
            dnadiff_file_in,
            sample1_file_in,
            sample2_file_in,
            vcf_file1_in,
            vcf_file2_in,
            vcf_reference_file,
            tmp_out,
            flank_length=5,
            discard_ref_calls=False,
        )
        verifier.run()

        self.assertTrue(filecmp.cmp(exp_out, tmp_out + ".stats.tsv", shallow=False))
        self.assertTrue(
            filecmp.cmp(exp_gt_conf, tmp_out + ".gt_conf_hist.tsv", shallow=False)
        )
        self.assertTrue(
            filecmp.cmp(exp_summary, tmp_out + ".summary.tsv", shallow=False)
        )
        for f in glob.glob(tmp_out + "*"):
            os.unlink(f)
