import filecmp
import os
import unittest

from minos import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "utils")


class TestUtils(unittest.TestCase):
    def test_fasta_to_upper_and_ACGT_only(self):
        """test fasta_to_upper_and_ACGT_only"""
        infile = os.path.join(data_dir, "fasta_to_upper_and_ACGT_only.in.fa")
        tmp_file = "tmp.fasta_to_upper_and_ACGT_only.fa"
        utils.rm_rf(tmp_file)
        utils.fasta_to_upper_and_ACGT_only(infile, tmp_file)
        expect = os.path.join(data_dir, "fasta_to_upper_and_ACGT_only.expect.fa")
        self.assertTrue(filecmp.cmp(tmp_file, expect, shallow=False))
        os.unlink(tmp_file)

    def test_estimate_max_read_length_and_read_error_rate_from_qual_scores_fastq_file(
        self,
    ):
        """test estimate_max_read_length_and_read_error_rate_from_qual_scores fastq file"""
        tmp_file = (
            "tmp.estimate_max_read_length_and_read_error_rate_from_qual_scores.fq"
        )
        with open(tmp_file, "w") as f:
            print("@1", "ACGT", "+", "IIHH", sep="\n", file=f)
            print("@2", "ACGTAG", "+", "IHGGFF", sep="\n", file=f)

        (
            got_length,
            got_qual,
        ) = utils.estimate_max_read_length_and_read_error_rate_from_qual_scores(
            tmp_file, number_of_reads=1
        )
        expect_qual = pow(10, -39.5 / 10)
        self.assertAlmostEqual(expect_qual, got_qual)
        self.assertEqual(4, got_length)
        (
            got_length,
            got_qual,
        ) = utils.estimate_max_read_length_and_read_error_rate_from_qual_scores(
            tmp_file, number_of_reads=2
        )
        expect_qual = pow(10, -38.7 / 10)
        self.assertAlmostEqual(expect_qual, got_qual)
        self.assertEqual(6, got_length)
        os.unlink(tmp_file)

    def test_estimate_max_read_length_and_read_error_rate_from_qual_scores_sam_file(
        self,
    ):
        """test estimate_max_read_length_and_read_error_rate_from_qual_scores sam file"""
        tmp_file = (
            "tmp.estimate_max_read_length_and_read_error_rate_from_qual_scores.sam"
        )
        with open(tmp_file, "w") as f:
            print("@SQ\tSN:ref\tLN:1000", file=f)
            print(
                1, 0, "ref", 42, 43, "4M", "*", 0, 0, "ACGT", "IIHH", sep="\t", file=f
            )
            print(
                2, 0, "ref", 42, 43, "4M", "*", 0, 0, "ACGT", "GGFF", sep="\t", file=f
            )

        (
            got_length,
            got_qual,
        ) = utils.estimate_max_read_length_and_read_error_rate_from_qual_scores(
            tmp_file, number_of_reads=1
        )
        expect_qual = pow(10, -39.5 / 10)
        self.assertAlmostEqual(expect_qual, got_qual)
        self.assertEqual(4, got_length)
        (
            got_length,
            got_qual,
        ) = utils.estimate_max_read_length_and_read_error_rate_from_qual_scores(
            tmp_file, number_of_reads=2
        )
        expect_qual = pow(10, -38.5 / 10)
        self.assertAlmostEqual(expect_qual, got_qual)
        self.assertEqual(4, got_length)
        os.unlink(tmp_file)

    def test_estimate_max_read_length_and_read_error_rate_from_qual_scores_sam_file_no_quals(
        self,
    ):
        """test estimate_max_read_length_and_read_error_rate_from_qual_scores sam file with no quals"""
        tmp_file = (
            "tmp.estimate_max_read_length_and_read_error_rate_from_qual_scores.sam"
        )
        with open(tmp_file, "w") as f:
            print("@SQ\tSN:ref\tLN:1000", file=f)
            print(1, 0, "ref", 42, 43, "4M", "*", 0, 0, "ACGT", "*", sep="\t", file=f)
            print(2, 0, "ref", 42, 43, "5M", "*", 0, 0, "ACGTA", "*", sep="\t", file=f)

        (
            got_length,
            got_qual,
        ) = utils.estimate_max_read_length_and_read_error_rate_from_qual_scores(
            tmp_file, number_of_reads=1
        )
        self.assertEqual(None, got_qual)
        self.assertEqual(4, got_length)
        (
            got_length,
            got_qual,
        ) = utils.estimate_max_read_length_and_read_error_rate_from_qual_scores(
            tmp_file, number_of_reads=2
        )
        self.assertEqual(None, got_qual)
        self.assertEqual(5, got_length)
        os.unlink(tmp_file)

    def test_estimate_max_read_length_and_read_error_rate_from_qual_scores_fasta_file(
        self,
    ):
        """test estimate_max_read_length_and_read_error_rate_from_qual_scores fasta file"""
        tmp_file = (
            "tmp.estimate_max_read_length_and_read_error_rate_from_qual_scores.fa"
        )
        with open(tmp_file, "w") as f:
            print(">1", "ACGT", sep="\n", file=f)
            print(">2", "ACGT", sep="\n", file=f)

        (
            got_length,
            got_qual,
        ) = utils.estimate_max_read_length_and_read_error_rate_from_qual_scores(
            tmp_file
        )
        self.assertEqual(None, got_qual)
        self.assertEqual(4, got_length)
        os.unlink(tmp_file)

    def test_remove_vars_from_vcf_at_contig_ends(self):
        ref_fa = os.path.join(data_dir, "remove_vars_from_vcf_at_contig_ends.ref.fa")
        vcf_in = os.path.join(data_dir, "remove_vars_from_vcf_at_contig_ends.in.vcf")
        expect_vcf = os.path.join(
            data_dir, "remove_vars_from_vcf_at_contig_ends.expect.vcf"
        )
        vcf_out = "tmp.remove_vars_from_vcf_at_contig_ends.vcf"
        utils.remove_vars_from_vcf_at_contig_ends(vcf_in, vcf_out, ref_fasta=ref_fa)
        self.assertTrue(filecmp.cmp(vcf_out, expect_vcf, shallow=False))
        utils.syscall(f"cp {vcf_in} {vcf_out}")
        ref_lengths = {"ref1": 4, "ref2": 2}
        utils.remove_vars_from_vcf_at_contig_ends(
            vcf_out, vcf_out, ref_lengths=ref_lengths
        )
