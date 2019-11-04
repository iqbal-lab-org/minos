import filecmp
import os
import unittest

from minos import vcf_file_split_deletions

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "vcf_file_split_deletions")


class TestVcfFileSplitDeletions(unittest.TestCase):
    def test_run(self):
        """test run"""
        infile = os.path.join(data_dir, "run.in.vcf")
        expect_5_small = os.path.join(data_dir, "run.out.cutoff5.small.vcf")
        expect_5_big = os.path.join(data_dir, "run.out.cutoff5.big.vcf")
        expect_10_small = os.path.join(data_dir, "run.out.cutoff10.small.vcf")
        expect_10_big = os.path.join(data_dir, "run.out.cutoff10.big.vcf")
        out_small = "tmp.vcf_file_split_deletions.out.small.vcf"
        out_big = "tmp.vcf_file_split_deletions.out.big.vcf"

        splitter = vcf_file_split_deletions.VcfFileSplitDeletions(
            infile, out_small, out_big, min_large_ref_length=5
        )
        splitter.run()
        self.assertTrue(filecmp.cmp(expect_5_small, out_small, shallow=False))
        self.assertTrue(filecmp.cmp(expect_5_big, out_big, shallow=False))
        os.unlink(out_small)
        os.unlink(out_big)

        splitter = vcf_file_split_deletions.VcfFileSplitDeletions(
            infile, out_small, out_big, min_large_ref_length=10
        )
        splitter.run()
        self.assertTrue(filecmp.cmp(expect_10_small, out_small, shallow=False))
        self.assertTrue(filecmp.cmp(expect_10_big, out_big, shallow=False))
        os.unlink(out_small)
        os.unlink(out_big)
