import os
import unittest

from minos import bam_read_extract

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "bam_read_extract")


def read_names_match(bam1, bam2):
    names1 = bam_read_extract.get_read_names(bam1)
    names2 = bam_read_extract.get_read_names(bam2)
    return names1 == names2


class TestBamReadExtract(unittest.TestCase):
    def test_get_readnames(self):
        """test get_readnames"""
        infile = os.path.join(data_dir, "all_reads.bam")
        expected = ["read." + str(i) for i in [5, 0, 1, 2, 3, 4]]
        got = bam_read_extract.get_read_names(infile)
        self.assertEqual(expected, got)

    def test_get_unmapped_reads(self):
        """test get_unmapped_reads"""
        infile = os.path.join(data_dir, "all_reads.bam")
        expected_bam = os.path.join(data_dir, "unmapped_reads.bam")
        tmp_out = "tmp.bam_read_extract.unmapped.bam"
        bam_read_extract.get_unmapped_reads(infile, tmp_out)
        self.assertTrue(read_names_match(expected_bam, tmp_out))
        os.unlink(tmp_out)

    def test_get_region(self):
        """test get_region"""
        infile = os.path.join(data_dir, "all_reads.bam")
        tmp_out = "tmp.bam_read_extract.get_region.bam"

        expected_bam = os.path.join(data_dir, "region.1.60-181.bam")
        bam_read_extract.get_region(infile, "1", 59, 180, tmp_out)
        self.assertTrue(read_names_match(expected_bam, tmp_out))
        os.unlink(tmp_out)

        expected_bam = os.path.join(data_dir, "region.1.61-180.bam")
        bam_read_extract.get_region(infile, "1", 60, 179, tmp_out)
        self.assertTrue(read_names_match(expected_bam, tmp_out))
        os.unlink(tmp_out)
