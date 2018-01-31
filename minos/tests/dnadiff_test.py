import shutil
import os
import unittest

from minos import dnadiff

modules_dir = os.path.dirname(os.path.abspath(dnadiff.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'dnadiff')

class TestDnadiff(unittest.TestCase):
    def test_run_dnadiff(self):
        # Just test that dnadiff runs and we get the
        # output files we need. Not check their contents.
        ref_fasta = os.path.join(data_dir, 'run_dnadiff.ref.fa')
        qry_fasta = os.path.join(data_dir, 'run_dnadiff.qry.fa')
        tmp_prefix = 'tmp.run_dnadiff'
        dnadiff.Dnadiff.run_dnadiff(ref_fasta, qry_fasta, tmp_prefix)

        expected_extensions = [
            '1coords',
            '1delta',
            'delta',
            'mcoords',
            'mdelta',
            'qdiff',
            'rdiff',
            'report',
            'snps',
        ]
        expected_files = [tmp_prefix + '.' + x for x in expected_extensions]

        for filename in expected_files:
            self.assertTrue(os.path.exists(filename))
            os.unlink(filename)
