import filecmp
import os
import shutil
import unittest

from minos import multi_sample_pipeline

modules_dir = os.path.dirname(os.path.abspath(multi_sample_pipeline.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'multi_sample_pipeline')

class TestMultiSamplePipeline(unittest.TestCase):
    def test_load_input_data_tsv(self):
        '''test _load_input_data_tsv'''
        # Using absolute paths, so can't easily have
        # pre-made input files. Hence make them now
        prefix = 'tmp.multi_sample_pipeline_load_input_data_tsv'
        tsv_file = prefix + '.data.tsv'
        vcf1 = prefix + '.sample1.vcf'
        vcf2 = prefix + '.sample2.vcf'
        reads1 = prefix + '.sample1.reads'
        reads2 = prefix + '.sample2.reads'
        for filename in [tsv_file, vcf1, vcf2, reads1, reads2]:
            try:
                os.unlink(filename)
            except:
                pass

        with open(tsv_file, 'w') as f:
            print(vcf1, reads1, sep='\t', file=f)
            print(vcf2, reads2, sep='\t', file=f)

        for filename in [vcf1, vcf2, reads1]:
            with open(filename, 'w'):
                pass
            with self.assertRaises(multi_sample_pipeline.Error):
                multi_sample_pipeline.MultiSamplePipeline._load_input_data_tsv(tsv_file)

        with open(reads2, 'w'):
            pass

        got = multi_sample_pipeline.MultiSamplePipeline._load_input_data_tsv(tsv_file)
        expected = [
            (os.path.abspath(vcf1), os.path.abspath(reads1)),
            (os.path.abspath(vcf2), os.path.abspath(reads2)),
        ]

        self.assertEqual(expected, got)

        for filename in [tsv_file, vcf1, vcf2, reads1, reads2]:
            os.unlink(filename)


    def test_write_nextflow_data_tsv(self):
        '''test _write_nextflow_data_tsv'''
        outfile = 'tmp.multi_sample_pipeline_test_write_nextflow_data_tsv'
        expected_file = os.path.join(data_dir, 'write_nextflow_data_tsv.expect.tsv')
        data = [
            ('vcf1', 'reads1'),
            ('vcf2', 'reads2'),
        ]

        multi_sample_pipeline.MultiSamplePipeline._write_nextflow_data_tsv(data, outfile)
        self.assertTrue(filecmp.cmp(expected_file, outfile, shallow=False))
        os.unlink(outfile)

