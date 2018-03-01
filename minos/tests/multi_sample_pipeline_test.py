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
        vcf3 = prefix + '.sample3.vcf'
        reads1 = prefix + '.sample1.reads'
        reads2 = prefix + '.sample2.reads'
        reads3_1 = prefix + '.sample3_1.reads'
        reads3_2 = prefix + '.sample3_2.reads'
        for filename in [tsv_file, vcf1, vcf2, reads1, reads2, reads3_1, reads3_2]:
            try:
                os.unlink(filename)
            except:
                pass

        with open(tsv_file, 'w') as f:
            print(vcf1, reads1, sep='\t', file=f)
            print(vcf2, reads2, sep='\t', file=f)
            print(vcf3, reads3_1, reads3_2, sep='\t', file=f)

        for filename in [vcf1, vcf2, vcf3, reads1, reads3_1, reads3_2]:
            with open(filename, 'w'):
                pass
            with self.assertRaises(multi_sample_pipeline.Error):
                multi_sample_pipeline.MultiSamplePipeline._load_input_data_tsv(tsv_file)

        with open(reads2, 'w'):
            pass

        got = multi_sample_pipeline.MultiSamplePipeline._load_input_data_tsv(tsv_file)
        expected = [
            (os.path.abspath(vcf1), [os.path.abspath(reads1)]),
            (os.path.abspath(vcf2), [os.path.abspath(reads2)]),
            (os.path.abspath(vcf3), [os.path.abspath(reads3_1), os.path.abspath(reads3_2)]),
        ]

        self.assertEqual(expected, got)

        for filename in [tsv_file, vcf1, vcf2, vcf3, reads1, reads2, reads3_1, reads3_2]:
            os.unlink(filename)


    def test_nextflow_helper_process_input_vcf_file(self):
        '''test _nextflow_helper_process_input_vcf_file'''
        infile = os.path.join(data_dir, 'nextflow_helper_process_input_vcf_file.in.vcf')
        expect_small = os.path.join(data_dir, 'nextflow_helper_process_input_vcf_file.out.small.vcf')
        expect_big = os.path.join(data_dir, 'nextflow_helper_process_input_vcf_file.out.big.vcf')
        expect_sample = os.path.join(data_dir, 'nextflow_helper_process_input_vcf_file.out.sample.txt')
        out_small = 'tmp.nextflow_helper_process_input_vcf_file.out.small.vcf'
        out_big = 'tmp.nextflow_helper_process_input_vcf_file.out.big.vcf'
        out_sample = 'tmp.nextflow_helper_process_input_vcf_file.out.sample.txt'
        multi_sample_pipeline.MultiSamplePipeline._nextflow_helper_process_input_vcf_file(infile, out_small, out_big, out_sample, 5)
        self.assertTrue(expect_small, out_small)
        self.assertTrue(expect_big, out_big)
        self.assertTrue(expect_sample, out_sample)
        os.unlink(out_small)
        os.unlink(out_big)
        os.unlink(out_sample)


    def test_write_nextflow_data_tsv(self):
        '''test _write_nextflow_data_tsv'''
        outfile = 'tmp.multi_sample_pipeline_test_write_nextflow_data_tsv'
        expected_file = os.path.join(data_dir, 'write_nextflow_data_tsv.expect.tsv')
        data = [
            ('vcf1', ['reads1']),
            ('vcf2', ['reads2.1', 'reads2.2']),
        ]

        multi_sample_pipeline.MultiSamplePipeline._write_nextflow_data_tsv(data, outfile)
        self.assertTrue(filecmp.cmp(expected_file, outfile, shallow=False))
        os.unlink(outfile)


    def test_prepare_nextflow_input_files(self):
        '''test _prepare_nextflow_input_files'''
        # Contents of the files is checked elsewhere.
        # We'll just check that the files exist
        ref_fasta = 'tmp.prepare_nextflow_input_files.in.ref.fa'
        outdir = 'tmp.prepare_nextflow_input_files.outdir'
        data_tsv = 'tmp.prepare_nextflow_input_files.in.tsv'
        vcf_file = 'tmp.prepare_nextflow_input_files.in.vcf'
        reads_file = 'tmp.prepare_nextflow_input_files.in.reads'
        with open(ref_fasta, 'w') as f:
            pass
        with open(data_tsv, 'w') as f:
            print(vcf_file, reads_file, sep='\t', file=f)
        with open(vcf_file, 'w'), open(reads_file, 'w'):
            pass
        if os.path.exists(outdir):
            shutil.rmtree(outdir)
        pipeline = multi_sample_pipeline.MultiSamplePipeline(ref_fasta, data_tsv, outdir)
        pipeline._make_output_dir()
        pipeline._prepare_nextflow_input_files()
        self.assertTrue(os.path.exists(outdir))
        self.assertTrue(os.path.exists(pipeline.nextflow_input_tsv))
        shutil.rmtree(outdir)
        os.unlink(ref_fasta)
        os.unlink(data_tsv)
        os.unlink(vcf_file)
        os.unlink(reads_file)


    def test_run(self):
        '''test run'''
        input_tsv = 'tmp.multi_sample_pipeline.run.in.tsv'
        ref_fasta = os.path.join(data_dir, 'run.ref.0.fa')
        with open(input_tsv, 'w') as f:
            for i in '1', '2':
                reads1 = os.path.join(data_dir, 'run.reads.' + i + '.1.fq')
                reads2 = os.path.join(data_dir, 'run.reads.' + i + '.2.fq')
                vcf = os.path.join(data_dir, 'run.calls.' + i + '.vcf')
                print(vcf, reads1, reads2, sep='\t', file=f)

        outdir = 'tmp.multi_sample_pipeline.run.out'

        pipeline = multi_sample_pipeline.MultiSamplePipeline(ref_fasta, input_tsv, outdir, min_large_ref_length=10)
        pipeline.run()

        shutil.rmtree(outdir)
        os.unlink(input_tsv)
