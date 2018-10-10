import shutil
import os
import unittest

from minos import adjudicator

modules_dir = os.path.dirname(os.path.abspath(adjudicator.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'adjudicator')

class TestAdjudicator(unittest.TestCase):
    def test_get_gramtools_kmer_size(self):
        '''test _get_gramtools_kmer_size'''
        build_dir = os.path.join(data_dir, 'get_gramtools_kmer_size.build')
        self.assertEqual(42, adjudicator.Adjudicator._get_gramtools_kmer_size(build_dir, None))
        self.assertEqual(42, adjudicator.Adjudicator._get_gramtools_kmer_size(build_dir, 20))
        self.assertEqual(20, adjudicator.Adjudicator._get_gramtools_kmer_size(None, 20))
        self.assertEqual(10, adjudicator.Adjudicator._get_gramtools_kmer_size(None, None))


    def test_run(self):
        '''test run'''
        # We're just testing that it doesn't crash.
        # Check the output files exist, but not their contents.
        # First run using splitting of VCF file.
        # Then run without splitting.
        outdir = 'tmp.adjudicator.out'
        if os.path.exists(outdir):
            shutil.rmtree(outdir)
        ref_fasta = os.path.join(data_dir, 'run.ref.fa')
        reads_file = os.path.join(data_dir, 'run.bwa.bam')
        vcf_files =  [os.path.join(data_dir, x) for x in ['run.calls.1.vcf', 'run.calls.2.vcf']]
        adj = adjudicator.Adjudicator(outdir, ref_fasta, [reads_file], vcf_files, variants_per_split=3, clean=False, gramtools_kmer_size=5)
        adj.run()
        self.assertTrue(os.path.exists(outdir))
        self.assertTrue(os.path.exists(adj.log_file))
        self.assertTrue(os.path.exists(adj.final_vcf))
        self.assertTrue(os.path.exists(adj.clustered_vcf))

        # Clean up and then run without splitting
        shutil.rmtree(outdir)
        adj = adjudicator.Adjudicator(outdir, ref_fasta, [reads_file], vcf_files, clean=False, gramtools_kmer_size=5)
        adj.run()
        self.assertTrue(os.path.exists(outdir))
        self.assertTrue(os.path.exists(adj.log_file))
        self.assertTrue(os.path.exists(adj.final_vcf))
        self.assertTrue(os.path.exists(adj.gramtools_build_dir))
        self.assertTrue(os.path.exists(adj.gramtools_quasimap_dir))
        self.assertTrue(os.path.exists(adj.clustered_vcf))
        self.assertTrue(os.path.exists(adj.plots_prefix + '.data.tsv'))
        self.assertTrue(os.path.exists(adj.plots_prefix + '.dp_hist.pdf'))
        self.assertTrue(os.path.exists(adj.plots_prefix + '.gt_conf_dp_scatter.pdf'))
        self.assertTrue(os.path.exists(adj.plots_prefix + '.gt_conf_hist.pdf'))

        # Now we've run the adjudicator, we have a gramtools
        # build directory. Rerun, but this time use the build
        # directory, so we test the gramtools_build_dir option
        outdir2 = 'tmp.adjudicator.out.2'
        gramtools_build_dir = adj.gramtools_build_dir
        if os.path.exists(outdir2):
            shutil.rmtree(outdir2)
        ref_fasta = os.path.join(data_dir, 'run.ref.fa')
        reads_file = os.path.join(data_dir, 'run.bwa.bam')
        # When gramtools build dir supplied, the Adjudicator assumes
        # one clsutered VCF file that matches the gramtools build run.
        # This is the clustered VCF made by the Adjudicator, so we
        # use that instead of the list of original VCF files
        vcf_files = [adj.clustered_vcf]
        adj = adjudicator.Adjudicator(outdir2, ref_fasta, [reads_file], vcf_files, gramtools_build_dir=gramtools_build_dir, clean=False, gramtools_kmer_size=5)
        adj.run()
        self.assertTrue(os.path.exists(outdir2))
        self.assertTrue(os.path.exists(adj.log_file))
        self.assertTrue(os.path.exists(adj.final_vcf))
        self.assertTrue(os.path.exists(adj.gramtools_build_dir))
        self.assertTrue(os.path.exists(adj.gramtools_quasimap_dir))
        self.assertTrue(os.path.exists(adj.clustered_vcf))
        self.assertTrue(os.path.exists(adj.plots_prefix + '.data.tsv'))
        self.assertTrue(os.path.exists(adj.plots_prefix + '.dp_hist.pdf'))
        self.assertTrue(os.path.exists(adj.plots_prefix + '.gt_conf_dp_scatter.pdf'))
        self.assertTrue(os.path.exists(adj.plots_prefix + '.gt_conf_hist.pdf'))
        self.assertFalse(os.path.exists(os.path.join(outdir2, 'gramtools.build')))
        shutil.rmtree(outdir)
        shutil.rmtree(outdir2)


    def test_run_empty_vcf_input_files(self):
        '''test run when input files have no variants'''
        outdir = 'tmp.adjudicator.out'
        if os.path.exists(outdir):
            shutil.rmtree(outdir)

        ref_fasta = os.path.join(data_dir, 'run.ref.fa')
        reads_file = os.path.join(data_dir, 'run.bwa.bam')
        vcf_files =  [os.path.join(data_dir, x) for x in ['run.calls.empty.1.vcf', 'run.calls.empty.2.vcf']]
        adj = adjudicator.Adjudicator(outdir, ref_fasta, [reads_file], vcf_files, clean=False, gramtools_kmer_size=5)
        with self.assertRaises(adjudicator.Error):
            adj.run()
        self.assertTrue(os.path.exists(outdir))
        self.assertTrue(os.path.exists(adj.log_file))
        self.assertFalse(os.path.exists(adj.final_vcf))
        self.assertFalse(os.path.exists(adj.gramtools_build_dir))
        self.assertFalse(os.path.exists(adj.gramtools_quasimap_dir))
        self.assertFalse(os.path.exists(adj.plots_prefix + '.data.tsv'))
        self.assertFalse(os.path.exists(adj.plots_prefix + '.dp_hist.pdf'))
        self.assertFalse(os.path.exists(adj.plots_prefix + '.gt_conf_dp_scatter.pdf'))
        self.assertFalse(os.path.exists(adj.plots_prefix + '.gt_conf_hist.pdf'))
        self.assertTrue(os.path.exists(adj.clustered_vcf))
        shutil.rmtree(outdir)
