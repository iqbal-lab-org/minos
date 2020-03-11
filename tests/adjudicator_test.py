import filecmp
import shutil
import os
import unittest

from minos import adjudicator

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "adjudicator")


class TestAdjudicator(unittest.TestCase):
    def test_get_gramtools_kmer_size(self):
        """test _get_gramtools_kmer_size"""
        build_dir = os.path.join(data_dir, "get_gramtools_kmer_size.build")
        self.assertEqual(
            42, adjudicator.Adjudicator._get_gramtools_kmer_size(build_dir, None)
        )
        self.assertEqual(
            42, adjudicator.Adjudicator._get_gramtools_kmer_size(build_dir, 20)
        )
        self.assertEqual(20, adjudicator.Adjudicator._get_gramtools_kmer_size(None, 20))
        self.assertEqual(
            10, adjudicator.Adjudicator._get_gramtools_kmer_size(None, None)
        )

    def test_run_clean_is_false(self):
        """test run when not cleaning up files afterwards"""
        # We're just testing that it doesn't crash.
        # Check the output files exist, but not their contents.
        # First run using splitting of VCF file.
        # Then run without splitting.
        outdir = "tmp.adjudicator.noclean.out"
        if os.path.exists(outdir):
            shutil.rmtree(outdir)
        ref_fasta = os.path.join(data_dir, "run.ref.fa")
        reads_file = os.path.join(data_dir, "run.bwa.bam")
        vcf_files = [
            os.path.join(data_dir, x) for x in ["run.calls.1.vcf", "run.calls.2.vcf"]
        ]
        adj = adjudicator.Adjudicator(
            outdir,
            ref_fasta,
            [reads_file],
            vcf_files,
            variants_per_split=3,
            clean=False,
            gramtools_kmer_size=5,
            genotype_simulation_iterations=1000,
        )
        adj.run()
        self.assertTrue(os.path.exists(outdir))
        self.assertTrue(os.path.exists(adj.log_file))
        self.assertTrue(os.path.exists(adj.final_vcf))
        self.assertTrue(os.path.exists(adj.clustered_vcf))

        # Clean up and then run without splitting
        shutil.rmtree(outdir)
        adj = adjudicator.Adjudicator(
            outdir,
            ref_fasta,
            [reads_file],
            vcf_files,
            clean=False,
            gramtools_kmer_size=5,
            genotype_simulation_iterations=1000,
        )
        adj.run()
        self.assertTrue(os.path.exists(outdir))
        self.assertTrue(os.path.exists(adj.log_file))
        self.assertTrue(os.path.exists(adj.final_vcf))
        self.assertTrue(os.path.exists(adj.gramtools_build_dir))
        self.assertTrue(os.path.exists(adj.gramtools_quasimap_dir))
        self.assertTrue(os.path.exists(adj.clustered_vcf))
        #self.assertTrue(os.path.exists(adj.plots_prefix + ".data.tsv"))
        #self.assertTrue(os.path.exists(adj.plots_prefix + ".dp_hist.pdf"))
        #self.assertTrue(os.path.exists(adj.plots_prefix + ".gt_conf_dp_scatter.pdf"))
        #self.assertTrue(os.path.exists(adj.plots_prefix + ".gt_conf_hist.pdf"))

        # Now we've run the adjudicator, we have a gramtools
        # build directory. Rerun, but this time use the build
        # directory, so we test the gramtools_build_dir option
        outdir2 = "tmp.adjudicator.out.2"
        gramtools_build_dir = adj.gramtools_build_dir
        if os.path.exists(outdir2):
            shutil.rmtree(outdir2)
        ref_fasta = os.path.join(data_dir, "run.ref.fa")
        reads_file = os.path.join(data_dir, "run.bwa.bam")
        # When gramtools build dir supplied, the Adjudicator assumes
        # one clsutered VCF file that matches the gramtools build run.
        # This is the clustered VCF made by the Adjudicator, so we
        # use that instead of the list of original VCF files
        vcf_files = [adj.clustered_vcf]
        adj = adjudicator.Adjudicator(
            outdir2,
            ref_fasta,
            [reads_file],
            vcf_files,
            gramtools_build_dir=gramtools_build_dir,
            clean=False,
            gramtools_kmer_size=5,
            genotype_simulation_iterations=1000,
        )
        adj.run()
        self.assertTrue(os.path.exists(outdir2))
        self.assertTrue(os.path.exists(adj.log_file))
        self.assertTrue(os.path.exists(adj.final_vcf))
        self.assertTrue(os.path.exists(adj.gramtools_build_dir))
        self.assertTrue(os.path.exists(adj.gramtools_quasimap_dir))
        self.assertTrue(os.path.exists(adj.clustered_vcf))
        #self.assertTrue(os.path.exists(adj.plots_prefix + ".data.tsv"))
        #self.assertTrue(os.path.exists(adj.plots_prefix + ".dp_hist.pdf"))
        #self.assertTrue(os.path.exists(adj.plots_prefix + ".gt_conf_dp_scatter.pdf"))
        #self.assertTrue(os.path.exists(adj.plots_prefix + ".gt_conf_hist.pdf"))
        self.assertFalse(os.path.exists(os.path.join(outdir2, "gramtools.build")))
        shutil.rmtree(outdir)
        shutil.rmtree(outdir2)

    def test_run_clean_is_true(self):
        """test run when we do clean files afterwards"""
        # We're just testing that it doesn't crash.
        # Check the output files exist, but not their contents.
        # First run using splitting of VCF file.
        # Then run without splitting.
        outdir = "tmp.adjudicator.clean.out"
        if os.path.exists(outdir):
            shutil.rmtree(outdir)
        ref_fasta = os.path.join(data_dir, "run.ref.fa")
        reads_file = os.path.join(data_dir, "run.bwa.bam")
        vcf_files = [
            os.path.join(data_dir, x) for x in ["run.calls.1.vcf", "run.calls.2.vcf"]
        ]
        adj = adjudicator.Adjudicator(
            outdir,
            ref_fasta,
            [reads_file],
            vcf_files,
            clean=True,
            gramtools_kmer_size=5,
            genotype_simulation_iterations=1000,
        )
        adj.run()
        self.assertTrue(os.path.exists(outdir))
        self.assertTrue(os.path.exists(os.path.join(outdir, "final.vcf")))
        shutil.rmtree(outdir)

    def test_run_empty_vcf_input_files(self):
        """test run when input files have no variants"""
        outdir = "tmp.adjudicator.out"
        if os.path.exists(outdir):
            shutil.rmtree(outdir)

        ref_fasta = os.path.join(data_dir, "run.ref.fa")
        reads_file = os.path.join(data_dir, "run.bwa.bam")
        vcf_files = [
            os.path.join(data_dir, x)
            for x in ["run.calls.empty.1.vcf", "run.calls.empty.2.vcf"]
        ]
        adj = adjudicator.Adjudicator(
            outdir,
            ref_fasta,
            [reads_file],
            vcf_files,
            clean=False,
            gramtools_kmer_size=5,
            genotype_simulation_iterations=1000,
        )
        with self.assertRaises(Exception):
            adj.run()
        self.assertTrue(os.path.exists(outdir))
        self.assertTrue(os.path.exists(adj.log_file))
        self.assertFalse(os.path.exists(adj.final_vcf))
        self.assertFalse(os.path.exists(adj.gramtools_build_dir))
        self.assertFalse(os.path.exists(adj.gramtools_quasimap_dir))
        #self.assertFalse(os.path.exists(adj.plots_prefix + ".data.tsv"))
        #self.assertFalse(os.path.exists(adj.plots_prefix + ".dp_hist.pdf"))
        #self.assertFalse(os.path.exists(adj.plots_prefix + ".gt_conf_dp_scatter.pdf"))
        #self.assertFalse(os.path.exists(adj.plots_prefix + ".gt_conf_hist.pdf"))
        self.assertTrue(os.path.exists(adj.clustered_vcf))
        shutil.rmtree(outdir)

    def test_add_gt_conf_percentile_and_filters_to_vcf_file(self):
        """test _add_gt_conf_percentile_and_filters_to_vcf_file"""
        original_file = os.path.join(
            data_dir, "add_gt_conf_percentile_to_vcf_file.in.vcf"
        )
        tmp_file = "tmp.adjudicator.add_gt_conf_percentile_to_vcf_file.vcf"
        expect_file = os.path.join(
            data_dir, "add_gt_conf_percentile_to_vcf_file.expect.vcf"
        )
        shutil.copyfile(original_file, tmp_file)
        error_rate = 0.00026045894282438386
        adjudicator.Adjudicator._add_gt_conf_percentile_and_filters_to_vcf_file(
            tmp_file, 60, 100, error_rate, iterations=1000, min_dp=2, min_gcp=2.5, min_frs=0.9
        )
        self.assertTrue(filecmp.cmp(tmp_file, expect_file, shallow=False))
        os.unlink(tmp_file)

    def test_0MeanDepth_stillRuns(self):
        """
        When mean depth is 0, we can get math errors: math.log(0) in genotype likelihood computation,
        and division by 0 in genotype confidence simulation.

        Note the former only actually occurs if there is a variant site with non-zero coverage;
        in this case, mean depth can get set to 0 due to rounding imprecision. This is tested in genotyper unit tests.
        """

        outdir = "tmp.adjudicator.out"
        if os.path.exists(outdir):
            shutil.rmtree(outdir)
        ref_fasta = os.path.join(data_dir, "run.ref.fa")
        reads_file = os.path.join(data_dir, "no_map_reads.fastq")
        vcf_files = [os.path.join(data_dir, "run.calls.1.vcf")]

        adj = adjudicator.Adjudicator(
            outdir,
            ref_fasta,
            [reads_file],
            vcf_files,
            clean=False,
            gramtools_kmer_size=5,
        )
        adj.run()
        # Make sure the coverage is 0
        self.assertEqual(adj.mean_depths[0], 0)
        # And also the test passes if it raises no math related errors.
