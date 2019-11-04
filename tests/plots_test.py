import filecmp
import shutil
import os
import unittest

import pandas as pd

from minos import plots

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, 'data', 'plots')

class TestPlots(unittest.TestCase):
    def test_load_dp_and_gt_conf_data_from_file(self):
        expect_dict = {
            'DP': [3, 6, 10, 42],
            'GT_CONF': [5.0, 10.0, 9.0, 42],
        }
        expect_df = pd.DataFrame(data=expect_dict)
        infile = os.path.join(data_dir, 'load_dp_and_gt_conf_data_from_file.tsv')
        got = plots.load_dp_and_gt_conf_data_from_file(infile)
        pd.util.testing.assert_frame_equal(expect_df, got)

        expect_dict['TP_OR_FP'] = ['FP', 'TP', 'TP', 'TP']
        expect_df = pd.DataFrame(data=expect_dict)
        infile = os.path.join(data_dir, 'load_dp_and_gt_conf_data_from_file.with_tp_fp.tsv')
        got = plots.load_dp_and_gt_conf_data_from_file(infile)
        pd.util.testing.assert_frame_equal(expect_df, got)


    def test_scatter_plot_gt_conf_vs_dp(self):
        infile = os.path.join(data_dir, 'dp_and_gt_conf_data.tsv')
        data = plots.load_dp_and_gt_conf_data_from_file(infile)
        tmpfile = 'tmp.test.dp_and_gt_conf_data.tsv.pdf'
        plots.scatter_plot_gt_conf_vs_dp(data, tmpfile)
        self.assertTrue(os.path.exists(tmpfile))
        self.assertNotEqual(0, os.stat(tmpfile).st_size)
        os.unlink(tmpfile)


    def test_scatter_plot_gt_conf_vs_dp_colour_by_tp_fp(self):
        infile = os.path.join(data_dir, 'dp_and_gt_conf_data.with_tp_fp.tsv')
        data = plots.load_dp_and_gt_conf_data_from_file(infile)
        tmpfile = 'tmp.test.scatter_plot_gt_conf_vs_dp_colour_by_tp_fp.pdf'
        plots.scatter_plot_gt_conf_vs_dp_colour_by_tp_fp(data, tmpfile, {'TP', 'FP'})
        self.assertTrue(os.path.exists(tmpfile))
        self.assertNotEqual(0, os.stat(tmpfile).st_size)
        os.unlink(tmpfile)


    def test_histogram_of_one_dataframe_column(self):
        infile = os.path.join(data_dir, 'dp_and_gt_conf_data.tsv')
        data = plots.load_dp_and_gt_conf_data_from_file(infile)
        tmpfile = 'tmp.test.histogram_of_one_dataframe_column.tsv.pdf'
        plots.histogram_of_one_dataframe_column(data, 'DP', tmpfile)
        self.assertTrue(os.path.exists(tmpfile))
        self.assertNotEqual(0, os.stat(tmpfile).st_size)
        os.unlink(tmpfile)


    def test_histogram_of_one_dataframe_column_color_by_tp_fp(self):
        infile = os.path.join(data_dir, 'dp_and_gt_conf_data.with_tp_fp.tsv')
        data = plots.load_dp_and_gt_conf_data_from_file(infile)
        tmpfile = 'tmp.test.histogram_of_one_dataframe_column.tsv.pdf'
        plots.histogram_of_one_dataframe_column_color_by_tp_fp(data, 'DP', tmpfile)
        self.assertTrue(os.path.exists(tmpfile))
        self.assertNotEqual(0, os.stat(tmpfile).st_size)
        os.unlink(tmpfile)


    def test_minos_vcf_to_plot_data(self):
        infile = os.path.join(data_dir, 'minos_vcf_to_plot_data.in.vcf')
        expect_file = os.path.join(data_dir, 'minos_vcf_to_plot_data.expect.tsv')
        tmpfile = 'tmp.test.minos_vcf_to_plot_data.tsv'
        self.assertEqual({'TP', 'FP'}, plots.minos_vcf_to_plot_data(infile, tmpfile))
        self.assertTrue(filecmp.cmp(expect_file, tmpfile, shallow=False))
        os.unlink(tmpfile)

        infile = os.path.join(data_dir, 'minos_vcf_to_plot_data.partial_check_geno.in.vcf')
        expect_file = os.path.join(data_dir, 'minos_vcf_to_plot_data.partial_check_geno.expect.tsv')
        self.assertEqual({'TP'}, plots.minos_vcf_to_plot_data(infile, tmpfile))
        self.assertTrue(filecmp.cmp(expect_file, tmpfile, shallow=False))
        os.unlink(tmpfile)

        infile = os.path.join(data_dir, 'minos_vcf_to_plot_data.no_check_geno.in.vcf')
        expect_file = os.path.join(data_dir, 'minos_vcf_to_plot_data.no_check_geno.expect.tsv')
        self.assertEqual(set(), plots.minos_vcf_to_plot_data(infile, tmpfile))
        self.assertTrue(filecmp.cmp(expect_file, tmpfile, shallow=False))
        os.unlink(tmpfile)

        infile = os.path.join(data_dir, 'minos_vcf_to_plot_data.partial_dp_and_gt_conf.in.vcf')
        expect_file = os.path.join(data_dir, 'minos_vcf_to_plot_data.partial_dp_and_gt_conf.expect.tsv')
        self.assertEqual(set(), plots.minos_vcf_to_plot_data(infile, tmpfile))
        self.assertTrue(filecmp.cmp(expect_file, tmpfile, shallow=False))
        os.unlink(tmpfile)

        infile = os.path.join(data_dir, 'minos_vcf_to_plot_data.no_dp_and_gt_conf.in.vcf')
        self.assertIsNone(plots.minos_vcf_to_plot_data(infile, tmpfile))
        self.assertFalse(os.path.exists(tmpfile))


    def test_plots_from_minos_vcf(self):
        outprefix = 'tmp.plots_from_minos_vcf.out'
        expect_files = [outprefix + x + '.pdf' for x in ['.gt_conf_dp_scatter', '.dp_hist', '.gt_conf_hist']]
        expect_files.append(outprefix + '.data.tsv')

        infile = os.path.join(data_dir, 'minos_vcf_to_plot_data.in.vcf')
        plots.plots_from_minos_vcf(infile, outprefix)
        for f in expect_files:
            self.assertTrue(os.path.exists(f))
            self.assertNotEqual(0, os.stat(f).st_size)
            os.unlink(f)

        infile = os.path.join(data_dir, 'minos_vcf_to_plot_data.no_check_geno.in.vcf')
        plots.plots_from_minos_vcf(infile, outprefix)
        for f in expect_files:
            self.assertTrue(os.path.exists(f))
            self.assertNotEqual(0, os.stat(f).st_size)
            os.unlink(f)

        infile = os.path.join(data_dir, 'minos_vcf_to_plot_data.fp_only.vcf')
        plots.plots_from_minos_vcf(infile, outprefix)
        for f in expect_files:
            self.assertTrue(os.path.exists(f))
            self.assertNotEqual(0, os.stat(f).st_size)
            os.unlink(f)

        infile = os.path.join(data_dir, 'minos_vcf_to_plot_data.tp_only.vcf')
        plots.plots_from_minos_vcf(infile, outprefix)
        for f in expect_files:
            self.assertTrue(os.path.exists(f))
            self.assertNotEqual(0, os.stat(f).st_size)
            os.unlink(f)

        infile = os.path.join(data_dir, 'minos_vcf_to_plot_data.no_dp_and_gt_conf.in.vcf')
        plots.plots_from_minos_vcf(infile, outprefix)
        for f in expect_files:
            self.assertFalse(os.path.exists(f))


