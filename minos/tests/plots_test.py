import shutil
import os
import unittest

import pandas as pd

from minos import plots

modules_dir = os.path.dirname(os.path.abspath(plots.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'plots')

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
        plots.scatter_plot_gt_conf_vs_dp_colour_by_tp_fp(data, tmpfile)
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

