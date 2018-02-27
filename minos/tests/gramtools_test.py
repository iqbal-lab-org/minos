import datetime
import shutil
import os
import unittest

from cluster_vcf_records import vcf_file_read, vcf_record

from minos import gramtools
from minos import __version__ as minos_version

modules_dir = os.path.dirname(os.path.abspath(gramtools.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'gramtools')


class TestGramtools(unittest.TestCase):
    def test_run_gramtools_build(self):
        '''test run_gramtools_build'''
        tmp_out_build = 'tmp.run_gramtools.out.build'
        if os.path.exists(tmp_out_build):
            shutil.rmtree(tmp_out_build)
        vcf_file = os.path.join(data_dir, 'run_gramtools.calls.vcf')
        ref_file = os.path.join(data_dir, 'run_gramtools.ref.fa')
        gramtools.run_gramtools_build(tmp_out_build, vcf_file, ref_file, 150)
        self.assertTrue(os.path.exists(tmp_out_build))
        shutil.rmtree(tmp_out_build)


    def test_run_gramtools(self):
        '''test run_gramtools'''
        tmp_out_build = 'tmp.run_gramtools.out.build'
        tmp_out_quasimap = 'tmp.run_gramtools.out.quasimap'
        if os.path.exists(tmp_out_build):
            shutil.rmtree(tmp_out_build)
        if os.path.exists(tmp_out_quasimap):
            shutil.rmtree(tmp_quasimap)
        vcf_file = os.path.join(data_dir, 'run_gramtools.calls.vcf')
        ref_file = os.path.join(data_dir, 'run_gramtools.ref.fa')
        reads_file = os.path.join(data_dir, 'run_gramtools.reads.fq')
        gramtools.run_gramtools(tmp_out_build, tmp_out_quasimap, vcf_file, ref_file, reads_file, 150)
        # We're trusing gramtools output for this test. The point here is to check
        # that gramtools can run. Parsing its output is checked elsewhere.
        self.assertTrue(os.path.exists(tmp_out_build))
        self.assertTrue(os.path.exists(tmp_out_quasimap))
        self.assertTrue(os.path.exists(os.path.join(tmp_out_quasimap, 'allele_base_coverage.json')))
        self.assertTrue(os.path.exists(os.path.join(tmp_out_quasimap, 'grouped_allele_counts_coverage.json')))
        shutil.rmtree(tmp_out_build)
        shutil.rmtree(tmp_out_quasimap)


    def test_run_gramtools_fails(self):
        '''test run_gramtools when fails'''
        # Don't trust error code. Instead, we check
        # that gramtools wrote the quesimap files we expected, as this
        # is a good proxy for success. One way to stop these files
        # from being written is to have no variants in the input VCF,
        # so that's what we do here
        tmp_out_build = 'tmp.run_gramtools.fail.out.build'
        tmp_out_quasimap = 'tmp.run_gramtools.fail.out.quasimap'
        if os.path.exists(tmp_out_build):
            shutil.rmtree(tmp_out_build)
        if os.path.exists(tmp_out_quasimap):
            shutil.rmtree(tmp_quasimap)
        vcf_file = os.path.join(data_dir, 'run_gramtools.empty.vcf')
        ref_file = os.path.join(data_dir, 'run_gramtools.ref.fa')
        reads_file = os.path.join(data_dir, 'run_gramtools.reads.fq')
        with self.assertRaises(gramtools.Error):
            gramtools.run_gramtools(tmp_out_build, tmp_out_quasimap, vcf_file, ref_file, reads_file, 150)
        shutil.rmtree(tmp_out_build)
        shutil.rmtree(tmp_out_quasimap)


    def test_run_gramtools_two_reads_files(self):
        '''test run_gramtools'''
        tmp_out_build = 'tmp.run_gramtools.2files.out.build'
        tmp_out_quasimap = 'tmp.run_gramtools.2files.out.quasimap'
        if os.path.exists(tmp_out_build):
            shutil.rmtree(tmp_out_build)
        if os.path.exists(tmp_out_quasimap):
            shutil.rmtree(tmp_quasimap)
        vcf_file = os.path.join(data_dir, 'run_gramtools.calls.vcf')
        ref_file = os.path.join(data_dir, 'run_gramtools.ref.fa')
        reads_file1 = os.path.join(data_dir, 'run_gramtools.reads_1.fq')
        reads_file2 = os.path.join(data_dir, 'run_gramtools.reads_2.fq')
        gramtools.run_gramtools(tmp_out_build, tmp_out_quasimap, vcf_file, ref_file, [reads_file1, reads_file2], 150)
        # We're trusing gramtools output for this test. The point here is to check
        # that gramtools can run. Parsing its output is checked elsewhere.
        self.assertTrue(os.path.exists(tmp_out_build))
        self.assertTrue(os.path.exists(tmp_out_quasimap))
        self.assertTrue(os.path.exists(os.path.join(tmp_out_quasimap, 'allele_base_coverage.json')))
        self.assertTrue(os.path.exists(os.path.join(tmp_out_quasimap, 'grouped_allele_counts_coverage.json')))
        shutil.rmtree(tmp_out_build)
        shutil.rmtree(tmp_out_quasimap)


    def test_load_gramtools_vcf_and_allele_coverage_files(self):
        '''test load_gramtools_vcf_and_allele_coverage_files'''
        vcf_file = os.path.join(data_dir, 'load_gramtools_vcf_and_allele_coverage.vcf')
        quasimap_dir = os.path.join(data_dir, 'load_gramtools_vcf_and_allele_coverage_files.quasimap')
        got_mean_depth, got_vcf_header, got_vcf_records, got_allele_coverage, got_allele_groups  = gramtools.load_gramtools_vcf_and_allele_coverage_files(vcf_file, quasimap_dir)

        expected_header, expected_vcf_records = vcf_file_read.vcf_file_to_list(vcf_file)
        self.assertEqual(expected_header, got_vcf_header)
        self.assertEqual(expected_vcf_records, got_vcf_records)
        self.assertEqual(10.500, got_mean_depth)

        # now test bad files cause error to be raised
        vcf_file = os.path.join(data_dir, 'load_gramtools_vcf_and_allele_coverage.short.vcf')
        with self.assertRaises(gramtools.Error):
            gramtools.load_gramtools_vcf_and_allele_coverage_files(vcf_file, quasimap_dir)

        vcf_file = os.path.join(data_dir, 'load_gramtools_vcf_and_allele_coverage.long.vcf')
        with self.assertRaises(gramtools.Error):
            gramtools.load_gramtools_vcf_and_allele_coverage_files(vcf_file, quasimap_dir)

        vcf_file = os.path.join(data_dir, 'load_gramtools_vcf_and_allele_coverage.bad_allele_count.vcf')
        with self.assertRaises(gramtools.Error):
            gramtools.load_gramtools_vcf_and_allele_coverage_files(vcf_file, quasimap_dir)


    def test_update_vcf_record_using_gramtools_allele_depths_heterozygous(self):
        '''test update_using_gramtools_allele_depths heterozygous'''
        record = vcf_record.VcfRecord('ref\t4\t.\tT\tG,TC\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        allele_combination_cov = {'1': 9, '2': 7, '3': 1}
        allele_groups_dict = {'1': {0}, '2': {1}, '3': {1,2}}
        allele_per_base_cov = [[9],[7],[1,0]]
        expected = vcf_record.VcfRecord('ref\t4\t.\tT\tG,TC\t.\t.\tDP=17\tGT:COV:GT_CONF\t0/1:9,7,0:47.69')
        mean_depth = 15
        error_rate = 0.001
        gramtools.update_vcf_record_using_gramtools_allele_depths(record, allele_combination_cov, allele_per_base_cov, allele_groups_dict, mean_depth, error_rate)
        self.assertEqual(expected, record)


    def test_update_vcf_record_using_gramtools_allele_depths_homozygous(self):
        '''test update_using_gramtools_allele_depths homozygous'''
        record = vcf_record.VcfRecord('ref\t4\t.\tT\tG,TC\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        allele_depths = {'T': 1, 'G': 80}
        allele_combination_cov = {'1': 1, '2': 80}
        allele_groups_dict = {'1': {0}, '2': {1}, '3': {1,2}}
        allele_per_base_cov = [[1],[80],[0,0]]
        expected = vcf_record.VcfRecord('ref\t4\t.\tT\tG,TC\t.\t.\tDP=81\tGT:COV:GT_CONF\t1/1:1,80,0:44.79')
        mean_depth = 85
        error_rate = 0.001
        gramtools.update_vcf_record_using_gramtools_allele_depths(record, allele_combination_cov, allele_per_base_cov, allele_groups_dict, mean_depth, error_rate)
        self.assertEqual(expected, record)


    def test_write_vcf_annotated_using_coverage_from_gramtools(self):
        '''test write_vcf_annotated_using_coverage_from_gramtools'''
        vcf_file_in = os.path.join(data_dir, 'write_vcf_annotated_using_coverage_from_gramtools.in.vcf')
        quasimap_dir = os.path.join(data_dir, 'write_vcf_annotated_using_coverage_from_gramtools.quasimap')
        mean_depth, vcf_header, vcf_records, allele_coverage, allele_groups  = gramtools.load_gramtools_vcf_and_allele_coverage_files(vcf_file_in, quasimap_dir)
        tmp_outfile = 'tmp.gramtools.write_vcf_annotated_using_coverage_from_gramtools.vcf'
        error_rate = 0.001
        gramtools.write_vcf_annotated_using_coverage_from_gramtools(mean_depth, vcf_records, allele_coverage, allele_groups, error_rate, tmp_outfile, sample_name='sample_42')
        expected_vcf = os.path.join(data_dir, 'write_vcf_annotated_using_coverage_from_gramtools.out.vcf')
        # Today's date and the verison of minos get added to the header.
        # We'll have to take account
        # of those by fixing what we get from the expected file
        expected_header, expected_vcf_records = vcf_file_read.vcf_file_to_list(expected_vcf)
        got_header, got_vcf_records = vcf_file_read.vcf_file_to_list(tmp_outfile)
        for i in range(len(expected_header)):
            if expected_header[i].startswith('##fileDate='):
                expected_header[i] = '##fileDate=' + str(datetime.date.today())
            elif expected_header[i].startswith('##source=minos'):
                expected_header[i] = '##source=minos, version ' + minos_version

        self.assertEqual(expected_header, got_header)
        self.assertEqual(expected_vcf_records, got_vcf_records)
        os.unlink(tmp_outfile)


    def test_load_allele_files(self):
        '''test load_allele_files'''
        expected_counts_list = [
            ({"0": 10, "1": 3, "2": 9}, [[1, 1, 1], [1, 1, 0]]),
            ({"0": 30, "3": 2, "4": 1}, [[0, 0, 1], [2, 2, 0], [2, 2, 0, 1, 3]]),
        ]

        expected_groups_dict = {
            "0": {0},
            "1": {1},
            "2": {0, 1},
            "3": {2},
            "4": {1, 2},
        }

        allele_base_counts_file = os.path.join(data_dir, 'load_allele_files.allele_base_counts.json')
        grouped_allele_counts_file = os.path.join(data_dir, 'load_allele_files.grouped_allele_counts.json')
        got_counts_list, got_groups_dict = gramtools.load_allele_files(allele_base_counts_file, grouped_allele_counts_file)
        self.assertEqual(expected_counts_list, got_counts_list)
        self.assertEqual(expected_groups_dict, got_groups_dict)

