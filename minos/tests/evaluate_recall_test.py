import filecmp
import os
import unittest
import glob

import pyfastaq

from cluster_vcf_records import vcf_file_read

from minos import evaluate_recall

modules_dir = os.path.dirname(os.path.abspath(evaluate_recall.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'evaluate_recall')


class TestEvaluateRecall(unittest.TestCase):
    def test_write_vars_plus_flanks_to_fasta_alt(self):
        vcfref_file_in = os.path.join(data_dir, 'vcfref.fa')
        sample_file_in = os.path.join(data_dir, 'sample1a.vcf')
        tmp_out = 'tmp.write_vars_plus_flanks_to_fasta.alt.fa'
        expected_out = os.path.join(data_dir, 'sample1a.plusflanks.fa')

        vcf_header, vcf_records = vcf_file_read.vcf_file_to_dict(sample_file_in, sort=True,
                                                                 remove_useless_start_nucleotides=True)
        vcf_ref_seqs = {}
        pyfastaq.tasks.file_to_dict(vcfref_file_in, vcf_ref_seqs)
        flank = 5
        evaluate_recall.EvaluateRecall._write_vars_plus_flanks_to_fasta(tmp_out, vcf_records, vcf_ref_seqs, flank, alt_only=True)

        self.assertTrue(filecmp.cmp(expected_out, tmp_out, shallow=False))
        os.unlink(tmp_out)

    def test_parse_sam_file_and_vcf(self):
        samfile = os.path.join(data_dir, 'mapped.sam')
        vcffile = os.path.join(data_dir, 'sample1a.vcf')
        flank = 5
        allow_mismatches = False
        evaluate_recall.EvaluateRecall._index_vcf(vcffile)
        found, gt_conf, allele = evaluate_recall.EvaluateRecall._parse_sam_file_and_vcf(samfile, vcffile + ".gz",
                                                                                        flank,
                                                                                        allow_mismatches)

        exp_found = ['1', '1', '0', '0', '1']  # nb doesn't currently handle '.' alleles
        exp_gt_conf = [42, 42, 0, 0, 0]
        exp_allele = ['0', '1', '0', '0', '0']
        self.assertEqual(exp_found, found)
        self.assertEqual(exp_gt_conf, gt_conf)
        self.assertEqual(exp_allele, allele)

    def test_run_with_filter_cluster_include_ref_alleles(self):
        '''test run with filtering and clustering'''
        truth_vcf_in = os.path.join(data_dir, 'truth.vcf')
        truth_ref_in = os.path.join(data_dir, 'truth_ref.fa')
        query_vcf_in = os.path.join(data_dir, 'sample1a.vcf')
        query_ref_in = os.path.join(data_dir, 'vcfref.fa')
        exp_out = os.path.join(data_dir, 'exp_stats.tsv')
        exp_gt_conf = os.path.join(data_dir, 'exp_gt_conf.tsv')

        tmp_out = 'tmp.evaluate_recall.out'
        verifier = evaluate_recall.EvaluateRecall(truth_vcf_in, truth_ref_in, query_vcf_in, query_ref_in, tmp_out,
                                                  flank_length=5, discard_ref_calls=False)
        verifier.run()

        self.assertTrue(filecmp.cmp(exp_out, tmp_out + '.stats.tsv', shallow=False))
        self.assertTrue(filecmp.cmp(exp_gt_conf, tmp_out + '.gt_conf_hist.tsv', shallow=False))
        for f in glob.glob(tmp_out + '*'):
            os.unlink(f)
