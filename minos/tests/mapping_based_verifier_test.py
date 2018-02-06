import filecmp
import os
import unittest

from cluster_vcf_records import vcf_record
import pyfastaq

from minos import mapping_based_verifier

modules_dir = os.path.dirname(os.path.abspath(mapping_based_verifier.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'mapping_based_verifier')


class TestMappingBasedVerifier(unittest.TestCase):
    def test_write_vars_plus_flanks_to_fasta(self):
        '''test _write_vars_plus_flanks_to_fasta'''
        ref_seqs = {
            #                                         12345678901234567890
            'ref1': pyfastaq.sequences.Fasta('ref1', 'AGTGGATAGCTAGCTAGAGA'),
            'ref2': pyfastaq.sequences.Fasta('ref2', 'AGGAGAGAGAGAGAGAA'),
            'ref3': pyfastaq.sequences.Fasta('ref3', 'AGCTTCATAGAGAGGTTTA'),
        }

        vcf_records = {
            'ref1': [
                vcf_record.VcfRecord('ref1\t3\tid_1\tT\tC,AG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80'),
                vcf_record.VcfRecord('ref1\t10\tid_2\tCT\tA\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80'),
            ],
            'ref3': [
                vcf_record.VcfRecord('ref3\t4\tid_3\tT\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80'),
            ],
        }

        tmp_file = 'tmp.mapping_based_verifier.write_vars_plus_flanks_to_fasta.fa'
        mapping_based_verifier.MappingBasedVerifier._write_vars_plus_flanks_to_fasta(tmp_file, vcf_records, ref_seqs, 3)

        expected_file = os.path.join(data_dir, 'write_vars_plus_flanks_to_fasta.fa')
        self.assertTrue(filecmp.cmp(expected_file, tmp_file, shallow=False))
        os.unlink(tmp_file)


    def test_map_seqs_to_ref(self):
        '''test _map_seqs_to_ref'''
        ref = os.path.join(data_dir, 'map_seqs_to_ref.ref.fa')
        reads = os.path.join(data_dir, 'map_seqs_to_ref.reads.fa')
        tmp_out = 'tmp.map_seqs_to_ref.sam'
        if os.path.exists(tmp_out):
            os.unlink(tmp_out)
        mapping_based_verifier.MappingBasedVerifier._map_seqs_to_ref(reads, ref, tmp_out)
        # trust the mapping. Just check output was written
        self.assertTrue(os.path.exists(tmp_out))
        os.unlink(tmp_out)


    def test_check_called_genotype(self):
        '''test _check_called_genotype'''
        vcf = vcf_record.VcfRecord('ref\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT:MINOS_CHECK_PASS\t0/0:0,0,0')
        expect_vcf = vcf_record.VcfRecord('ref\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT:MINOS_CHECK_PASS:MINOS_CHECK_GENOTYPE\t0/0:0,0,0:0')
        mapping_based_verifier.MappingBasedVerifier._check_called_genotype(vcf)
        self.assertEqual(expect_vcf, vcf)

        vcf = vcf_record.VcfRecord('ref\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT:MINOS_CHECK_PASS\t0/0:0,0,1')
        expect_vcf = vcf_record.VcfRecord('ref\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT:MINOS_CHECK_PASS:MINOS_CHECK_GENOTYPE\t0/0:0,0,1:0')
        mapping_based_verifier.MappingBasedVerifier._check_called_genotype(vcf)
        self.assertEqual(expect_vcf, vcf)

        vcf = vcf_record.VcfRecord('ref\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT:MINOS_CHECK_PASS\t0/0:1,0,0')
        expect_vcf = vcf_record.VcfRecord('ref\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT:MINOS_CHECK_PASS:MINOS_CHECK_GENOTYPE\t0/0:1,0,0:1')
        mapping_based_verifier.MappingBasedVerifier._check_called_genotype(vcf)
        self.assertEqual(expect_vcf, vcf)

        vcf = vcf_record.VcfRecord('ref\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT:MINOS_CHECK_PASS\t0/1:1,0,0')
        expect_vcf = vcf_record.VcfRecord('ref\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT:MINOS_CHECK_PASS:MINOS_CHECK_GENOTYPE\t0/1:1,0,0:HET')
        mapping_based_verifier.MappingBasedVerifier._check_called_genotype(vcf)
        self.assertEqual(expect_vcf, vcf)

        vcf = vcf_record.VcfRecord('ref\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT:MINOS_CHECK_PASS\t1/1:0,1,0')
        expect_vcf = vcf_record.VcfRecord('ref\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT:MINOS_CHECK_PASS:MINOS_CHECK_GENOTYPE\t1/1:0,1,0:1')
        mapping_based_verifier.MappingBasedVerifier._check_called_genotype(vcf)
        self.assertEqual(expect_vcf, vcf)

        vcf = vcf_record.VcfRecord('ref\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tMINOS_CHECK_PASS\t0,0,0')
        expect_vcf = vcf_record.VcfRecord('ref\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tMINOS_CHECK_PASS:MINOS_CHECK_GENOTYPE\t0,0,0:UNKNOWN_NO_GT')
        mapping_based_verifier.MappingBasedVerifier._check_called_genotype(vcf)
        self.assertEqual(expect_vcf, vcf)


    def test_get_missing_vcf_records(self):
        '''test _get_missing_vcf_records'''
        vcfs_should_find = {
            'ref.1': [
                vcf_record.VcfRecord('ref.1\t42\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT\t1/1'),
                vcf_record.VcfRecord('ref.1\t100\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT\t2/2'),
                vcf_record.VcfRecord('ref.1\t200\t.\tC\tA,G\t42.0\t.\tDP4=42\tGT\t2/2'),
            ],
            'ref.2': [
                vcf_record.VcfRecord('ref.2\t50\t.\tT\tA\t42.0\t.\tDP4=42\tGT\t1/1'),
            ],
        }

        vcfs_to_check = {
            'ref.1': [
                vcf_record.VcfRecord('ref.1\t12\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT\t1/1'),
                vcf_record.VcfRecord('ref.1\t42\t.\tT\tG,A\t42.0\t.\tDP4=42\tGT\t2/2'),
                vcf_record.VcfRecord('ref.1\t200\t.\tC\tA,G\t42.0\t.\tDP4=42\tGT\t1/1'),
            ],
            'ref.3': [
                vcf_record.VcfRecord('ref.3\t25\t.\tA\tT\t42.0\t.\tDP4=42\tGT\t1/1'),
            ],
        }

        expected = {
            'ref.1': [
                vcf_record.VcfRecord('ref.1\t100\t.\tT\tA,G\t42.0\t.\tDP4=42\tGT\t2/2'),
                vcf_record.VcfRecord('ref.1\t200\t.\tC\tA,G\t42.0\t.\tDP4=42\tGT\t2/2'),
            ],
            'ref.2': [
                vcf_record.VcfRecord('ref.2\t50\t.\tT\tA\t42.0\t.\tDP4=42\tGT\t1/1'),
            ],
        }
        got = mapping_based_verifier.MappingBasedVerifier._get_missing_vcf_records(vcfs_to_check, vcfs_should_find)
        self.assertEqual(expected, got)


    def test_get_total_length_of_expected_regions_called(self):
        '''test _get_total_length_of_expected_regions_called'''
        expected_regions = {
            'ref.1': [
                pyfastaq.intervals.Interval(101, 200), # 100 long, 92 get called
                pyfastaq.intervals.Interval(251, 260), # 10 long, none get called
            ],
            'ref.2': [
                pyfastaq.intervals.Interval(42,43), # 2 long, none get called
            ],
        }

        called_vcf_records = {
            'ref.1': [
                vcf_record.VcfRecord('ref.1\t100\t.\tACGTACTGTA\tA,G\t42.0\t.\tDP4=42\tGT\t2/2'),
            ],
        }

        got_all, got_called  = mapping_based_verifier.MappingBasedVerifier._get_total_length_of_expected_regions_called(expected_regions, called_vcf_records)
        self.assertEqual(112, got_all)
        self.assertEqual(8, got_called)


    def test_run(self):
        '''test run'''
        vcf_file_in = os.path.join(data_dir, 'run.calls.vcf')
        vcf_reference_file = os.path.join(data_dir, 'run.ref.fa')
        verify_reference_file = os.path.join(data_dir, 'run.ref.mutated.fa')
        tmp_out = 'tmp.mapping_based_verifier.out'
        verifier = mapping_based_verifier.MappingBasedVerifier(vcf_file_in, vcf_reference_file, verify_reference_file, tmp_out, flank_length=31)
        verifier.run()
        expected_out = os.path.join(data_dir, 'run.out')
        for suffix in ['.stats.tsv', '.vcf']:
            expected_file = expected_out + suffix
            got_file = tmp_out + suffix
            self.assertTrue(filecmp.cmp(expected_file, got_file, shallow=False))
            os.unlink(got_file)
        samfile = tmp_out + '.sam'
        self.assertTrue(os.path.exists(samfile))
        os.unlink(samfile)


class TestOther(unittest.TestCase):
    def test_sam_reader(self):
        '''test sam_reader'''
        samfile = os.path.join(data_dir, 'sam_reader.sam')
        expected = [
            ['read1.0.0.0', 'read1.0.0.1'],
            ['read1.10.1.0'],
            ['read2.100.0.0'],
            ['read3.42.0.0', 'read3.42.0.1', 'read3.42.0.2'],
        ]

        sread = mapping_based_verifier.sam_reader(samfile)
        i = 0
        for got_list in sread:
            got_names = [x.query_name for x in got_list]
            self.assertEqual(expected[i], got_names)
            i += 1


