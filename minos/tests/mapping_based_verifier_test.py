import filecmp
import os
import unittest

from cluster_vcf_records import vcf_record
import pysam
import pyfastaq

from minos import mapping_based_verifier

modules_dir = os.path.dirname(os.path.abspath(mapping_based_verifier.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'mapping_based_verifier')


class TestMappingBasedVerifier(unittest.TestCase):
    def test_needleman_wunsch(self):
        '''test _needleman_wunsch'''
        seq1 = 'ACGTGTCACAG'
        seq2 = 'AGTCTGACATG'
        aln1, aln2 = mapping_based_verifier.MappingBasedVerifier._needleman_wunsch(seq1, seq2)
        expect1 = 'ACGTGTCACA-G'
        expect2 = 'A-GTCTGACATG'
        self.assertEqual(expect1, aln1)
        self.assertEqual(expect2, aln2)


    def test_edit_distance_from_alignment_strings(self):
        '''test _edit_distance_from_alignment_strings'''
        self.assertEqual(0, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('A', 'A'))
        self.assertEqual(1, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('A', 'C'))
        self.assertEqual(1, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('A-G', 'ACG'))
        self.assertEqual(1, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('A--G', 'ACTG'))
        self.assertEqual(1, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('ACG', 'A-G'))
        self.assertEqual(1, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('ACTG', 'A--G'))
        self.assertEqual(2, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('ACTC', 'A--G'))
        self.assertEqual(2, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('A--G', 'ACTC'))
        self.assertEqual(1, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('A-TG', 'AC-G'))
        self.assertEqual(1, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('AC-G', 'A-TG'))
        self.assertEqual(1, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('A--G', 'A-CG'))
        self.assertEqual(1, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('A-CG', 'A--G'))
        self.assertEqual(2, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('A--TCAC-G', 'AGT--ACTG'))
        self.assertEqual(2, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('AGT--ACTG', 'A--TCAC-G'))
        self.assertEqual(3, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('A--TCAC-GG', 'AGT--ACTG-'))
        self.assertEqual(3, mapping_based_verifier.MappingBasedVerifier._edit_distance_from_alignment_strings('AGT--ACTG-', 'A--TCAC-GG'))


    def test_edit_distance_between_seqs(self):
        '''test _edit_distance_between_seqs'''
        self.assertEqual(0, mapping_based_verifier.MappingBasedVerifier._edit_distance_between_seqs('A', 'A'))
        self.assertEqual(1, mapping_based_verifier.MappingBasedVerifier._edit_distance_between_seqs('A', 'C'))
        self.assertEqual(1, mapping_based_verifier.MappingBasedVerifier._edit_distance_between_seqs('AG', 'ACG'))
        self.assertEqual(2, mapping_based_verifier.MappingBasedVerifier._edit_distance_between_seqs('AGTGCAT', 'ACGTGCGT'))


    def test_load_exclude_regions_bed_file(self):
        '''test _load_exclude_regions_bed_file'''
        i1 = pyfastaq.intervals.Interval(42, 43)
        i2 = pyfastaq.intervals.Interval(50, 60)
        i3 = pyfastaq.intervals.Interval(100, 102)

        expected = {
            'ref1': [i1, i2],
            'ref2': [i3],
        }

        tmp_bed = 'tmp.load_exclude_regions_bed_file.bed'
        with open(tmp_bed, 'w') as f:
            print('ref1',  43, 44, sep='\t', file=f)
            print('ref1',  51, 55, sep='\t', file=f)
            print('ref1',  56, 61, sep='\t', file=f)
            print('ref2',  101, 103, sep='\t', file=f)

        self.assertEqual({}, mapping_based_verifier.MappingBasedVerifier._load_exclude_regions_bed_file(None))
        got = mapping_based_verifier.MappingBasedVerifier._load_exclude_regions_bed_file(tmp_bed)
        self.assertEqual(expected, got)
        os.unlink(tmp_bed)


    def test_interval_intersects_an_interval_in_list(self):
        '''test _interval_intersects_an_interval_in_list'''
        interval_list = [
            pyfastaq.intervals.Interval(10, 15),
            pyfastaq.intervals.Interval(20, 30),
            pyfastaq.intervals.Interval(40, 50),
            pyfastaq.intervals.Interval(60, 70),
        ]

        tests = [
            (pyfastaq.intervals.Interval(1, 9), False),
            (pyfastaq.intervals.Interval(1, 10), True),
            (pyfastaq.intervals.Interval(1, 20), True),
            (pyfastaq.intervals.Interval(11, 12), True),
            (pyfastaq.intervals.Interval(16, 19), False),
            (pyfastaq.intervals.Interval(15, 75), True),
            (pyfastaq.intervals.Interval(71, 100), False),
        ]

        for interval, true_or_false in tests:
            self.assertEqual(true_or_false, mapping_based_verifier.MappingBasedVerifier._interval_intersects_an_interval_in_list(interval, interval_list))


    def test_filter_vcf_for_clustering(self):
        '''test _filter_vcf_for_clustering'''
        vcf_in = os.path.join(data_dir, 'filter_vcf_for_clustering.in.vcf')
        expected_vcf = os.path.join(data_dir, 'filter_vcf_for_clustering.expect.vcf')
        tmp_out = 'tmp.filter_vcf_for_clustering.out.vcf'
        mapping_based_verifier.MappingBasedVerifier._filter_vcf_for_clustering(vcf_in, tmp_out)
        self.assertTrue(filecmp.cmp(expected_vcf, tmp_out, shallow=False))
        os.unlink(tmp_out)


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


    def test_check_if_sam_match_is_good(self):
        '''test _check_if_sam_match_is_good'''
        ref_fasta = os.path.join(data_dir, 'check_if_sam_match_is_good.ref.fa')
        ref_seqs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_seqs)
        sam_in = os.path.join(data_dir, 'check_if_sam_match_is_good.ref.sam')
        samfile = pysam.AlignmentFile(sam_in, "r")
        for sam_record in samfile.fetch(until_eof=True):
            expected = {'yes': True, 'no': False}[sam_record.query_name.split('.')[-1]]
            self.assertEqual(expected, mapping_based_verifier.MappingBasedVerifier._check_if_sam_match_is_good(sam_record, ref_seqs, 29, allow_mismatches=True))

        for sam_record in samfile.fetch(until_eof=True):
            expected = {'yes': True, 'no': False}[sam_record.query_name.split('.')[-1]]
            self.assertEqual('no', mapping_based_verifier.MappingBasedVerifier._check_if_sam_match_is_good(sam_record, ref_seqs, 29, allow_mismatches=False))


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


    def test_run_no_filter_cluster(self):
        '''test run without filtering and clustering'''
        vcf_file_in = os.path.join(data_dir, 'run.calls.vcf')
        vcf_reference_file = os.path.join(data_dir, 'run.ref.fa')
        verify_reference_file = os.path.join(data_dir, 'run.ref.mutated.fa')
        tmp_out = 'tmp.mapping_based_verifier.out.no_filter_cluster'
        verifier = mapping_based_verifier.MappingBasedVerifier(vcf_file_in, vcf_reference_file, verify_reference_file, tmp_out, flank_length=31, filter_and_cluster_vcf=False)
        verifier.run()
        expected_out = os.path.join(data_dir, 'run.out.no_filter_cluster')
        for suffix in ['.false_negatives.vcf', '.stats.tsv', '.vcf', '.gt_conf_hist.TP.tsv', '.gt_conf_hist.FP.tsv']:
            expected_file = expected_out + suffix
            got_file = tmp_out + suffix
            self.assertTrue(filecmp.cmp(expected_file, got_file, shallow=False))
            os.unlink(got_file)
        samfile = tmp_out + '.sam'
        self.assertTrue(os.path.exists(samfile))
        os.unlink(samfile)
        for suffix in ['.dnadiff.merged.vcf', '.dnadiff.qdiff', '.dnadiff.raw.vcf', '.dnadiff.snps']:
            os.unlink(tmp_out + suffix)


    def test_run_no_filter_cluster_with_exclude_region(self):
        '''test run without filtering and clustering but with region excluded'''
        vcf_file_in = os.path.join(data_dir, 'run.calls.vcf')
        vcf_reference_file = os.path.join(data_dir, 'run.ref.fa')
        verify_reference_file = os.path.join(data_dir, 'run.ref.mutated.fa')
        exclude_regions_bed_file = os.path.join(data_dir, 'run.exclude.bed')
        tmp_out = 'tmp.mapping_based_verifier.out.no_filter_cluster'
        verifier = mapping_based_verifier.MappingBasedVerifier(vcf_file_in, vcf_reference_file, verify_reference_file, tmp_out, flank_length=31, filter_and_cluster_vcf=False, exclude_regions_bed_file=exclude_regions_bed_file)
        verifier.run()
        expected_out = os.path.join(data_dir, 'run.out.no_filter_cluster_with_exclude')
        for suffix in ['.false_negatives.vcf', '.stats.tsv', '.vcf', '.gt_conf_hist.TP.tsv', '.gt_conf_hist.FP.tsv']:
            expected_file = expected_out + suffix
            got_file = tmp_out + suffix
            self.assertTrue(filecmp.cmp(expected_file, got_file, shallow=False))
            os.unlink(got_file)
        samfile = tmp_out + '.sam'
        self.assertTrue(os.path.exists(samfile))
        os.unlink(samfile)
        for suffix in ['.dnadiff.merged.vcf', '.dnadiff.qdiff', '.dnadiff.raw.vcf', '.dnadiff.snps']:
            os.unlink(tmp_out + suffix)


    def test_run_with_filter_cluster(self):
        '''test run with filtering and clustering'''
        vcf_file_in = os.path.join(data_dir, 'run.calls.vcf')
        vcf_reference_file = os.path.join(data_dir, 'run.ref.fa')
        verify_reference_file = os.path.join(data_dir, 'run.ref.mutated.fa')
        tmp_out = 'tmp.mapping_based_verifier.out.with_filter_cluster'
        verifier = mapping_based_verifier.MappingBasedVerifier(vcf_file_in, vcf_reference_file, verify_reference_file, tmp_out, flank_length=31, filter_and_cluster_vcf=True)
        verifier.run()
        expected_out = os.path.join(data_dir, 'run.out.with_filter_cluster')
        for suffix in ['.false_negatives.vcf', '.stats.tsv', '.vcf', '.gt_conf_hist.TP.tsv', '.gt_conf_hist.FP.tsv']:
            expected_file = expected_out + suffix
            got_file = tmp_out + suffix
            self.assertTrue(filecmp.cmp(expected_file, got_file, shallow=False))
            os.unlink(got_file)
        samfile = tmp_out + '.sam'
        self.assertTrue(os.path.exists(samfile))
        os.unlink(samfile)
        for suffix in ['.dnadiff.merged.vcf', '.dnadiff.qdiff', '.dnadiff.raw.vcf', '.dnadiff.snps', '.filter.vcf', '.filter.cluster.vcf']:
            os.unlink(tmp_out + suffix)


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


