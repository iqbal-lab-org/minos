import copy
import os
import random
import shutil
import unittest

import pyfastaq
from cluster_vcf_records import vcf_file_read, vcf_record

from minos import dnadiff


def write_fasta_files(ref_fa, qry_fa):
    random.seed(1)
    nucleotides = ['A', 'C', 'G', 'T']
    ref_seq = [random.choice(nucleotides) for _ in range(5000)]
    # add in one of each type of variant that is
    # reported by dnadiff in the {q,r}diff files
    qry_seq = copy.copy(ref_seq)
    # gap in each of the seqs (will be in .qdiff)
    ref_seq_deletion = ''.join(ref_seq[999:1500])
    ref_seq = ref_seq[:1000] + ref_seq[1500:]
    qry_seq_deletion = ''.join(qry_seq[3999:4500])
    qry_seq = qry_seq[:4000] + qry_seq[4500:]
    ref_seqs = {'ref.gap': pyfastaq.sequences.Fasta('ref.gap', ''.join(ref_seq))}
    qry_seqs = {'qry.gap': pyfastaq.sequences.Fasta('qry.gap', ''.join(qry_seq))}

    # duplication in one of the seqs (will be in .qdiff)
    random.seed(2)
    repeat_seq1 = [random.choice(nucleotides) for _ in range(100)]
    repeat_seq2 = [random.choice(nucleotides) for _ in range(100)]
    left_flank = [random.choice(nucleotides) for _ in range(500)]
    middle_bit = [random.choice(nucleotides) for _ in range(500)]
    right_flank = [random.choice(nucleotides) for _ in range(500)]
    ref_seq = left_flank + repeat_seq1 + middle_bit + 10 * repeat_seq2 + right_flank
    qry_seq = left_flank + 2 * repeat_seq1 + middle_bit + repeat_seq2 + right_flank
    ref_seqs['ref.dup'] = pyfastaq.sequences.Fasta('ref.dup', ''.join(ref_seq))
    qry_seqs['qry.dup'] = pyfastaq.sequences.Fasta('qry.dup', ''.join(qry_seq))

    # long indel of 100bp / 200bp (will be in .qdiff)
    random.seed(3)
    ref_seq = [random.choice(nucleotides) for _ in range(5000)]
    insertion_seq = [random.choice(nucleotides) for _ in range(100)]
    qry_seq = ref_seq[:2500] + insertion_seq + ref_seq[2700:]
    ref_seqs['ref.indel'] = pyfastaq.sequences.Fasta('ref.indel', ''.join(ref_seq))
    qry_seqs['qry.indel'] = pyfastaq.sequences.Fasta('qry.indel', ''.join(qry_seq))

    # snp and short indel (should not end up in qdiff file,
    # but instead in .snps file)
    random.seed(4)
    ref_seq = [random.choice(nucleotides) for _ in range(500)]
    ref_seq[99] = 'A'
    qry_seq = ref_seq[:195] + ref_seq[200:400] + list('ACGT') + ref_seq[400:]
    ref_seq[99] = 'G'
    ref_seqs['ref.snp_indel'] = pyfastaq.sequences.Fasta('ref.snp_indel', ''.join(ref_seq))
    qry_seqs['qry.snp_indel'] = pyfastaq.sequences.Fasta('qry.snp_indel', ''.join(qry_seq))

    # Another pair of sequences with a SNP. Just to check that
    # the dictionary made by _load_snps_file is OK
    random.seed(5)
    ref_seq = [random.choice(nucleotides) for _ in range(500)]
    ref_seq[249] = 'A'
    qry_seq = copy.copy(ref_seq)
    ref_seq[249] = 'G'
    ref_seqs['ref.snp_indel.2'] = pyfastaq.sequences.Fasta('ref.snp_indel.2', ''.join(ref_seq))
    qry_seqs['qry.snp_indel.2'] = pyfastaq.sequences.Fasta('qry.snp_indel.2', ''.join(qry_seq))

    # Test when query is reverse complemented.
    # When this happens, the snps file from mummer has -1 in column 10, and it
    # reports the reverse complement of the qry sequence. eg here we're adding a SNP at
    # (1-based) position 100 in the ref, changing the G to a A. But then reverse complementing
    # means there is a T at that position in the query sequence. The SNP file reports ref -> qry change as
    # G -> A, but we're making VCF in terms of the query. We want T -> C.
    random.seed(6)
    ref_seq = [random.choice(nucleotides) for _ in range(500)]
    ref_seq[99] = 'A'
    ref_seq[195] = 'A'
    ref_seq[196] = 'A'
    ref_seq[197] = 'G'
    ref_seq[198] = 'A'
    ref_seq[199] = 'G'
    ref_seq[200] = 'C'
    qry_seq = ref_seq[:195] + ref_seq[200:400] + list('AGTC') + ref_seq[400:]
    ref_seq[99] = 'G'
    ref_seqs['ref.snp_indel_qry_rev'] = pyfastaq.sequences.Fasta('ref.snp_indel_qry_rev', ''.join(ref_seq))
    qry_seqs['qry.snp_indel_qry_rev'] = pyfastaq.sequences.Fasta('qry.snp_indel_qry_rev', ''.join(qry_seq))
    qry_seqs['qry.snp_indel_qry_rev'].revcomp()

    # Test when ref is reverse complemented. This is like the query being revcomped, tested
    # previously..Turns out that mummer always has the ref on the forward strand.
    # This SNP is effectively C in the ref, changed to a T in the query. We see
    # C -> T in the dnadiff file
    random.seed(7)
    ref_seq = [random.choice(nucleotides) for _ in range(500)]
    ref_seq[99] = 'A'
    qry_seq = ref_seq[:195] + ref_seq[200:400] + list('AGTC') + ref_seq[400:]
    ref_seq[99] = 'G'
    ref_seqs['ref.snp_indel_ref_rev'] = pyfastaq.sequences.Fasta('ref.snp_indel_ref_rev', ''.join(ref_seq))
    ref_seqs['ref.snp_indel_ref_rev'].revcomp()
    qry_seqs['qry.snp_indel_ref_rev'] = pyfastaq.sequences.Fasta('qry.snp_indel_ref_rev', ''.join(qry_seq))


    name_suffixes = ['gap', 'dup', 'indel', 'snp_indel', 'snp_indel.2', 'snp_indel_qry_rev', 'snp_indel_ref_rev']
    with open(ref_fa, 'w') as f_ref, open(qry_fa, 'w') as f_qry:
        for suffix in name_suffixes:
            print(ref_seqs['ref.' + suffix], file=f_ref)
            print(qry_seqs['qry.' + suffix], file=f_qry)


    # This is confusing with all the reverse complementing. Made perfect reads from the ref, and used them
    # to call variants with bwa+samtools against the query. Results of that used to make the
    # following expected records. Results were:
    # qry.snp_indel	100	.	A	G,<*>	0	.	DP=28;I16=0,0,10,15,0,0,1000,40000,0,0,1500,90000,0,0,379,6411;QS=0,1,0;VDB=0.0364506;SGB=-0.692914;MQSB=1;MQ0F=0	PL	255,75,0,255,75,255
    #qry.snp_indel	195	.	A	AGATTC	0	.	INDEL;IDV=13;IMF=0.65;DP=20;I16=0,0,8,6,0,0,1800,241472,0,0,786,44244,0,0,199,3323;QS=0,1;VDB=0.193256;SGB=-0.686358;MQSB=0.991701;MQ0F=0	PL	255,42,0
    #qry.snp_indel	395	.	AACGT	A	0	.	INDEL;IDV=29;IMF=0.90625;DP=32;I16=0,0,13,16,0,0,2500,216028,0,0,1718,101804,0,0,456,7930;QS=0,1;VDB=0.0162957;SGB=-0.693079;MQSB=0.839255;MQ0F=0	PL	255,87,0
    #qry.snp_indel.2	250	.	A	G,<*>	0	.	DP=40;I16=0,0,16,20,0,0,1388,54342,0,0,2160,129600,0,0,471,7875;QS=0,1,0;VDB=0.328997;SGB=-0.693139;MQSB=1;MQ0F=0	PL	255,108,0,255,108,255
    #qry.snp_indel_qry_rev	100	.	GGACT	G	0	.	INDEL;IDV=17;IMF=0.772727;DP=22;I16=0,0,8,9,0,0,1547,140777,0,0,1002,59076,0,0,311,6025;QS=0,1;VDB=0.000460773;SGB=-0.690438;MQSB=0.981652;MQ0F=0	PL	255,51,0
    #qry.snp_indel_qry_rev	304	.	G	GCTCTT	0	.	INDEL;IDV=10;IMF=0.769231;DP=13;I16=0,0,10,1,0,0,1392,185984,0,0,606,33444,0,0,161,2777;QS=0,1;VDB=0.184302;SGB=-0.676189;MQSB=1;MQ0F=0	PL	255,33,0
    #qry.snp_indel_qry_rev	400	.	T	C,<*>	0	.	DP=33;I16=0,0,13,17,0,0,1200,48000,0,0,1800,108000,0,0,468,8084;QS=0,1,0;VDB=0.010605;SGB=-0.693097;MQSB=1;MQ0F=0PL	255,90,0,255,90,255
    #qry.snp_indel_ref_rev	100	.	A	G,<*>	0	.	DP=29;I16=0,0,15,10,0,0,978,38554,0,0,1500,90000,0,0,313,4825;QS=0,1,0;VDB=0.540121;SGB=-0.692914;MQSB=1;MQ0F=0	PL	255,75,0,255,75,255
    #qry.snp_indel_ref_rev	195	.	TG	TGCGTGG	0	.	INDEL;IDV=16;IMF=0.761905;DP=21;I16=0,0,9,9,0,0,2208,289344,0,0,1008,56592,0,0,253,4165;QS=0,1;VDB=0.184493;SGB=-0.691153;MQSB=0.28276;MQ0F=0	PL	255,54,0
    #qry.snp_indel_ref_rev	395	.	TAGTC	T	0	.	INDEL;IDV=22;IMF=0.709677;DP=31;I16=0,0,16,6,0,0,1914,166518,0,0,1300,76840,0,0,395,7577;QS=0,1;VDB=0.000139106;SGB=-0.692562;MQSB=0.97584;MQ0F=0	PL	255,66,0
    expected_vcf_records = {
        'qry.snp_indel': [
            vcf_record.VcfRecord('qry.snp_indel\t100\t.\tA\tG\t.\t.\tSVTYPE=DNADIFF_SNP\tGT\t1/1'),
            vcf_record.VcfRecord('qry.snp_indel\t195\t.\tA\tAGATTC\t.\t.\tSVTYPE=DNADIFF_INS\tGT\t1/1'),
            vcf_record.VcfRecord('qry.snp_indel\t395\t.\tAACGT\tA\t.\t.\tSVTYPE=DNADIFF_DEL\tGT\t1/1'),
        ],
        'qry.snp_indel.2': [
            vcf_record.VcfRecord('qry.snp_indel.2\t250\t.\tA\tG\t.\t.\tSVTYPE=DNADIFF_SNP\tGT\t1/1'),
        ],
        'qry.snp_indel_qry_rev' : [
            vcf_record.VcfRecord('qry.snp_indel_qry_rev\t100\t.\tGGACT\tG\t.\t.\tSVTYPE=DNADIFF_DEL\tGT\t1/1'),
            vcf_record.VcfRecord('qry.snp_indel_qry_rev\t304\t.\tG\tGCTCTT\t.\t.\tSVTYPE=DNADIFF_INS\tGT\t1/1'),
            vcf_record.VcfRecord('qry.snp_indel_qry_rev\t400\t.\tT\tC\t.\t.\tSVTYPE=DNADIFF_SNP\tGT\t1/1'),
        ],
        'qry.snp_indel_ref_rev' : [
            vcf_record.VcfRecord('qry.snp_indel_ref_rev\t100\t.\tA\tG\t.\t.\tSVTYPE=DNADIFF_SNP\tGT\t1/1'),
            vcf_record.VcfRecord('qry.snp_indel_ref_rev\t196\t.\tG\tGCGTGG\t.\t.\tSVTYPE=DNADIFF_INS\tGT\t1/1'),
            vcf_record.VcfRecord('qry.snp_indel_ref_rev\t395\t.\tTAGTC\tT\t.\t.\tSVTYPE=DNADIFF_DEL\tGT\t1/1'),
        ]
    }

    expected_regions = {
        'qry.indel': [pyfastaq.intervals.Interval(2500, 2599)],
        'qry.gap': [
            pyfastaq.intervals.Interval(1000, 1498),
            pyfastaq.intervals.Interval(3999, 4000),
        ],
        'qry.dup': [
            pyfastaq.intervals.Interval(599, 603),
            pyfastaq.intervals.Interval(1199, 1300),
        ],
    }

    return ref_seqs, qry_seqs, expected_vcf_records, expected_regions


dnadiff_output_extensions = [
    '1coords',
    '1delta',
    'delta',
    'mcoords',
    'mdelta',
    'qdiff',
    'rdiff',
    'report',
    'snps',
]


modules_dir = os.path.dirname(os.path.abspath(dnadiff.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'dnadiff')

class TestDnadiff(unittest.TestCase):
    def test_load_seq_file(self):
        '''test _load_seq_file'''
        seq1 = pyfastaq.sequences.Fasta('seq.1', 'ACGTA')
        seq2 = pyfastaq.sequences.Fasta('seq.2', 'AAA')
        tmp_file = 'tmp.dnadiff.test_load_seq_file.fa'
        with open(tmp_file, 'w') as f:
            print(seq1, file=f)
            print(seq2, file=f)
        got_names, got_seqs = dnadiff.Dnadiff._load_seq_file(tmp_file)
        expected_names = ['seq.1', 'seq.2']
        expected_seqs = {'seq.1': seq1, 'seq.2': seq2}
        self.assertEqual(expected_names, got_names)
        self.assertEqual(expected_seqs, got_seqs)
        os.unlink(tmp_file)


    def test_run_dnadiff(self):
        '''test _run_dnadiff'''
        # Just test that dnadiff runs and we get the
        # output files we need. Not check their contents.
        ref_fasta = os.path.join(data_dir, 'run_dnadiff.ref.fa')
        qry_fasta = os.path.join(data_dir, 'run_dnadiff.qry.fa')
        tmp_prefix = 'tmp.run_dnadiff'
        dnadiff.Dnadiff._run_dnadiff(ref_fasta, qry_fasta, tmp_prefix)
        self.assertTrue(os.path.exists(tmp_prefix + '.snps'))
        self.assertTrue(os.path.exists(tmp_prefix + '.qdiff'))
        dnadiff.Dnadiff.clean_dnadiff_files(tmp_prefix)



    def test_load_qdiff_file(self):
        '''test _load_qdiff_file'''
        # Easy (for me at least!) to get confused between
        # insertions and deletions bewteen ref and query.
        # Generate ref and qry sequences here, where we know the differences
        # and then run dnadiff to make diff files.
        # Then test the actual loading of them.
        # Note: dnadiff does not always get the coords exact, so
        # when we make the expceted coords, this is made from looking at the
        # output of dnadiff (manually checking that it's "close enough").
        ref_fa = 'tmp.test_load_qdiff_file.ref.fa'
        qry_fa = 'tmp.test_load_qdiff_file.qry.fa'
        ref_seqs, qry_seqs, expected_vcf_records, expected_regions = write_fasta_files(ref_fa, qry_fa)

        outprefix = 'tmp.test_load_qdiff_file.dnadiff'
        dnadiff.Dnadiff._run_dnadiff(ref_fa, qry_fa, outprefix)
        qdiff_file = outprefix + '.qdiff'
        got_regions = dnadiff.Dnadiff._load_qdiff_file(qdiff_file)
        self.assertEqual(expected_regions, got_regions)

        dnadiff.Dnadiff.clean_dnadiff_files(outprefix)
        os.unlink(ref_fa)
        os.unlink(qry_fa)


    def test_snps_file_file_to_unmerged_vcf(self):
        '''test _snps_file_file_to_unmerged_vcf'''
        ref_fa = 'tmp.test_snps_file_file_to_unmerged_vcf.ref.fa'
        qry_fa = 'tmp.test_snps_file_file_to_unmerged_vcf.qry.fa'
        ref_seqs, qry_seqs, expected_vcf_records, expected_regions = write_fasta_files(ref_fa, qry_fa)
        outprefix = 'tmp.test_snps_file_file_to_unmerged_vcf.dnadiff'
        vcf_out = outprefix + '.out.vcf'
        dnadiff.Dnadiff._run_dnadiff(ref_fa, qry_fa, outprefix)
        dnadiff.Dnadiff._snps_file_file_to_unmerged_vcf(outprefix + '.snps', qry_seqs, vcf_out)
        header, got = vcf_file_read.vcf_file_to_dict(vcf_out)
        self.assertTrue(expected_vcf_records, got)

        dnadiff.Dnadiff.clean_dnadiff_files(outprefix)
        os.unlink(ref_fa)
        os.unlink(qry_fa)
        os.unlink(vcf_out)


    def test_make_all_variants_intervals(self):
        '''test _make_all_variants_intervals'''
        variants = {
            'seq.1': [
                vcf_record.VcfRecord('seq.1\t15\t.\tAGTTGTC\tA\t.\t.\tSVTYPE=DEL'),
                vcf_record.VcfRecord('seq.1\t100\t.\tT\tA\t.\t.\tSVTYPE=SNP'),
            ],
            'seq.2': [
                vcf_record.VcfRecord('seq.1\t43\t.\tA\tACGTA\t.\t.\tSVTYPE=INS'),
            ],
        }
        big_variant_intervals = {
            'seq.1': [
                pyfastaq.intervals.Interval(9, 19),
                pyfastaq.intervals.Interval(50, 60),
            ],
            'seq.3': [
                pyfastaq.intervals.Interval(42, 45),
            ],
        }
        got = dnadiff.Dnadiff._make_all_variants_intervals(variants, big_variant_intervals)
        expected = {
            'seq.1': [
                pyfastaq.intervals.Interval(9, 20),
                pyfastaq.intervals.Interval(50, 60),
                pyfastaq.intervals.Interval(99, 99),
            ],
            'seq.2': [
                pyfastaq.intervals.Interval(42, 42),
            ],
            'seq.3': [
                pyfastaq.intervals.Interval(42, 45),
            ],
        }
        self.assertEqual(expected, got)


    def test_run(self):
        '''test run'''
        ref_fa = 'tmp.test_run_file.ref.fa'
        qry_fa = 'tmp.test_run_file.qry.fa'
        ref_seqs, qry_seqs, expected_vcf_records, expected_regions = write_fasta_files(ref_fa, qry_fa)
        outprefix = 'tmp.test_run_file.dnadiff'
        dnadiffer = dnadiff.Dnadiff(ref_fa, qry_fa, outprefix)
        dnadiffer.run()
        self.assertEqual(expected_vcf_records['qry.snp_indel_ref_rev'], dnadiffer.variants['qry.snp_indel_ref_rev'])
        self.assertEqual(expected_regions, dnadiffer.big_variant_intervals)
        expected_all_variant_intervals = dnadiff.Dnadiff._make_all_variants_intervals(dnadiffer.variants, dnadiffer.big_variant_intervals)
        self.assertEqual(expected_all_variant_intervals, dnadiffer.all_variant_intervals)
        dnadiff.Dnadiff.clean_dnadiff_files(outprefix)
        os.unlink(ref_fa)
        os.unlink(qry_fa)
        os.unlink(dnadiffer.unmerged_vcf)
        os.unlink(dnadiffer.merged_vcf)
