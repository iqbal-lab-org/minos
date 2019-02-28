import copy
from operator import attrgetter
import os

import pyfastaq
import pymummer
from cluster_vcf_records import vcf_clusterer, vcf_file_read, vcf_record

from minos import utils

class Error (Exception): pass

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
    'unqry',
    'unref',
]

class Dnadiff:
    def __init__(self, ref_fasta, query_fasta, outprefix):
        self.ref_fasta = ref_fasta
        self.query_fasta = query_fasta
        self.outprefix = outprefix
        self.ref_seq_names, self.ref_seqs = Dnadiff._load_seq_file(self.ref_fasta)
        self.query_seq_names, self.query_seqs = Dnadiff._load_seq_file(self.query_fasta)


    @classmethod
    def clean_dnadiff_files(cls, prefix):
        for extension in dnadiff_output_extensions:
            # not all files get written, hence try except pass
            try:
                os.unlink(prefix + '.' + extension)
            except:
                pass


    @classmethod
    def _load_seq_file(cls, infile):
        '''Returns tuple of:
        1. list of seqnames in same order as in file, and
        2. dictionary of seqs'''
        seqreader = pyfastaq.sequences.file_reader(infile)
        seqs = {}
        seqnames = []
        for seq in seqreader:
            seqnames.append(seq.id)
            seqs[seq.id] = copy.copy(seq)

        return seqnames, seqs


    @classmethod
    def _run_dnadiff(cls, ref_fasta, query_fasta, outprefix):
        command = ' '.join(['dnadiff -p', outprefix, ref_fasta, query_fasta])
        utils.syscall(command)


    @classmethod
    def _load_qdiff_file(cls, qdiff_file):
        '''Loads the out.qdiff file made by dnadiff.
        Returns dict of seq name -> list of sorted non-overlapping
        regions that contain variants compared to the ref'''
        diff_regions = {}

        with open(qdiff_file) as f:
            for line in f:
                fields = line.rstrip().split()
                start, end = sorted([int(fields[2]) - 1, int(fields[3]) - 1])
                interval = pyfastaq.intervals.Interval(start, end)
                if fields[0] not in diff_regions:
                    diff_regions[fields[0]] = []
                diff_regions[fields[0]].append(interval)

        for interval_list in diff_regions.values():
            pyfastaq.intervals.merge_overlapping_in_list(interval_list)

        return diff_regions


    @classmethod
    def _snps_file_file_to_unmerged_vcf(cls, snps_file, query_seqs, outfile):
        '''Loads the .snps file made by dnadiff.
        query_seqs = dictionary of the query sequences.
        ref_seqs, qry_seqs = dictionaries of the genome sequences.
        Writes a new VCF file unmerged records.'''
        vcf_records = {}
        variants = pymummer.snp_file.get_all_variants(snps_file)

        for variant in variants:
            # If the variant is reversed, it means that either the ref or query had to be
            # reverse complemented when aligned by mummer. Need to do the appropriate
            # reverse (complement) fixes so the VCF has the correct REF and ALT sequences
            if variant.reverse:
                qry_seq = pyfastaq.sequences.Fasta('x', variant.qry_base)
                qry_seq.revcomp()
                variant.qry_base = ''.join(reversed(qry_seq.seq))
                ref_seq = pyfastaq.sequences.Fasta('x', variant.ref_base)
                ref_seq.revcomp()
                variant.ref_base = ref_seq.seq

            if variant.var_type == pymummer.variant.SNP:
                new_record = vcf_record.VcfRecord('\t'.join([
                    variant.qry_name,
                    str(variant.qry_start + 1),
                    '.',
                    variant.qry_base,
                    variant.ref_base,
                    '.',
                    '.',
                    'SVTYPE=DNADIFF_SNP',
                    'GT',
                    '1/1',
                ]))
            elif variant.var_type == pymummer.variant.DEL:
                # The query has sequence missing, compared to the
                # reference. We're making VCF records w.r.t. the
                # query, so this is an insertion. So need to
                # get the nucleotide before the insertion as well.
                new_record = vcf_record.VcfRecord('\t'.join([
                    variant.qry_name,
                    str(variant.qry_start + 1),
                    '.',
                    query_seqs[variant.qry_name][variant.qry_start],
                    query_seqs[variant.qry_name][variant.qry_start] + variant.ref_base,
                    '.',
                    '.',
                    'SVTYPE=DNADIFF_INS',
                    'GT',
                    '1/1',
                ]))
            elif variant.var_type == pymummer.variant.INS:
                # The ref has sequence missing, compared to the
                # query. We're making VCF records w.r.t. the
                # query, so this is a deletion. So need to
                # get the nucleotide before the deletion as well.
                new_record = vcf_record.VcfRecord('\t'.join([
                    variant.qry_name,
                    str(variant.qry_start),
                    '.',
                    query_seqs[variant.qry_name][variant.qry_start - 1] + variant.qry_base,
                    query_seqs[variant.qry_name][variant.qry_start - 1],
                    '.',
                    '.',
                    'SVTYPE=DNADIFF_DEL',
                    'GT',
                    '1/1',
                ]))
            else:
                raise Error('Unknown variant type: ' + str(variant))

            assert new_record.REF == query_seqs[new_record.CHROM][new_record.POS:new_record.POS + len(new_record.REF)]

            if new_record.CHROM not in vcf_records:
                vcf_records[new_record.CHROM] = []

            vcf_records[new_record.CHROM].append(new_record)

        for vcf_list in vcf_records.values():
            vcf_list.sort(key=attrgetter('POS'))

        with open(outfile, 'w') as f:
            for key, vcf_list in sorted(vcf_records.items()):
                for record in vcf_list:
                    print(record, file=f)


    @classmethod
    def _make_all_variants_intervals(cls, variants, big_variant_intervals):
        '''Makes union of all positions where there are variants,
        by combining output of _load_qdiff_file() and variants from merged VCF file'''
        all_ref_seqs = set(variants.keys()).union(set(big_variant_intervals.keys()))
        intervals = {x: [] for x in all_ref_seqs}

        for ref_seq in intervals:
            if ref_seq in variants:
                for variant in variants[ref_seq]:
                    if len(variant.REF) > 1 or len(variant.ALT) > 1:
                        intervals[ref_seq].append(pyfastaq.intervals.Interval(variant.POS + 1, variant.ref_end_pos()))
                    else:
                        intervals[ref_seq].append(pyfastaq.intervals.Interval(variant.POS, variant.ref_end_pos()))

            if ref_seq in big_variant_intervals:
                intervals[ref_seq].extend(big_variant_intervals[ref_seq])

            pyfastaq.intervals.merge_overlapping_in_list(intervals[ref_seq])

        return intervals


    def run(self):
        snps_file = self.outprefix + '.snps'
        qdiff_file = self.outprefix + '.qdiff'
        self.unmerged_vcf = self.outprefix + '.raw.vcf'
        self.merged_vcf = self.outprefix + '.merged.vcf'
        for filename in [snps_file, qdiff_file]:
            if os.path.exists(filename):
                os.unlink(filename)
        tmp_prefix = self.outprefix + '.tmp'

        for ref_name, query_name in zip(self.ref_seq_names, self.query_seq_names):
            ref_fasta = tmp_prefix + '.ref.fa'
            query_fasta = tmp_prefix + '.query.fa'
            with open(ref_fasta, 'w') as f:
                print(self.ref_seqs[ref_name], file=f)
            with open(query_fasta, 'w') as f:
                print(self.query_seqs[query_name], file=f)
            Dnadiff._run_dnadiff(ref_fasta, query_fasta, tmp_prefix)
            utils.syscall('cat ' + tmp_prefix + '.snps >> ' + snps_file)
            utils.syscall('cat ' + tmp_prefix + '.qdiff >> ' + qdiff_file)
            Dnadiff.clean_dnadiff_files(tmp_prefix)
            os.unlink(ref_fasta)
            os.unlink(query_fasta)

        Dnadiff._snps_file_file_to_unmerged_vcf(self.outprefix + '.snps', self.query_seqs, self.unmerged_vcf)
        clusterer = vcf_clusterer.VcfClusterer([self.unmerged_vcf], self.query_fasta, self.merged_vcf, merge_method='simple')
        clusterer.run()
        header, self.variants = vcf_file_read.vcf_file_to_dict(self.merged_vcf, remove_useless_start_nucleotides=True)
        self.big_variant_intervals = Dnadiff._load_qdiff_file(self.outprefix + '.qdiff')
        self.all_variant_intervals = Dnadiff._make_all_variants_intervals(self.variants, self.big_variant_intervals)

