import copy
from operator import attrgetter

import pyfastaq
import pymummer
from cluster_vcf_records import vcf_record

from minos import utils

class Error (Exception): pass


class Dnadiff:
    def __init__(self, ref_fasta, query_fasta, outprefix, query_seqs=None):
        self.ref_fasta = ref_fasta
        self.query_fasta = query_fasta
        self.outprefix = outprefix
        if query_seqs is None:
            self.query_seqs = {}
            pyfastaq.tasks.file_to_dict(self.query_fasta, self.query_seqs)
        else:
            self.query_seqs = query_seqs


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
    def _load_snps_file(cls, snps_file, query_seqs):
        '''Loads the .snps file made by dnadiff.
        query_seqs = dictionary of the query sequences.
        ref_seqs, qry_seqs = dictionaries of the genome sequences.
        Returns a dictionary of qry name -> list of VCF records'''
        vcf_records = {}
        variants = pymummer.snp_file.get_all_variants(snps_file)

        for variant in variants:
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
                ]))
            elif variant.var_type == pymummer.variant.DEL:
                # The query has sequence missing, compared to the
                # reference. We're making VCF records w.r.t. the
                # query, so this is an insertion. So need to
                # get the nucleotide before the insertion as well.
                new_record = vcf_record.VcfRecord('\t'.join([
                    variant.qry_name,
                    str(variant.qry_start),
                    '.',
                    query_seqs[variant.qry_name][variant.qry_start - 1],
                    query_seqs[variant.qry_name][variant.qry_start - 1] + variant.ref_base,
                    '.',
                    '.',
                    'SVTYPE=DNADIFF_INS',
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
                ]))
            else:
                raise Error('Unknown variant type: ' + str(variant))

            if new_record.CHROM not in vcf_records:
                vcf_records[new_record.CHROM] = []

            vcf_records[new_record.CHROM].append(new_record)

        for vcf_list in vcf_records.values():
            vcf_list.sort(key=attrgetter('POS'))

        return vcf_records


    @classmethod
    def _make_all_variants_intervals(cls, variants, big_variant_intervals):
        '''Makes union of all positions where there are variants,
        by combining output of _load_qdiff_file() and _load_snps_file()'''
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
        Dnadiff._run_dnadiff(self.ref_fasta, self.query_fasta, self.outprefix)
        self.variants = Dnadiff._load_snps_file(self.outprefix + '.snps', self.query_seqs)
        self.big_variant_intervals = Dnadiff._load_qdiff_file(self.outprefix + '.qdiff')
        self.all_variant_intervals = Dnadiff._make_all_variants_intervals(self.variants, self.big_variant_intervals)

