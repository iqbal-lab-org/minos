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


    def run(self):
        Dnadiff._run_dnadiff(self.ref_fasta, self.query_fasta, self.outprefix)
        self.variants = Dnadiff._load_snps_file(self.outprefix + '.snps', self.query_seqs)
        self.big_variant_intervals = Dnadiff._load_qdiff_file(self.outprefix + '.qdiff')

