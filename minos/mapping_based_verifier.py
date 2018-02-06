import os

import pyfastaq
import pysam

from cluster_vcf_records import vcf_file_read

from minos import dependencies, utils

class Error (Exception): pass


def sam_reader(infile):
    '''Reads SAM file made by _MappingBasedVerifier._map_seqs_to_ref().
    yields list of all sam records corresponding to each VCF record.'''
    samfile = pysam.AlignmentFile(infile, "r")
    current_records = []
    current_ref_start_record_index = None

    for sam_record in samfile.fetch(until_eof=True):
        ref_name, expected_start, vcf_record_index, allele_index = sam_record.query_name.rsplit('.', maxsplit=3)
        ref_start_record_index = (ref_name, expected_start, vcf_record_index)

        if current_ref_start_record_index in [None, ref_start_record_index]:
            current_records.append(sam_record)
        else:
            yield current_records
            current_records = [sam_record]

        current_ref_start_record_index = ref_start_record_index

    if len(current_records) > 0:
        yield current_records


class MappingBasedVerifier:
    '''vcf_file_in = input VCF file, to be verified.
    vcf_reference_file = reference sequence file that was used to make vcf_file_in.
    verify_reference_file = reference file that we assume is the truth.
    outprefix = prefix of output files.
    flank_length = length of "truth" sequence to take before/after alleles when mapping.

    Writes 3 files:
    outprefix.sam = mapping of seqs + flanks to verify_reference_file.
    outprefix.vcf = VCF file with annotations for validation
    outprefix.stats.tsv = summary stats (see dict output by
                          _parse_sam_file_and_update_vcf_records_and_gather_stats()
                          for a description)'''
    def __init__(self, vcf_file_in, vcf_reference_file, verify_reference_file, outprefix, flank_length=31, expected_variants_vcf=None, run_dnadiff=True):
        self.vcf_file_in = os.path.abspath(vcf_file_in)
        self.vcf_reference_file = os.path.abspath(vcf_reference_file)
        self.verify_reference_file = os.path.abspath(verify_reference_file)
        self.vcf_file_out = os.path.abspath(outprefix + '.vcf')
        self.sam_file_out = os.path.abspath(outprefix + '.sam')
        self.seqs_out = os.path.abspath(outprefix + '.fa')
        self.stats_out = os.path.abspath(outprefix + '.stats.tsv')
        self.flank_length = flank_length
        self.expected_variants_vcf = expected_variants_vcf
        self.run_dnadiff = run_dnadiff
        if self.run_dnadiff and self.expected_variants_vcf is not None:
            raise Error('Error! Incompatible options. expected_variants_vcf file provided, and run_dnadiff is True')


    @classmethod
    def _write_vars_plus_flanks_to_fasta(cls, outfile, vcf_records, ref_seqs, flank_length):
        '''Given a dict of vcf records made by vcf_file_read.vcf_file_to_dict(),
        and its correcsponding file of reference sequences, writes a new fasta file
        of each ref seq and inferred variant sequence plus flank_length nucleotides added to
        its start and end. Calls each sequence:
            ref_name.start_position.vcf_list_index.allele_number
        where allele_numbers in same order as VCF, with ref seq = allele 0.'''
        with open(outfile, 'w') as f:
            for ref_name in sorted(vcf_records):
                for i, vcf_record in enumerate(vcf_records[ref_name]):
                    start_position, alleles = vcf_record.inferred_var_seqs_plus_flanks(ref_seqs[ref_name], flank_length)
                    for allele_index, allele_seq in enumerate(alleles):
                        seq_name = '.'.join([ref_name, str(start_position + 1), str(i), str(allele_index)])
                        print('>' + seq_name, allele_seq, sep='\n', file=f)


    @classmethod
    def _map_seqs_to_ref(cls, seqs_file, ref_file, outfile):
        '''Map seqs_file to ref_file using BWA MEM.
        Output is SAM file written to outfile'''
        bwa_binary = dependencies.find_binary('bwa')
        command = ' '.join([
            bwa_binary, 'mem',
            '-a', # report all mappings
            '-Y', # use soft clipping for supplementary alignments
            ref_file,
            seqs_file,
            '>', outfile,
        ])
        utils.syscall(command)


    @classmethod
    def _parse_sam_file_and_update_vcf_records_and_gather_stats(cls, infile, vcf_records):
        '''Input is SAM file made by _map_seqs_to_ref(), and corresponding dict
        of VCF records made by vcf_file_read.file_to_dict.
        Adds validation info to each VCF record. Returns a dict of stats that
        has keys => values:
            total => total number of VCF records
            ref_pass => number of records where REF allele found in reference
            alt_pass => number of records where (at least) one ALT allele found in reference'''
        stats = {x: 0 for x in ['total', 'ref_pass', 'alt_pass']}
        sreader = sam_reader(infile)

        for sam_list in sreader:
            stats['total'] += 1

            # Split the sam_list (Which has all records for all alleles from one
            # VCf line) into a list per allele.
            sam_records_by_allele = []
            current_allele_index = None
            ref_info_tuples = set()
            for sam_record in sam_list:
                ref_name, expected_start, vcf_record_index, allele_index = sam_record.query_name.rsplit('.', maxsplit=3)
                ref_info_tuples.add((ref_name, expected_start, vcf_record_index))
                if current_allele_index is None or current_allele_index != allele_index:
                    sam_records_by_allele.append([sam_record])
                else:
                    sam_records_by_allele[-1].append(sam_record)
                current_allele_index = allele_index

            assert len(ref_info_tuples) == 1
            ref_name, expected_start, vcf_record_index = ref_info_tuples.pop()
            expected_start = int(expected_start) - 1
            vcf_record_index = int(vcf_record_index)
            results = {x: [] for x in ['MINOS_CHECK_PASS', 'MINOS_CHECK_BEST_HITS', 'MINOS_CHECK_NM', 'MINOS_CHECK_CIGAR', 'MINOS_CHECK_MD']}

            for allele_sam_list in sam_records_by_allele:
                indexes_with_matched_flanks = []
                results['MINOS_CHECK_PASS'].append('0')

                for i, hit in enumerate(allele_sam_list):
                    if hit.is_unmapped:
                        continue

                    try:
                        rs = hit.get_reference_sequence().upper()
                        start_ok = rs[:31]==hit.query_sequence[:31]
                        end_ok = rs[-31:]==hit.query_sequence[-31:]
                    except:
                        continue

                    if start_ok and end_ok:
                        indexes_with_matched_flanks.append(i)

                if len(indexes_with_matched_flanks) > 0:
                    min_nm = min([allele_sam_list[i].get_tag('NM') for i in indexes_with_matched_flanks])
                    best_nm_hits_indexes = [i for i in indexes_with_matched_flanks if allele_sam_list[i].get_tag('NM') == min_nm]

                    for i in best_nm_hits_indexes:
                        if allele_sam_list[i].query_alignment_length == allele_sam_list[i].infer_query_length() and min_nm == 0:
                            results['MINOS_CHECK_PASS'][-1] = '1'
                            break

                    results['MINOS_CHECK_BEST_HITS'].append(str(len(best_nm_hits_indexes)))
                    results['MINOS_CHECK_NM'].append(str(min_nm))
                    cigars = '_'.join([allele_sam_list[i].cigarstring for i in best_nm_hits_indexes])
                    results['MINOS_CHECK_CIGAR'].append(cigars)
                    md = '_'.join([allele_sam_list[i].get_tag('MD') for i in best_nm_hits_indexes])
                    results['MINOS_CHECK_MD'].append(md)
                else:
                    results['MINOS_CHECK_BEST_HITS'].append('0')
                    results['MINOS_CHECK_NM'].append('NA')
                    results['MINOS_CHECK_CIGAR'].append('NA')
                    results['MINOS_CHECK_MD'].append('NA')

            vcf_record = vcf_records[ref_name][vcf_record_index]

            for key in sorted(results):
                vcf_record.set_format_key_value(key, ','.join(results[key]))

            if '1' in results['MINOS_CHECK_PASS'][1:]:
                stats['alt_pass'] += 1
            if results['MINOS_CHECK_PASS'][0] == '1':
                stats['ref_pass'] += 1
        return stats

    def run(self):
        vcf_header, vcf_records = vcf_file_read.vcf_file_to_dict(self.vcf_file_in, sort=True)
        vcf_ref_seqs = {}
        pyfastaq.tasks.file_to_dict(self.vcf_reference_file, vcf_ref_seqs)
        MappingBasedVerifier._write_vars_plus_flanks_to_fasta(self.seqs_out, vcf_records, vcf_ref_seqs, self.flank_length)
        MappingBasedVerifier._map_seqs_to_ref(self.seqs_out, self.verify_reference_file, self.sam_file_out)
        os.unlink(self.seqs_out)
        stats = MappingBasedVerifier._parse_sam_file_and_update_vcf_records_and_gather_stats(self.sam_file_out, vcf_records)

        with open(self.vcf_file_out, 'w') as f:
            print(*vcf_header, sep='\n', file=f)
            for r in vcf_records:
                for v in vcf_records[r]:
                    print(v, file=f)

        with open(self.stats_out, 'w') as f:
            keys = ['total', 'ref_pass', 'alt_pass']
            print(*keys, sep='\t', file=f)
            print(*[stats[x] for x in keys], sep='\t', file=f)

