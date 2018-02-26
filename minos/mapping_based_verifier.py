import logging
import os

import pyfastaq
import pysam

from cluster_vcf_records import vcf_clusterer, vcf_file_read

from minos import dependencies, dnadiff, utils

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
    def __init__(self, vcf_file_in, vcf_reference_file, verify_reference_file, outprefix, flank_length=31, expected_variants_vcf=None, run_dnadiff=True, filter_and_cluster_vcf=True):
        self.vcf_file_in = os.path.abspath(vcf_file_in)
        self.vcf_reference_file = os.path.abspath(vcf_reference_file)
        self.verify_reference_file = os.path.abspath(verify_reference_file)
        self.vcf_file_out = os.path.abspath(outprefix + '.vcf')
        self.sam_file_out = os.path.abspath(outprefix + '.sam')
        self.filtered_vcf = os.path.abspath(outprefix + '.filter.vcf')
        self.clustered_vcf = os.path.abspath(outprefix + '.filter.cluster.vcf')
        self.seqs_out = os.path.abspath(outprefix + '.fa')
        self.stats_out = os.path.abspath(outprefix + '.stats.tsv')
        self.gt_conf_hists_filenames = {
            'TP': os.path.abspath(outprefix + '.gt_conf_hist.TP.tsv'),
            'FP': os.path.abspath(outprefix + '.gt_conf_hist.FP.tsv'),
        }
        self.flank_length = flank_length
        self.expected_variants_vcf = expected_variants_vcf
        self.run_dnadiff = run_dnadiff
        if self.run_dnadiff and self.expected_variants_vcf is not None:
            raise Error('Error! Incompatible options. expected_variants_vcf file provided, and run_dnadiff is True')
        self.dnadiff_outprefix = os.path.abspath(outprefix + '.dnadiff')
        self.vcf_false_negatives_file_out = os.path.abspath(outprefix + '.false_negatives.vcf')
        self.filter_and_cluster_vcf = filter_and_cluster_vcf

        if self.filter_and_cluster_vcf:
            self.vcf_to_check = self.clustered_vcf
        else:
            self.vcf_to_check = self.vcf_file_in


    @classmethod
    def _filter_vcf_for_clustering(cls, infile, outfile):
        header_lines, vcf_records = vcf_file_read.vcf_file_to_dict(infile, sort=True, homozygous_only=False, remove_asterisk_alts=True, remove_useless_start_nucleotides=True)

        with open(outfile, 'w') as f:
            print(*header_lines, sep='\n', file=f)
            for ref_name in vcf_records:
                for vcf_record in vcf_records[ref_name]:
                    if vcf_record.FILTER == 'MISMAPPED_UNPLACEABLE':
                        continue
                    if vcf_record.FORMAT is None or 'GT' not in vcf_record.FORMAT:
                        logging.warning('No GT in vcf record:' + str(vcf_record))
                        continue

                    genotype = vcf_record.FORMAT['GT']
                    genotypes = genotype.split('/')
                    called_alleles = set(genotypes)
                    if len(called_alleles) != 1 or called_alleles == {'0'} or '.' in called_alleles:
                        continue

                    vcf_record.set_format_key_value('GT', '1/1')

                    try:
                        vcf_record.ALT = [vcf_record.ALT[int(genotypes[0]) - 1]]
                    except:
                        raise Error('BAD VCf line:' + str(vcf_record))

                    print(vcf_record, file=f)


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
    def _check_called_genotype(cls, vcf_record):
        if 'GT' in vcf_record.FORMAT:
            called_alts = list(set(vcf_record.FORMAT['GT'].split('/')))
            if len(called_alts) > 1:
                vcf_record.set_format_key_value('MINOS_CHECK_GENOTYPE', 'HET')
            else:
                called_alt = int(called_alts[0])
                pass_per_alt = vcf_record.FORMAT['MINOS_CHECK_PASS'].split(',')
                vcf_record.set_format_key_value('MINOS_CHECK_GENOTYPE',pass_per_alt[called_alt])
        else:
            vcf_record.set_format_key_value('MINOS_CHECK_GENOTYPE', 'UNKNOWN_NO_GT')


    @classmethod
    def _check_if_sam_match_is_good(cls, sam_record, ref_seqs, flank_length, query_sequence=None):
        if sam_record.is_unmapped:
            return False

        # don't allow too many soft clipped bases
        if (sam_record.cigartuples[0] == 4 and sam_record.cigartuples[0][1] > 3) or (sam_record.cigartuples[-1][0] == 4 and sam_record.cigartuples[-1][1] > 3):
            return False

        if query_sequence is None:
            query_sequence = sam_record.query_sequence
        assert query_sequence is not None

        assert sam_record.reference_name in ref_seqs

        # if the query is short, which happens when the variant we
        # are checking is too near the start or end of the ref sequence
        if len(query_sequence) < 2 * flank_length + 1:
            # This is an edge case. We don't really know which part
            # of the query seq we're looking for, so guess
            length_diff = 2 * flank_length - len(query_sequence)

            if sam_record.query_alignment_start < 5:
                alt_seq_end = len(query_sequence) - flank_length - 1
                alt_seq_start = min(alt_seq_end, flank_length - length_diff)
            else:
                alt_seq_start = flank_length
                alt_seq_end = max(alt_seq_start, length_diff + len(query_sequence) - flank_length - 1)
        else:
            alt_seq_start = flank_length
            alt_seq_end = len(query_sequence) - flank_length - 1

        aligned_pairs = sam_record.get_aligned_pairs()
        wanted_aligned_pairs = []
        current_pos = 0

        i = 0
        while i < len(query_sequence):
            if aligned_pairs[i][0] is None:
                if alt_seq_start - 1 <= current_pos <= alt_seq_end + 1:
                    wanted_aligned_pairs.append(aligned_pairs[i])
            elif current_pos > alt_seq_end:
                break
            else:
                current_pos = aligned_pairs[i][0]
                if alt_seq_start - 1 <= current_pos <= alt_seq_end + 1:
                    wanted_aligned_pairs.append(aligned_pairs[i])

            i += 1

        assert len(wanted_aligned_pairs) > 0

        for pair in wanted_aligned_pairs:
            if None in pair or query_sequence[pair[0]] != ref_seqs[sam_record.reference_name][pair[1]]:
                return False

        return True


    @classmethod
    def _parse_sam_file_and_update_vcf_records_and_gather_stats(cls, infile, vcf_records, flank_length, ref_seqs):
        '''Input is SAM file made by _map_seqs_to_ref(), and corresponding dict
        of VCF records made by vcf_file_read.file_to_dict.
        Adds validation info to each VCF record. Returns a dict of stats that
        has counts up the values in MINOS_CHECK_GENOTYPE in the output file:
            total => total number of VCF records
            HET => number of records with a heterozygous call
            'gt_wrong' => number of records where genotype is incorrect
            'gt_correct' => number of records where genotype is correct'''
        stats = {x: 0 for x in ['total', '0', '1', 'HET', 'UNKNOWN_NO_GT']}
        gt_conf_hists = {'TP': {}, 'FP': {}}
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
            results = {x: [] for x in ['MINOS_CHECK_PASS']}

            for allele_sam_list in sam_records_by_allele:
                # Important! only the first hit actually has the sequence!
                # So give that to MappingBasedVerifier._check_if_sam_match_is_good()
                indexes_of_good_matches = [i for i in range(len(allele_sam_list)) if MappingBasedVerifier._check_if_sam_match_is_good(allele_sam_list[i], ref_seqs, flank_length, query_sequence=allele_sam_list[0].query_sequence)]

                if len(indexes_of_good_matches) > 0:
                    results['MINOS_CHECK_PASS'].append('1')
                else:
                    results['MINOS_CHECK_PASS'].append('0')

            vcf_record = vcf_records[ref_name][vcf_record_index]

            for key in sorted(results):
                vcf_record.set_format_key_value(key, ','.join(results[key]))

            MappingBasedVerifier._check_called_genotype(vcf_record)
            stats[vcf_record.FORMAT['MINOS_CHECK_GENOTYPE']] += 1

            if 'GT_CONF' in vcf_record.FORMAT and vcf_record.FORMAT['MINOS_CHECK_GENOTYPE'] in {'0', '1'}:
                tp_or_fp = {'1': 'TP', '0': 'FP'}[vcf_record.FORMAT['MINOS_CHECK_GENOTYPE']]
                gt_conf = int(vcf_record.FORMAT['GT_CONF'])
                gt_conf_hists[tp_or_fp][gt_conf] = gt_conf_hists[tp_or_fp].get(gt_conf, 0) + 1

        stats['gt_wrong'] = stats['0']
        del stats['0']
        stats['gt_correct'] = stats['1']
        del stats['1']
        return stats, gt_conf_hists


    @classmethod
    def _get_missing_vcf_records(cls, called_vcfs, expected_vcfs):
        missing_records = {}

        for ref_seq in expected_vcfs:
            if ref_seq not in called_vcfs:
                missing_records[ref_seq] = expected_vcfs[ref_seq]
                continue

            i = 0
            called_vcf_list = called_vcfs[ref_seq]
            missing_records[ref_seq] = []

            for expected_record in expected_vcfs[ref_seq]:
                while i < len(called_vcf_list) and called_vcf_list[i].ref_end_pos() < expected_record.POS:
                    i += 1

                if i == len(called_vcf_list):
                    missing_records[ref_seq].append(expected_record)
                else:
                    expected_alts = expected_record.called_alts_from_genotype()
                    called_alts = called_vcf_list[i].called_alts_from_genotype()
                    if expected_record.POS != called_vcf_list[i].POS or None in [expected_alts, called_alts] or expected_alts != called_alts:
                        missing_records[ref_seq].append(expected_record)

        return missing_records


    @classmethod
    def _get_total_length_of_expected_regions_called(cls, expected_regions, called_vcf_records):
        '''Returns total length of expected regions, and total length that was called'''
        total_length = 0
        called_length = 0

        for ref_seq, expected_list in expected_regions.items():
            if ref_seq in called_vcf_records:
                called_list = [pyfastaq.intervals.Interval(x.POS, x.ref_end_pos()) for x in called_vcf_records[ref_seq]]
                pyfastaq.intervals.merge_overlapping_in_list(called_list)
                total_length += pyfastaq.intervals.length_sum_from_list(expected_list)
                intersect_regions = pyfastaq.intervals.intersection(called_list, expected_list)
                called_length += pyfastaq.intervals.length_sum_from_list(intersect_regions)
            else:
                total_length += pyfastaq.intervals.length_sum_from_list(expected_list)

        return total_length, called_length


    def run(self):
        if self.filter_and_cluster_vcf:
            MappingBasedVerifier._filter_vcf_for_clustering(self.vcf_file_in, self.filtered_vcf)
            clusterer = vcf_clusterer.VcfClusterer([self.filtered_vcf], self.vcf_reference_file, self.clustered_vcf, merge_method='simple', max_distance_between_variants=self.flank_length)
            clusterer.run()

        vcf_header, vcf_records = vcf_file_read.vcf_file_to_dict(self.vcf_to_check, sort=True, remove_useless_start_nucleotides=True)
        sample_from_header = vcf_file_read.get_sample_name_from_vcf_header_lines(vcf_header)
        if sample_from_header is None:
            sample_from_header = 'sample'
        vcf_ref_seqs = {}
        pyfastaq.tasks.file_to_dict(self.vcf_reference_file, vcf_ref_seqs)
        verify_ref_seqs = {}
        pyfastaq.tasks.file_to_dict(self.verify_reference_file, verify_ref_seqs)

        MappingBasedVerifier._write_vars_plus_flanks_to_fasta(self.seqs_out, vcf_records, vcf_ref_seqs, self.flank_length)
        MappingBasedVerifier._map_seqs_to_ref(self.seqs_out, self.verify_reference_file, self.sam_file_out)
        os.unlink(self.seqs_out)
        stats, gt_conf_hists = MappingBasedVerifier._parse_sam_file_and_update_vcf_records_and_gather_stats(self.sam_file_out, vcf_records, self.flank_length, verify_ref_seqs)

        with open(self.vcf_file_out, 'w') as f:
            print(*vcf_header, sep='\n', file=f)
            for r in vcf_records:
                for v in vcf_records[r]:
                    print(v, file=f)

        # false negative stats, if possible
        stats['variant_regions_total'] = 'NA'
        stats['called_variant_regions'] = 'NA'

        if self.run_dnadiff:
            dnadiffer = dnadiff.Dnadiff(
                self.verify_reference_file,
                self.vcf_reference_file,
                self.dnadiff_outprefix,
            )
            dnadiffer.run()
            stats['variant_regions_total'], stats['called_variant_regions'] = MappingBasedVerifier._get_total_length_of_expected_regions_called(dnadiffer.all_variant_intervals, vcf_records)
            expected_variants = dnadiffer.variants
        elif self.expected_variants_vcf is not None:
            header, expected_variants = vcf_file_read.vcf_file_to_dict(self.expected_variants_vcf, sort=True, remove_useless_start_nucleotides=True)
        else:
            expected_variants = None

        if expected_variants is None:
            stats['false_negatives'] = 'NA'
        else:
            missed_vcf_records = MappingBasedVerifier._get_missing_vcf_records(vcf_records, expected_variants)
            stats['false_negatives'] = 0
            with open(self.vcf_false_negatives_file_out, 'w') as f:
                print('##fileformat=VCFv4.2', file=f)
                print('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_from_header, sep='\t', file=f)
                for vcf_list in missed_vcf_records.values():
                    stats['false_negatives'] += len(vcf_list)
                    print(*vcf_list, sep='\n', file=f)

        # write stats file
        with open(self.stats_out, 'w') as f:
            keys = ['total', 'gt_correct', 'gt_wrong', 'HET', 'UNKNOWN_NO_GT', 'variant_regions_total', 'called_variant_regions', 'false_negatives']
            print(*keys, sep='\t', file=f)
            print(*[stats[x] for x in keys], sep='\t', file=f)


        # write GT_CONG histogram files
        for key, filename in self.gt_conf_hists_filenames.items():
            with open(filename, 'w') as f:
                print('GT_CONF\tCount', file=f)
                for gt_conf, count in sorted(gt_conf_hists[key].items()):
                    print(gt_conf, count, sep='\t', file=f)

