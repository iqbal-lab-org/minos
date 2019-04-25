import logging
import os

import pyfastaq
import pysam

from Bio import pairwise2

from cluster_vcf_records import vcf_clusterer, vcf_file_read

from minos import dependencies, dnadiff, plots, utils

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
    def __init__(self, vcf_file_in, vcf_reference_file, verify_reference_file, outprefix, flank_length=31, merge_length=None, expected_variants_vcf=None, run_dnadiff=True, filter_and_cluster_vcf=True, discard_ref_calls=True, allow_flank_mismatches=True, exclude_regions_bed_file=None, max_soft_clipped=3):
        self.vcf_file_in = os.path.abspath(vcf_file_in)
        self.vcf_reference_file = os.path.abspath(vcf_reference_file)
        self.verify_reference_file = os.path.abspath(verify_reference_file)
        self.vcf_file_out = os.path.abspath(outprefix + '.vcf')
        self.vcf_file_plots_out = os.path.abspath(outprefix + '.vcf.plots')
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
        self.merge_length = flank_length if merge_length is None else merge_length
        self.expected_variants_vcf = expected_variants_vcf
        self.run_dnadiff = run_dnadiff
        if self.run_dnadiff and self.expected_variants_vcf is not None:
            raise Error('Error! Incompatible options. expected_variants_vcf file provided, and run_dnadiff is True')
        self.dnadiff_outprefix = os.path.abspath(outprefix + '.dnadiff')
        self.vcf_false_negatives_file_out = os.path.abspath(outprefix + '.false_negatives.vcf')
        self.filter_and_cluster_vcf = filter_and_cluster_vcf
        self.discard_ref_calls = discard_ref_calls
        self.allow_flank_mismatches = allow_flank_mismatches

        if self.filter_and_cluster_vcf:
            self.vcf_to_check = self.clustered_vcf
        else:
            self.vcf_to_check = self.vcf_file_in

        self.exclude_regions = MappingBasedVerifier._load_exclude_regions_bed_file(exclude_regions_bed_file)
        self.max_soft_clipped = max_soft_clipped


    @classmethod
    def _needleman_wunsch(cls, seq1, seq2, match=1, mismatch=-3, gap_open=-11, gap_extend=-4):
        '''Returns global alignment strings from NM alignment of the
        two sequences. Dashes for gaps'''
        alignments = pairwise2.align.globalms(seq1, seq2, match, mismatch, gap_open, gap_extend)
        assert len(alignments[0][0]) == len(alignments[0][1])
        return alignments[0][0], alignments[0][1]


    @classmethod
    def _edit_distance_from_alignment_strings(cls, str1, str2):
        '''Input should be seqs output by _needleman_wunsch().
        Returns the edit distance between the sequences'''
        assert len(str1) == len(str2)
        edit_distance = 0
        in_gap = False

        for i, char1 in enumerate(str1):
            if char1 == '-' or str2[i] == '-':
                if not in_gap:
                    in_gap = True
                    edit_distance += 1
            else:
                in_gap = False

                if char1 != str2[i]:
                    edit_distance += 1

        return edit_distance


    @classmethod
    def _edit_distance_between_seqs(cls, seq1, seq2):
        '''Input is two strings). They are globally aligned
        and the edit distance is returned. An indel of any length
        is counted as one edit'''
        aln1, aln2 = MappingBasedVerifier._needleman_wunsch(seq1, seq2)
        return MappingBasedVerifier._edit_distance_from_alignment_strings(aln1, aln2)


    @classmethod
    def _add_edit_distances_to_vcf_record(cls, record):
        '''Adds EDIT_DIST => list of numbers, one for each allele, to the INFO column.
        Also adds GT_EDIT_DIST => float, which is the edit distance between the genotyped allele
        and the ref. If het call, takes the mean edit distance of the two calls'''
        edit_distances = [MappingBasedVerifier._edit_distance_between_seqs(record.REF, x) for x in record.ALT]
        record.INFO['EDIT_DIST'] = ','.join([str(x) for x in edit_distances])
        if 'GT' in record.FORMAT:
            genotypes = set([int(x) for x in record.FORMAT['GT'].split('/')])
            geno_edit_distances = []
            for genotype in genotypes:
                if genotype == 0:
                    geno_edit_distances.append(0)
                else:
                    geno_edit_distances.append(edit_distances[genotype - 1])
            gt_edit_dist = round(sum(geno_edit_distances) / len(geno_edit_distances), 1)
            record.set_format_key_value('GT_EDIT_DIST', str(gt_edit_dist))


    @classmethod
    def _load_exclude_regions_bed_file(cls, infile):
        regions = {}
        if infile is not None:
            with open(infile) as f:
                for line in f:
                    fields = line.rstrip().split('\t')
                    if fields[0] not in regions:
                        regions[fields[0]] = []
                    start = int(fields[1])
                    end = int(fields[2]) - 1
                    regions[fields[0]].append(pyfastaq.intervals.Interval(start, end))

            for ref_name in regions:
                pyfastaq.intervals.merge_overlapping_in_list(regions[ref_name])

        return regions


    @classmethod
    def _interval_intersects_an_interval_in_list(cls, interval, interval_list):
        # This could be faster by doing something like a binary search.
        # But we're looking for points in intervals, so fiddly to implement.
        # Not expecting a log interval list, so just do a simple check
        # from start to end for now
        i = 0
        while i < len(interval_list) and interval.start > interval_list[i].end:
            i += 1

        return i < len(interval_list) and interval.intersects(interval_list[i])


    @classmethod
    def _filter_vcf_for_clustering(cls, infile, outfile, discard_ref_calls=True):
        header_lines, vcf_records = vcf_file_read.vcf_file_to_dict(infile, sort=True, homozygous_only=False, remove_asterisk_alts=True, remove_useless_start_nucleotides=True)

        with open(outfile, 'w') as f:
            print(*header_lines, sep='\n', file=f)
            for ref_name in vcf_records:
                for vcf_record in vcf_records[ref_name]:
                    if 'MISMAPPED_UNPLACEABLE' in vcf_record.FILTER:
                        continue
                    if vcf_record.FORMAT is None or 'GT' not in vcf_record.FORMAT:
                        logging.warning('No GT in vcf record:' + str(vcf_record))
                        continue
                    if vcf_record.REF in [".", ""]:
                        continue

                    genotype = vcf_record.FORMAT['GT']
                    genotypes = genotype.split('/')
                    called_alleles = set(genotypes)
                    if len(called_alleles) != 1 or (discard_ref_calls and called_alleles == {'0'}) or '.' in called_alleles:
                        continue

                    if len(vcf_record.ALT) > 1:
                        if called_alleles != {'0'}:
                            vcf_record.set_format_key_value('GT', '1/1')
                            try:
                                vcf_record.ALT = [vcf_record.ALT[int(genotypes[0]) - 1]]
                            except:
                                raise Error('BAD VCf line:' + str(vcf_record))
                        else:
                            vcf_record.set_format_key_value('GT', '0/0')
                            vcf_record.ALT = [vcf_record.ALT[0]]
                    if vcf_record.ALT[0] in [".",""]:
                        continue

                    if vcf_record.FORMAT['GT'] == '0':
                        vcf_record.FORMAT['GT'] = '0/0'
                    elif vcf_record.FORMAT['GT'] == '1':
                        vcf_record.FORMAT['GT'] = '1/1'

                    if 'GL' in vcf_record.FORMAT.keys() and 'GT_CONF' not in vcf_record.FORMAT.keys():
                        likelihoods = vcf_record.FORMAT['GL'].split(',')
                        assert(len(likelihoods) > 2)
                        if called_alleles == {'0'}:
                            vcf_record.set_format_key_value('GT_CONF',str(float(likelihoods[0]) - float(likelihoods[1])))
                        else:
                            vcf_record.set_format_key_value('GT_CONF', str(float(likelihoods[int(genotypes[0])]) - float(likelihoods[0])))
                    if 'SupportFraction' in vcf_record.INFO.keys() and 'GT_CONF' not in vcf_record.FORMAT.keys():
                        vcf_record.set_format_key_value('GT_CONF',
                                                        str(float(vcf_record.INFO['SupportFraction'])*100))
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
                pass_per_alt = vcf_record.FORMAT['MINOS_CHECK_ALLELES'].split(',')
                if pass_per_alt[called_alt] == 'Pass':
                    check_geno = '1'
                elif pass_per_alt[called_alt] == 'Exclude':
                    check_geno = 'Exclude'
                else:
                    check_geno = '0'
                vcf_record.set_format_key_value('MINOS_CHECK_GENOTYPE', check_geno)
        else:
            vcf_record.set_format_key_value('MINOS_CHECK_GENOTYPE', 'UNKNOWN_NO_GT')


    @classmethod
    def _check_if_sam_match_is_good(cls, sam_record, ref_seqs, flank_length, query_sequence=None, allow_mismatches=True, max_soft_clipped=3):
        logging.debug(f'Checking SAM: {sam_record}')

        if sam_record.is_unmapped:
            return 'Unmapped'

        if not allow_mismatches:
            try:
                nm = sam_record.get_tag('NM')
            except:
                raise Error('No NM tag found in sam record:' + str(sam_record))

            all_mapped = len(sam_record.cigartuples) == 1 and sam_record.cigartuples[0][0] == 0
            if all_mapped and nm == 0:
                logging.debug('SAM record passed no mismatches allowed check')
                return 'Good'
            else:
                logging.debug('SAM record failed no mismatches allowed check')
                return 'Bad_mismatches'

        # don't allow too many soft clipped bases
        if (sam_record.cigartuples[0][0] == 4 and sam_record.cigartuples[0][1] > max_soft_clipped) \
                or (sam_record.cigartuples[-1][0] == 4 and sam_record.cigartuples[-1][1] > max_soft_clipped):
            logging.debug('SAM record failed soft clipping check')
            return 'Bad_soft_clipped'

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
        logging.debug(f'aligned_pairs: {aligned_pairs}')
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

        logging.debug(f'wanted_aligned_pairs: {wanted_aligned_pairs}')
        assert len(wanted_aligned_pairs) > 0

        for pair in wanted_aligned_pairs:
            if None in pair or query_sequence[pair[0]] != ref_seqs[sam_record.reference_name][pair[1]]:
                logging.debug('SAM record failed because mismatch in allele sequence plus 1bp either side')
                return 'Bad_allele_mismatch'

        logging.debug('SAM record passed all checks')
        return 'Good'


    @classmethod
    def _parse_sam_file_and_update_vcf_records_and_gather_stats(cls, infile, vcf_records, flank_length, ref_seqs, allow_mismatches=True, exclude_regions=None, max_soft_clipped=3):
        '''Input is SAM file made by _map_seqs_to_ref(), and corresponding dict
        of VCF records made by vcf_file_read.file_to_dict.
        Adds validation info to each VCF record. Returns a dict of stats that
        has counts up the values in MINOS_CHECK_GENOTYPE in the output file:
            total => total number of VCF records
            HET => number of records with a heterozygous call
            'gt_wrong' => number of records where genotype is incorrect
            'gt_correct' => number of records where genotype is correct
            'gt_excluded' => number of records excluded because matched to a position in exclude_regions'''
        if  exclude_regions is None:
            exclude_regions = {}

        stats = {x: 0 for x in ['total', '0', '1', 'HET', 'UNKNOWN_NO_GT', 'Exclude', 'tp_edit_dist', 'fp_edit_dist']}
        gt_conf_hists = {'TP': {}, 'FP': {}}
        sreader = sam_reader(infile)

        for sam_list in sreader:
            stats['total'] += 1

            # Split the sam_list (Which has all records for all alleles from one
            # VCf line) into a list per allele.
            sam_records_by_allele = []
            current_allele_index = None
            vcf_probe_info_tuples = set()
            for sam_record in sam_list:
                vcfref_name, expected_start, vcf_record_index, allele_index = sam_record.query_name.rsplit('.', maxsplit=3)
                vcf_probe_info_tuples.add((vcfref_name, expected_start, vcf_record_index))
                if current_allele_index is None or current_allele_index != allele_index:
                    sam_records_by_allele.append([sam_record])
                else:
                    sam_records_by_allele[-1].append(sam_record)
                current_allele_index = allele_index

            assert len(vcf_probe_info_tuples) == 1
            vcfref_name, expected_start, vcf_record_index = vcf_probe_info_tuples.pop()
            expected_start = int(expected_start) - 1
            vcf_record_index = int(vcf_record_index)
            results = {x: [] for x in ['MINOS_CHECK_ALLELES']}

            for allele_sam_list in sam_records_by_allele:
                # Important! only the first hit actually has the sequence!
                # So give that to MappingBasedVerifier._check_if_sam_match_is_good()
                match_result_types = set()
                for i in range(len(allele_sam_list)):
                    if allele_sam_list[i].is_unmapped:
                        match_result_types.add('Unmapped')
                    else:
                        exclude = False
                        truth_name = allele_sam_list[i].reference_name
                        if truth_name in exclude_regions:
                            start = allele_sam_list[i].reference_start
                            end = allele_sam_list[i].reference_end - 1
                            interval = pyfastaq.intervals.Interval(start, end)
                            exclude = MappingBasedVerifier._interval_intersects_an_interval_in_list(interval, exclude_regions[truth_name])

                        if exclude:
                            match_result_types.add('Exclude')
                        else:
                            match_type = MappingBasedVerifier._check_if_sam_match_is_good(allele_sam_list[i], ref_seqs, flank_length, query_sequence=allele_sam_list[0].query_sequence, allow_mismatches=allow_mismatches, max_soft_clipped=max_soft_clipped)
                            match_result_types.add(match_type)
                            #if good_match:
                            #    match_result_types.add('1')
                            #else:
                            #    match_result_types.add('0')


                #indexes_of_good_matches = [i for i in range(len(allele_sam_list)) if MappingBasedVerifier._check_if_sam_match_is_good(allele_sam_list[i], ref_seqs, flank_length, query_sequence=allele_sam_list[0].query_sequence, allow_mismatches=allow_mismatches)]

                #if len(indexes_of_good_matches) > 0:
                #    results['MINOS_CHECK_ALLELES'].append('1')
                #else:
                #    results['MINOS_CHECK_ALLELES'].append('0')

                logging.debug(f'match_result_types: {match_result_types}')
                if 'Good' in match_result_types:
                    results['MINOS_CHECK_ALLELES'].append('Pass')
                elif 'Exclude' in match_result_types:
                    results['MINOS_CHECK_ALLELES'].append('Exclude')
                else:
                    results['MINOS_CHECK_ALLELES'].append('Fail-' + '-'.join(sorted(list(match_result_types))))

            vcf_record = vcf_records[vcfref_name][vcf_record_index]

            for key in sorted(results):
                vcf_record.set_format_key_value(key, ','.join(results[key]))

            MappingBasedVerifier._check_called_genotype(vcf_record)
            stats[vcf_record.FORMAT['MINOS_CHECK_GENOTYPE']] += 1
            MappingBasedVerifier._add_edit_distances_to_vcf_record(vcf_record)

            if 'GT_CONF' in vcf_record.FORMAT and vcf_record.FORMAT['MINOS_CHECK_GENOTYPE'] in {'0', '1'}:
                tp_or_fp = {'1': 'TP', '0': 'FP'}[vcf_record.FORMAT['MINOS_CHECK_GENOTYPE']]
                gt_conf = int(float(vcf_record.FORMAT['GT_CONF']))
                gt_conf_hists[tp_or_fp][gt_conf] = gt_conf_hists[tp_or_fp].get(gt_conf, 0) + 1

            # if it's a homozygous call, update the edit distance totals
            if 'GT' in vcf_record.FORMAT and len(set(vcf_record.FORMAT['GT'].split('/'))) == 1:
                if vcf_record.FORMAT['MINOS_CHECK_GENOTYPE'] == '1':
                    stats['tp_edit_dist'] += int(float(vcf_record.FORMAT['GT_EDIT_DIST']))
                elif vcf_record.FORMAT['MINOS_CHECK_GENOTYPE'] == '0':
                    stats['fp_edit_dist'] += int(float(vcf_record.FORMAT['GT_EDIT_DIST']))


        stats['gt_wrong'] = stats['0']
        del stats['0']
        stats['gt_correct'] = stats['1']
        del stats['1']
        stats['gt_excluded'] = stats['Exclude']
        del stats['Exclude']
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
                elif called_vcf_list[i].POS > expected_record.POS or expected_record.ref_end_pos() > called_vcf_list[i].ref_end_pos():
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
            MappingBasedVerifier._filter_vcf_for_clustering(self.vcf_file_in, self.filtered_vcf, discard_ref_calls=self.discard_ref_calls)
            if self.discard_ref_calls:
                clusterer = vcf_clusterer.VcfClusterer([self.filtered_vcf], self.vcf_reference_file, self.clustered_vcf, merge_method='simple', max_distance_between_variants=self.merge_length)
            else:
                clusterer = vcf_clusterer.VcfClusterer([self.filtered_vcf], self.vcf_reference_file, self.clustered_vcf, merge_method='gt_aware', max_distance_between_variants=self.merge_length)
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
        stats, gt_conf_hists = MappingBasedVerifier._parse_sam_file_and_update_vcf_records_and_gather_stats(self.sam_file_out, vcf_records, self.flank_length, verify_ref_seqs, allow_mismatches=self.allow_flank_mismatches, exclude_regions=self.exclude_regions, max_soft_clipped=self.max_soft_clipped)

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
                    if len(vcf_list) > 0:
                        print(*vcf_list, sep='\n', file=f)

        # write stats file
        with open(self.stats_out, 'w') as f:
            keys = ['total', 'gt_correct', 'gt_wrong', 'gt_excluded', 'HET', 'tp_edit_dist', 'fp_edit_dist', 'UNKNOWN_NO_GT', 'variant_regions_total', 'called_variant_regions', 'false_negatives']
            print(*keys, sep='\t', file=f)
            print(*[stats[x] for x in keys], sep='\t', file=f)


        # write GT_CONG histogram files
        for key, filename in self.gt_conf_hists_filenames.items():
            with open(filename, 'w') as f:
                print('GT_CONF\tCount', file=f)
                for gt_conf, count in sorted(gt_conf_hists[key].items()):
                    print(gt_conf, count, sep='\t', file=f)

        plots.plots_from_minos_vcf(self.vcf_file_out, self.vcf_file_plots_out)

