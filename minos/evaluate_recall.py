import logging
import os
import math
import glob

import pyfastaq
import pysam
import pandas as pd

from Bio import pairwise2, SeqIO

from cluster_vcf_records import vcf_clusterer, vcf_file_read

from minos import dependencies, dnadiff, plots, utils

class Error (Exception): pass

class EvaluateRecall:
    '''tto write'''
    def __init__(self, truth_vcf_file, truth_vcf_ref, query_vcf_file, query_vcf_ref, outprefix, flank_length=31, merge_length=None, filter_and_cluster_vcf=True, discard_ref_calls=True, allow_flank_mismatches=True, exclude_regions_bed_file=None):
        self.truth_vcf_file = os.path.abspath(truth_vcf_file)
        self.truth_vcf_ref = os.path.abspath(truth_vcf_ref)
        self.query_vcf_file = os.path.abspath(query_vcf_file)
        self.query_vcf_ref = os.path.abspath(query_vcf_ref)
        self.sam_file_out = os.path.abspath(outprefix + '.sam')
        self.seqs_out_truth = os.path.abspath(outprefix + '.truth.fa')
        self.filtered_query_vcf = os.path.abspath(outprefix + '.query.filter.vcf')
        self.filtered_truth_vcf = os.path.abspath(outprefix + '.truth.filter.vcf')
        self.clustered_vcf_query = os.path.abspath(outprefix + '.query.filter.cluster.vcf')
        self.clustered_vcf_truth = os.path.abspath(outprefix + '.truth.filter.cluster.vcf')
        self.seqs_out_query = os.path.abspath(outprefix + '.query.fa')
        self.sam_summary = os.path.abspath(outprefix + '.summary.tsv')
        self.stats_out = os.path.abspath(outprefix + '.stats.tsv')
        self.gt_conf_hist_out = os.path.abspath(outprefix + '.gt_conf_hist.tsv')

        self.flank_length = flank_length
        self.merge_length = flank_length if merge_length is None else merge_length
        self.filter_and_cluster_vcf = filter_and_cluster_vcf
        self.discard_ref_calls = discard_ref_calls
        self.allow_flank_mismatches = allow_flank_mismatches
        print("flank_length = " + str(self.flank_length) + " and merge_length = " + str(self.merge_length))

        if self.filter_and_cluster_vcf:
            self.vcf_to_check_truth = self.clustered_vcf_truth
            self.vcf_to_check_query = self.clustered_vcf_query
        else:
            self.vcf_to_check_truth = self.truth_vcf_file
            self.vcf_to_check_query = self.query_vcf_file

        self.exclude_regions = EvaluateRecall._load_exclude_regions_bed_file(exclude_regions_bed_file)

    @classmethod
    def _load_exclude_regions_bed_file(cls, infile):
        regions = {}
        if infile is not None:
            with open(infile) as f:
                for line in f:
                    fields = line.rstrip().split('\t')
                    if fields[0] not in regions:
                        regions[fields[0]] = []
                    start = int(fields[1]) - 1
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
                    if vcf_record.FILTER == 'MISMAPPED_UNPLACEABLE':
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
    def _write_vars_plus_flanks_to_fasta(cls, outfile, vcf_records, ref_seqs, flank_length, ref_only=False):
        '''Given a dict of vcf records made by vcf_file_read.vcf_file_to_dict(),
        and its correcsponding file of reference sequences, writes a new fasta file
        of each ref seq and inferred variant sequence plus flank_length nucleotides added to
        its start and end. Calls each sequence:
            ref_name.start_position.vcf_list_index.allele_number
        where allele_numbers in same order as VCF, with ref seq = allele 0.'''
        prev_ref_name = None
        prev_ref_pos = None
        j = 0
        with open(outfile, 'w') as f:
            for ref_name in sorted(vcf_records):
                for i, vcf_record in enumerate(vcf_records[ref_name]):
                    start_position, alleles = vcf_record.inferred_var_seqs_plus_flanks(ref_seqs[ref_name], flank_length)

                    for allele_index, allele_seq in enumerate(alleles):
                        if not ref_only or allele_index == 0:
                            seq_name = '.'.join([ref_name, str(start_position + 1), str(j), str(i), str(allele_index)])
                            allele_seq = allele_seq.replace('.','')
                            print('>' + seq_name, allele_seq, sep='\n', file=f)
                    if prev_ref_name == ref_name and prev_ref_pos == start_position:
                        j += 1
                    else:
                        j = 0
                    prev_ref_name = ref_name
                    prev_ref_pos = start_position



    @classmethod
    def _map_seqs_to_seqs(cls, seqs_file_ref, seqs_file_query, outfile):
        '''Map seqs_file to ref_file using BWA MEM.
        Output is SAM file written to outfile'''
        bwa_binary = dependencies.find_binary('bwa')
        command = ' '.join([
            bwa_binary, 'index',
            seqs_file_ref,
        ])
        utils.syscall(command)

        command = ' '.join([
            bwa_binary, 'aln',
            seqs_file_ref,
            seqs_file_query,
            '>', outfile + ".tmp",
        ])
        utils.syscall(command)

        command = ' '.join([
            bwa_binary, 'samse',
            seqs_file_ref,
            outfile + ".tmp",
            seqs_file_query,
            '>', outfile,
        ])
        utils.syscall(command)
        #os.unlink(outfile + ".tmp")

    @classmethod
    def _check_if_sam_match_is_good(cls, sam_record, flank_length, query_sequence=None, allow_mismatches=True):
        if sam_record.is_unmapped:
            try:
                #NB some mapped things at the edge of reference sequence will get unmapped flag, so check if this is one
                nm = sam_record.get_tag('NM')
                all_mapped = len(sam_record.cigartuples) == 1 and sam_record.cigartuples[0][0] == 0
            except:
                #real unmapped, no cigar or tag
                print("unmapped and no nm or cigar tags")
                return False

        try:
            nm = sam_record.get_tag('NM')
        except:
            raise Error('No NM tag found in sam record:' + str(sam_record))

        if not allow_mismatches:
            all_mapped = len(sam_record.cigartuples) == 1 and sam_record.cigartuples[0][0] == 0
            print("no mismatches allowed: ", all_mapped, nm==0)
            return all_mapped and nm == 0

        # don't allow too many soft clipped bases
        if (sam_record.cigartuples[0][0] == 4 and sam_record.cigartuples[0][1] > 3) or (sam_record.cigartuples[-1][0] == 4 and sam_record.cigartuples[-1][1] > 3):
            print("too many soft clipped bases")
            return False

        if nm==0:
            print("no mismatches")
            return True

        if query_sequence is None:
            query_sequence = sam_record.query_sequence
        assert query_sequence is not None

        #assert sam_record.reference_name in ref_seqs

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
        print(aligned_pairs)
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

        #assert len(wanted_aligned_pairs) > 0

        #for pair in wanted_aligned_pairs:
        #    if None in pair or query_sequence[pair[0]] != ref_seqs[sam_record.reference_name][pair[1]]:
        #        return False

        return True


    @classmethod
    def _index_vcf(cls, vcffile):
        '''Index VCF file'''
        bgzip_binary = dependencies.find_binary('bgzip')
        command = ' '.join([
            bgzip_binary,
            '-c',
            vcffile,
            '>',
            vcffile + ".gz",
        ])
        utils.syscall(command)

        tabix_binary = dependencies.find_binary('tabix')
        command = ' '.join([
            tabix_binary,
            '-p',
            'vcf',
            vcffile + ".gz",
        ])
        utils.syscall(command)

    @classmethod
    def _parse_sam_file_and_vcf(cls, samfile, query_vcf_file, flank_length, allow_mismatches, exclude_regions=None):
        if  exclude_regions is None:
            exclude_regions = {}

        print("start _parse_sam_file_and_vcf")
        found = []
        gt_conf = []
        allele = []
        samfile_handle = pysam.AlignmentFile(samfile, "r")
        sam_previous_record_name = None
        for sam_record in samfile_handle.fetch(until_eof=True):
            print(str(sam_record))
            if sam_record.query_name == sam_previous_record_name:
                continue
            sam_previous_record_name = sam_record.query_name
            found_conf = False
            found_allele = False

            # see if excluded region in bed file
            ref, start, ref_num, var_num, allele_num = sam_record.query_name.rsplit('.', maxsplit=5)
            start = int(start) + flank_length
            exclude = False
            for ref_name in exclude_regions.keys():
                end = int(start) + 1
                interval = pyfastaq.intervals.Interval(start, end)
                exclude = EvaluateRecall._interval_intersects_an_interval_in_list(interval,
                                                                                  exclude_regions[ref_name])
            if exclude:
                found.append('E')
                gt_conf.append(0)
                allele.append('0')
                continue

            good_match = EvaluateRecall._check_if_sam_match_is_good(sam_record,
                                                                    flank_length,
                                                                    query_sequence=sam_record.query_sequence,
                                                                    allow_mismatches=allow_mismatches)
            alignment_start = str(sam_record).split("\t")[3]
            if good_match:
                ref_name, expected_start, vcf_pos_index, vcf_record_index, allele_index = sam_record.reference_name.rsplit('.', maxsplit=4)

                vcf_reader = pysam.VariantFile(query_vcf_file)
                for i, vcf_record in enumerate(vcf_reader.fetch(ref_name, int(expected_start) + int(alignment_start) + flank_length - 2, int(expected_start) + int(alignment_start) + flank_length)):
                    print(vcf_record)
                    if i == int(vcf_pos_index):
                        sample_name = vcf_record.samples.keys()[0]
                        if 'GT' in vcf_record.format.keys() and len(set(vcf_record.samples[sample_name]['GT'])) == 1:
                            if int(allele_index) == vcf_record.samples[sample_name]['GT'][0]:
                                found.append('1')
                                allele.append(str(allele_index))
                                found_allele = True
                                print(vcf_record.format.keys())
                                if 'GT_CONF' in vcf_record.format.keys():
                                    gt_conf.append(int(float(vcf_record.samples[sample_name]['GT_CONF'])))
                                    found_conf = True
            if not found_allele:
                found.append('0')
                allele.append('0')
            if not found_conf:
                gt_conf.append(0)
        assert len(found) == len(gt_conf)
        assert len(found) == len(allele)
        print(found)
        print(gt_conf)
        print(allele)
        return found, gt_conf, allele

    @classmethod
    def _parse_sam_files(cls, truth_vcf_file, samfile, query_vcf_file, outfile, flank_length, allow_mismatches=True, exclude_regions=None):
        '''Input is the original dnadiff snps file of sites we are searching for
        and 2 SAM files made by _map_seqs_to_seqs(), which show mappings of snp sites
        from from the dnadiff snps file to the vcf (i.e. searches if VCF contains an record
        with the appropriate sequence.
        Creates a tsv detailing whether the snp difference could be detected and at what
        GT_CONF threshold.
        '''
        header_lines, vcf_records = vcf_file_read.vcf_file_to_dict(truth_vcf_file, sort=True,
                                                                   homozygous_only=False,
                                                                   remove_asterisk_alts=True,
                                                                   remove_useless_start_nucleotides=True)
        id = []
        ref = []
        alt = []
        for ref_name in vcf_records:
            for record in vcf_records[ref_name]:
                id.append(record.POS)
                ref.append(record.REF)
                alt.append(record.ALT[0])
        query_found, query_conf, query_allele = EvaluateRecall._parse_sam_file_and_vcf(samfile, query_vcf_file,
                                                                                       flank_length, allow_mismatches,
                                                                                       exclude_regions)
        assert len(id) == len(query_found)
        out_df = pd.DataFrame({'id': id,
                               'ref': ref,
                               'alt': alt,
                               'query_found': query_found,
                               'query_conf': query_conf,
                               'query_allele': query_allele})
        out_df.to_csv(outfile, sep='\t')

    @classmethod
    def _gather_stats(cls, tsv_file):
        stats = {x: 0 for x in ['total', 'found_vars', 'missed_vars', 'excluded_vars']}
        gt_conf_hist = {}

        snps = pd.read_table(tsv_file, index_col=0)
        for line in snps.itertuples():
            stats['total'] += 1
            print(line)
            if (line[4] == 'E'):
                stats['excluded_vars'] += 1
            elif (line[4] == 1 or line[4] == '1'):
                print("found")
                stats['found_vars'] += 1
                gt_confs = [i for i in {line[5]} if not math.isnan(i)]
                gt_conf = None
                if len(gt_confs) > 0:
                    gt_conf = max(gt_confs)
                gt_conf_hist[gt_conf] = gt_conf_hist.get(gt_conf, 0) + 1
            else:
                stats['missed_vars'] += 1
        return stats, gt_conf_hist

    def run(self):
        # Cluster together variants in each vcf
        if self.filter_and_cluster_vcf:
            EvaluateRecall._filter_vcf_for_clustering(self.truth_vcf_file, self.filtered_truth_vcf, self.discard_ref_calls)
            EvaluateRecall._filter_vcf_for_clustering(self.query_vcf_file, self.filtered_query_vcf, self.discard_ref_calls)
            if self.discard_ref_calls:
                clusterer_query = vcf_clusterer.VcfClusterer([self.filtered_query_vcf], self.query_vcf_ref, self.clustered_vcf_query, merge_method='simple', max_distance_between_variants=self.merge_length)
                clusterer_truth = vcf_clusterer.VcfClusterer([self.filtered_truth_vcf], self.truth_vcf_ref, self.clustered_vcf_truth, merge_method='simple', max_distance_between_variants=self.merge_length)
            else:
                clusterer_query = vcf_clusterer.VcfClusterer([self.filtered_query_vcf], self.query_vcf_ref, self.clustered_vcf_query, merge_method='gt_aware', max_distance_between_variants=self.merge_length)
                clusterer_truth = vcf_clusterer.VcfClusterer([self.filtered_truth_vcf], self.truth_vcf_ref, self.clustered_vcf_truth, merge_method='gt_aware', max_distance_between_variants=self.merge_length)
            clusterer_query.run()
            clusterer_truth.run()

        vcf_header, vcf_records_truth = vcf_file_read.vcf_file_to_dict(self.vcf_to_check_truth, sort=True, remove_useless_start_nucleotides=True)
        vcf_header, vcf_records_query = vcf_file_read.vcf_file_to_dict(self.vcf_to_check_query, sort=True, remove_useless_start_nucleotides=True)
        sample_from_header = vcf_file_read.get_sample_name_from_vcf_header_lines(vcf_header)
        if sample_from_header is None:
            sample_from_header = 'sample'
        truth_vcf_ref_seqs = {}
        pyfastaq.tasks.file_to_dict(self.truth_vcf_ref, truth_vcf_ref_seqs)
        query_vcf_ref_seqs = {}
        pyfastaq.tasks.file_to_dict(self.query_vcf_ref, query_vcf_ref_seqs)

        EvaluateRecall._write_vars_plus_flanks_to_fasta(self.seqs_out_truth, vcf_records_truth, truth_vcf_ref_seqs, self.flank_length, ref_only=True)
        EvaluateRecall._write_vars_plus_flanks_to_fasta(self.seqs_out_query, vcf_records_query, query_vcf_ref_seqs, self.flank_length)
        EvaluateRecall._map_seqs_to_seqs(self.seqs_out_query, self.seqs_out_truth, self.sam_file_out)
        #for f in glob.glob(self.seqs_out_truth + '*'):
            #os.unlink(f)
        #for f in glob.glob(self.seqs_out_query + '*'):
            #os.unlink(f)

        EvaluateRecall._index_vcf(self.vcf_to_check_query)
        self.vcf_to_check_query = self.vcf_to_check_query + ".gz"
        EvaluateRecall._parse_sam_files(self.truth_vcf_file, self.sam_file_out, self.vcf_to_check_query,
                                        self.sam_summary, self.flank_length,
                                        allow_mismatches=self.allow_flank_mismatches,
                                        exclude_regions=self.exclude_regions)
        stats, gt_conf_hist = EvaluateRecall._gather_stats(self.sam_summary)
        #os.unlink(self.seqs_out_truth)
        #os.unlink(self.seqs_out_truth)
        #for f in glob.glob(self.vcf_to_check_truth + '*'):
        #    os.unlink(f)
        #for f in glob.glob(self.vcf_to_check_query + '*'):
        #    os.unlink(f)

        # write stats file
        with open(self.stats_out, 'w') as f:
            keys = stats.keys()
            print(*keys, sep='\t', file=f)
            print(*[stats[x] for x in keys], sep='\t', file=f)


        # write GT_CONF histogram files
        with open(self.gt_conf_hist_out, 'w') as f:
            print('GT_CONF\tCount', file=f)
            for gt_conf, count in sorted(gt_conf_hist.items()):
                print(gt_conf, count, sep='\t', file=f)