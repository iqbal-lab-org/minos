from collections import namedtuple
import logging
import os
import pickle
import json

import cluster_vcf_records

from minos import gramtools

class Error (Exception): pass


split_file_attributes = [
    'filename',
    'file_number',
    'chrom',
    'chrom_start',
    'chrom_end',
    'file_start_index',
    'file_end_index',
    'use_start_index',
    'use_end_index',
    'gramtools_build_dir',
]
SplitFile = namedtuple('SplitFile', split_file_attributes)


class VcfChunker:
    def __init__(self, outdir, vcf_infile=None, ref_fasta=None, variants_per_split=None, max_read_length=200, total_splits=100, flank_length=200, gramtools_kmer_size=15):
        self.outdir = os.path.abspath(outdir)
        self.metadata_pickle = os.path.join(self.outdir, 'data.pickle')

        if os.path.exists(self.outdir):
            self._load_existing_data()
        else:
            logging.info('Setting up new chunked VCF directory ' + self.outdir)
            assert vcf_infile is not None
            assert ref_fasta is not None
            self.vcf_infile = os.path.abspath(vcf_infile)
            self.ref_fasta = os.path.abspath(ref_fasta)
            self.variants_per_split = variants_per_split
            self.total_splits = total_splits
            self.flank_length = flank_length
            self.gramtools_kmer_size = gramtools_kmer_size
            self.max_read_length = max_read_length

            if not os.path.exists(self.vcf_infile):
                raise Error('VCF file not found: ' + self.vcf_infile)
            if not os.path.exists(self.ref_fasta):
                raise Error('Reference FASTA file not found: ' + self.ref_fasta)

            try:
                os.mkdir(self.outdir)
            except:
                raise Error('Error mkdir ' + self.outdir)

            self.vcf_split_files = {} # ref name -> list of SplitFile


    def _save_metadata(self):
        metadata = {
            'vcf_infile': self.vcf_infile,
            'ref_fasta': self.ref_fasta,
            'variants_per_split': self.variants_per_split,
            'total_splits': self.total_splits,
            'flank_length': self.flank_length,
            'gramtools_kmer_size': self.gramtools_kmer_size,
            'max_read_length': self.max_read_length,
            'total_split_files': self.total_split_files,
            'split_files': self.vcf_split_files,
            'total_input_records': self.total_input_records,
        }
        with open(self.metadata_pickle, 'wb') as f:
            pickle.dump(metadata, f, pickle.HIGHEST_PROTOCOL)


    def _load_existing_data(self):
        with open(self.metadata_pickle, 'rb') as f:
            metadata = pickle.load(f)

        self.vcf_infile = metadata['vcf_infile']
        self.ref_fasta = metadata['ref_fasta']
        self.variants_per_split = metadata['variants_per_split']
        self.total_splits = metadata['total_splits']
        self.flank_length = metadata['flank_length']
        self.gramtools_kmer_size = metadata['gramtools_kmer_size']
        self.max_read_length = metadata['max_read_length']
        self.total_split_files = metadata['total_split_files']
        self.vcf_split_files = metadata['split_files']
        self.total_input_records = metadata['total_input_records']
        logging.info('Loaded existing data from chunked VCF directory ' + self.outdir)


    @classmethod
    def _chunk_end_indexes_from_vcf_record_list(cls, record_list, start_index, total_sites, flank_length):
        '''Returns tuple of:
           1. last index of VCF record that we want to use for variant calling
           2. index of last variant in the chunk, which can't be used for variant calling
              but should end up in the gramtools graph'''
        file_start_index = start_index
        while file_start_index > 0:
            distance_to_previous_variant = record_list[start_index].POS - record_list[file_start_index - 1].ref_end_pos()
            if distance_to_previous_variant > flank_length:
                break
            file_start_index -= 1

        use_vcf_end_index = min(start_index + total_sites - 1, len(record_list) - 1)
        if use_vcf_end_index == len(record_list) - 1:
            return file_start_index, use_vcf_end_index, use_vcf_end_index

        file_end_index = use_vcf_end_index

        while file_end_index < len(record_list) - 1:
            distance_to_next_variant = record_list[file_end_index + 1].POS - record_list[use_vcf_end_index].ref_end_pos()
            if distance_to_next_variant > flank_length:
                break
            file_end_index += 1

        return file_start_index, use_vcf_end_index, min(file_end_index, len(record_list) - 1)


    @classmethod
    def _total_variants_in_vcf_dict(cls, vcf_dict):
        return sum([len(x) for x in vcf_dict.values()])


    def make_split_files(self):
        if len(self.vcf_split_files) > 0:
            return

        self.total_split_files = 0
        self.total_input_records = 0
        vcf_header_lines, vcf_records = cluster_vcf_records.vcf_file_read.vcf_file_to_dict(self.vcf_infile)
        if self.variants_per_split is None:
            total_records = VcfChunker._total_variants_in_vcf_dict(vcf_records)
            self.variants_per_split = 1 + int(total_records / self.total_splits)

        for ref_name, vcf_record_list in vcf_records.items():
            file_end_index = -1
            use_end_index = -1
            self.vcf_split_files[ref_name] = []
            self.total_input_records += len(vcf_record_list)

            while use_end_index < len(vcf_record_list) - 1:
                if file_end_index == -1:
                    use_start_index = 0
                else:
                    use_start_index = self.vcf_split_files[ref_name][-1].use_end_index + 1

                file_start_index, use_end_index, file_end_index = VcfChunker._chunk_end_indexes_from_vcf_record_list(vcf_record_list, use_start_index, self.variants_per_split, self.flank_length)
                split_file = SplitFile(
                    os.path.join(self.outdir, 'split.' + str(self.total_split_files) + '.in.vcf'),
                    self.total_split_files,
                    ref_name,
                    max(0, min(vcf_record_list[file_start_index].POS, vcf_record_list[use_start_index].POS - self.flank_length)),
                    max(vcf_record_list[file_end_index].ref_end_pos(), vcf_record_list[use_end_index].ref_end_pos() + self.flank_length),
                    file_start_index,
                    file_end_index,
                    use_start_index,
                    use_end_index,
                    os.path.join(self.outdir, 'split.' + str(self.total_split_files) + '.gramtools_build'),
                )

                self.vcf_split_files[ref_name].append(split_file)

                with open(split_file.filename, 'w') as f:
                    print(*vcf_header_lines, sep='\n', file=f)
                    for i in range(file_start_index, file_end_index + 1, 1):
                        print(vcf_record_list[i], file=f)

                gramtools.run_gramtools_build(split_file.gramtools_build_dir, split_file.filename, self.ref_fasta, self.max_read_length, self.gramtools_kmer_size)
                self.total_split_files += 1
                logging.info('Made split VCF file ' + split_file.filename + ' and ran gramtools build. Total split files: ' + str(self.total_split_files))

        self._save_metadata()


    def merge_files(self, files_to_merge, outfile):
        total_output_records = 0
        printed_header_lines = False
        logging.info('Making merged VCF file ' + outfile)

        with open(outfile, 'w') as f:
            for ref_name in self.vcf_split_files:
                assert ref_name in files_to_merge
                assert len(self.vcf_split_files[ref_name]) == len(files_to_merge[ref_name])
                for i, split_file in enumerate(self.vcf_split_files[ref_name]):
                    header_lines, records_to_merge = cluster_vcf_records.vcf_file_read.vcf_file_to_list(files_to_merge[ref_name][i])
                    if not printed_header_lines:
                        print(*header_lines, sep='\n', file=f)
                        printed_header_lines = True
                    start_i = split_file.use_start_index - split_file.file_start_index
                    end_i = start_i + split_file.use_end_index - split_file.use_start_index
                    for j in range(start_i, end_i + 1, 1):
                        total_output_records += 1
                        print(records_to_merge[j], file=f)

        if self.total_input_records != total_output_records:
            raise Error('Number of input VCF records = ' + str(self.total_input_records) + ' != ' + str(total_output_records) + ' = numnber of output VCF records. Cannot continue')

        logging.info('Finished making merged VCF file. Total records: ' + str(self.total_input_records))
