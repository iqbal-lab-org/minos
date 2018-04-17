from collections import namedtuple
import os

import cluster_vcf_records

class Error (Exception): pass


split_file_attributes = [
    'filename',
    'chrom',
    'chrom_start',
    'chrom_end',
    'file_start_index',
    'file_end_index',
    'use_start_index',
    'use_end_index',
]
SplitFile = namedtuple('SplitFile', split_file_attributes)


class VcfChunker:
    def __init__(self, vcf_infile, outdir, variants_per_split=None, total_splits=100, flank_length=200):
        self.vcf_infile = os.path.abspath(vcf_infile)
        self.outdir = os.path.abspath(outdir)
        self.variants_per_split = variants_per_split
        self.total_splits = total_splits
        self.flank_length = flank_length
        if not os.path.exists(self.vcf_infile):
            raise Error('VCF file not found: ' + self.vcf_infile)

        try:
            os.mkdir(self.outdir)
        except:
            raise Error('Error mkdir ' + self.outdir)

        self.vcf_split_files = {} # ref name -> list of SplitFile


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
        assert len(self.vcf_split_files) == 0
        total_files = 0
        vcf_header_lines, vcf_records = cluster_vcf_records.vcf_file_read.vcf_file_to_dict(self.vcf_infile)
        if self.variants_per_split is None:
            total_records = VcfChunker._total_variants_in_vcf_dict(vcf_records)
            self.variants_per_split = 1 + int(total_records / self.total_splits)

        for ref_name, vcf_record_list in vcf_records.items():
            file_end_index = -1
            self.vcf_split_files[ref_name] = []

            while file_end_index < len(vcf_record_list) - 1:
                if file_end_index == -1:
                    use_start_index = 0
                else:
                    use_start_index = self.vcf_split_files[ref_name][-1].use_end_index + 1

                file_start_index, use_end_index, file_end_index = VcfChunker._chunk_end_indexes_from_vcf_record_list(vcf_record_list, use_start_index, self.variants_per_split, self.flank_length)
                split_file = SplitFile(
                    os.path.join(self.outdir, 'split.' + str(total_files) + '.in.vcf'),
                    ref_name,
                    max(0, min(vcf_record_list[file_start_index].POS, vcf_record_list[use_start_index].POS - self.flank_length)),
                    max(vcf_record_list[file_end_index].ref_end_pos(), vcf_record_list[use_end_index].ref_end_pos() + self.flank_length),
                    file_start_index,
                    file_end_index,
                    use_start_index,
                    use_end_index,
                )

                self.vcf_split_files[ref_name].append(split_file)

                with open(split_file.filename, 'w') as f:
                    print(*vcf_header_lines, sep='\n', file=f)
                    for i in range(file_start_index, file_end_index + 1, 1):
                        print(vcf_record_list[i], file=f)

                total_files += 1

