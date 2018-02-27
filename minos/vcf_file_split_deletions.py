import logging
import os
import shutil

from cluster_vcf_records import vcf_file_read

class Error (Exception): pass

class VcfFileSplitDeletions:
    def __init__(self,
        infile,
        outfile_small_vars,
        outfile_long_deletions,
        min_large_ref_length=50,
    ):
        self.infile = infile
        self.outfile_small_vars = outfile_small_vars
        self.outfile_long_deletions = outfile_long_deletions
        self.min_large_ref_length = min_large_ref_length


    def run(self):
        header_lines, vcf_records = vcf_file_read.vcf_file_to_list(self.infile)

        with open(self.outfile_small_vars, 'w') as f_small, open(self.outfile_long_deletions, 'w') as f_big:
            print(*header_lines, sep='\n', file=f_small)
            print(*header_lines, sep='\n', file=f_big)

            for record in vcf_records:
                if len(record.REF) < self.min_large_ref_length:
                    print(record, file=f_small)
                else:
                    print(record, file=f_big)

