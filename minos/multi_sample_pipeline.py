import logging
import os
import shutil

from cluster_vcf_records import vcf_file_read

class Error (Exception): pass

class MultiSamplePipeline:
    def __init__(self,
        input_data_tsv,
        output_dir,
        nextflow_config_file,
        min_large_ref_length=50,
    ):
        self.input_data_tsv = os.apth.abspath(input_data_tsv)
        self.output_dir = os.path.abspath(output_dir)
        self.nextflow_config_file = None if nextflow_config_file is None else os.path.abspath(nextflow_config_file)
        self.min_large_ref_length = min_large_ref_length


    def _write_nextflow_script(self):
        pass


    def _prepare_nextflow_input_files(self):
        # parse input_data_tsv
        # check files exist
        # make files for input to nextflow
        pass


