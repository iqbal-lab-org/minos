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


    @classmethod
    def _load_input_data_tsv(cls, infile):
        logging.info('Start reading file ' + infile)
        data = []
        with open(infile) as f:
            for line in f:
                try:
                    vcf_file, *reads_files = line.rstrip().split('\t')
                except:
                    raise Error('Bad line in input TSV file: ' + line.rstrip())

                if not(os.path.exists(vcf_file)):
                    raise Error('VCF file not found: ' + vcf_file)
                for reads_file in reads_files:
                    if not(os.path.exists(reads_file)):
                        raise Error('Reads file not found: ' + reads_file)

                data.append((os.path.abspath(vcf_file), [os.path.abspath(x) for x in reads_files]))

        logging.info('Finish reading file ' + infile + '. Loaded ' + str(len(data)) + 'samples')
        return data


    @classmethod
    def _write_nextflow_data_tsv(cls, data, outfile):
        with open(outfile, 'w') as f:
            print('sample_id', 'vcf_file', 'reads_files', sep='\t', file=f)
            for i, (vcf_file, reads_files) in enumerate(data):
                print(i, vcf_file, ' '.join(reads_files), sep='\t', file=f)


    def _write_nextflow_script(self):
        pass


    def _prepare_nextflow_input_files(self):
        # parse input_data_tsv
        # check files exist
        # make files for input to nextflow
        pass


