import json
import logging
import os
import shutil

from cluster_vcf_records import vcf_clusterer, vcf_file_read

from minos import dependencies, gramtools, utils

class Error (Exception): pass

class Adjudicator:
    def __init__(self,
        outdir,
        ref_fasta,
        reads_files,
        vcf_files,
        max_read_length=None,
        read_error_rate=None,
        overwrite_outdir=False,
        max_alleles_per_cluster=5000,
        gramtools_build_dir=None,
        gramtools_kmer_size=None,
        sample_name=None,
    ):
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.reads_files = [os.path.abspath(x) for x in reads_files]
        self.vcf_files = [os.path.abspath(x) for x in vcf_files]
        self.overwrite_outdir = overwrite_outdir
        self.max_alleles_per_cluster = max_alleles_per_cluster
        self.sample_name = sample_name
        self.outdir = os.path.abspath(outdir)
        self.log_file = os.path.join(self.outdir, 'log.txt')
        self.clustered_vcf = os.path.join(self.outdir, 'gramtools.in.vcf')
        self.final_vcf = os.path.join(self.outdir, 'final.vcf')

        if gramtools_build_dir is None:
            self.gramtools_build_dir = os.path.join(self.outdir, 'gramtools.build')
        else:
            self.gramtools_build_dir = os.path.abspath(gramtools_build_dir)
            if not os.path.exists(self.gramtools_build_dir):
                raise Error('Error! gramtools_build_dir=' + self.gramtools_build_dir + ' used, but directory not found on disk. Cannot continue')

        self.gramtools_kmer_size = Adjudicator._get_gramtools_kmer_size(gramtools_build_dir, gramtools_kmer_size)
        self.gramtools_quasimap_dir = os.path.join(self.outdir, 'gramtools.quasimap')
        self.perl_generated_vcf = os.path.join(self.gramtools_build_dir, 'perl_generated_vcf')
        self.read_error_rate = read_error_rate
        self.max_read_length = max_read_length



    @classmethod
    def _get_gramtools_kmer_size(cls, build_dir, input_kmer_size):
        if build_dir is not None:
            json_file = os.path.join(build_dir, 'build_report.json')
            with open(json_file) as f:
                json_build_report = json.load(f)

            if 'kmer_size' in json_build_report:
                kmer_size_from_json = json_build_report['kmer_size']
            else:
                raise Error('kmer_size not found in graamtools build report' + json_file)

            if input_kmer_size is not None:
                logging.warning('gramtools_kmer_size option used, but going to ignore it because gramtools_build_dir used, which means that gramtools build has already been run')

            return kmer_size_from_json
        elif input_kmer_size is None:
            return 15
        else:
            return input_kmer_size


    def run(self):
        if os.path.exists(self.outdir) and self.overwrite_outdir:
            shutil.rmtree(self.outdir)

        try:
            os.mkdir(self.outdir)
        except:
            raise Error('Error making output directory ' + self.outdir)

        fh = logging.FileHandler(self.log_file, mode='w')
        log = logging.getLogger()
        formatter = logging.Formatter('[minos %(asctime)s %(levelname)s] %(message)s', datefmt='%d-%m-%Y %H:%M:%S')
        fh.setFormatter(formatter)
        log.addHandler(fh)
        dependencies.check_and_report_dependencies(programs=['gramtools'])
        logging.info('Dependencies look OK')

        if self.read_error_rate is None or self.max_read_length is None:
            logging.info('One or both of read_error_rate and max_read_length not known. Estimate from first 10,000 reads...')
            estimated_read_length, estimated_read_error_rate = utils.estimate_max_read_length_and_read_error_rate_from_qual_scores(self.reads_files[0])
            logging.info('Estaimted max_read_length=' + str(estimated_read_length) + ' and read_error_rate=' + str(estimated_read_length))

        self.read_error_rate = estimated_read_error_rate if self.read_error_rate is None else self.read_error_rate
        self.max_read_length = estimated_read_length if self.max_read_length is None else self.max_read_length
        logging.info('Using max_read_length=' + str(self.max_read_length) + ' and read_error_rate=' + str(self.read_error_rate))


        clusterer = vcf_clusterer.VcfClusterer(
            self.vcf_files,
            self.ref_fasta,
            self.clustered_vcf,
            max_distance_between_variants=1,
            max_alleles_per_cluster=self.max_alleles_per_cluster,
        )
        clusterer.run()

        logging.info('Clustered VCF files, to make one VCF input file for gramtools')

        if not vcf_file_read.vcf_file_has_at_least_one_record(self.clustered_vcf):
            error_message = 'No VCF records. Cannot continue. Please check that the input VCF files contained at least one variant'
            logging.error(error_message)
            raise Error(error_message)

        build_report, quasimap_report = gramtools.run_gramtools(
            self.gramtools_build_dir,
            self.gramtools_quasimap_dir,
            self.clustered_vcf,
            self.ref_fasta,
            self.reads_files,
            self.max_read_length,
            kmer_size=self.gramtools_kmer_size,
        )

        logging.info('Loading gramtools quasimap output files ' + self.gramtools_quasimap_dir)
        mean_depth, vcf_header, vcf_records, allele_coverage, allele_groups = gramtools.load_gramtools_vcf_and_allele_coverage_files(self.perl_generated_vcf, self.gramtools_quasimap_dir)
        logging.info('Finished loading gramtools files')
        if self.sample_name is None:
            sample_name = vcf_file_read.get_sample_name_from_vcf_header_lines(vcf_header)
        else:
            sample_name = self.sample_name
        assert sample_name is not None
        logging.info('Writing VCf output file ' + self.final_vcf)
        gramtools.write_vcf_annotated_using_coverage_from_gramtools(
            mean_depth,
            vcf_records,
            allele_coverage,
            allele_groups,
            self.read_error_rate,
            self.final_vcf,
            self.gramtools_kmer_size,
            sample_name=sample_name,
            max_read_length=self.max_read_length,
        )

        logging.info('All done! Thank you for using minos :)')
