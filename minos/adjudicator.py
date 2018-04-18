import json
import logging
import os
import shutil
import sys

from cluster_vcf_records import vcf_clusterer, vcf_file_read

from minos import bam_read_extract, dependencies, gramtools, utils, vcf_chunker

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
        gramtools_kmer_size=15,
        sample_name=None,
        variants_per_split=None,
        total_splits=None,
        clean=True,
    ):
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.reads_files = [os.path.abspath(x) for x in reads_files]
        self.vcf_files = [os.path.abspath(x) for x in vcf_files]
        self.overwrite_outdir = overwrite_outdir
        self.max_alleles_per_cluster = max_alleles_per_cluster
        self.sample_name = sample_name
        self.outdir = os.path.abspath(outdir)
        self.split_output_dir = os.path.join(self.outdir, 'split.out')
        self.log_file = os.path.join(self.outdir, 'log.txt')
        self.clustered_vcf = os.path.join(self.outdir, 'gramtools.in.vcf')
        self.unfiltered_vcf_file = os.path.join(self.outdir, 'debug.calls_with_zero_cov_alleles.vcf')
        self.final_vcf = os.path.join(self.outdir, 'final.vcf')

        if gramtools_build_dir is None:
            self.split_input_dir = os.path.join(self.outdir, 'split.in')
            self.gramtools_build_dir = os.path.join(self.outdir, 'gramtools.build')
            self.user_supplied_gramtools_build_dir = False
        else:
            self.gramtools_build_dir = os.path.abspath(gramtools_build_dir)
            self.split_input_dir = self.gramtools_build_dir
            self.user_supplied_gramtools_build_dir = True
            if not os.path.exists(self.gramtools_build_dir):
                raise Error('Error! gramtools_build_dir=' + self.gramtools_build_dir + ' used, but directory not found on disk. Cannot continue')

        self.gramtools_kmer_size = gramtools_kmer_size
        self.gramtools_quasimap_dir = os.path.join(self.outdir, 'gramtools.quasimap')
        self.perl_generated_vcf = os.path.join(self.gramtools_build_dir, 'perl_generated_vcf')
        self.read_error_rate = read_error_rate
        self.max_read_length = max_read_length
        self.variants_per_split = variants_per_split
        self.total_splits = total_splits

        if (self.total_splits is not None or self.variants_per_split is not None) and len(self.reads_files) != 1:
            raise Error('Error! If using splitting, must input one reads file (which is assumed to be a sorted indexed BAM file)')

        self.clean = clean


    @classmethod
    def _get_gramtools_kmer_size(cls, build_dir, input_kmer_size):
        if build_dir is not None and os.path.exists(build_dir):
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
        logging.info('Command run: ' + ' '.join(sys.argv))
        dependencies.check_and_report_dependencies(programs=['gramtools'])
        logging.info('Dependencies look OK')

        if self.read_error_rate is None or self.max_read_length is None:
            logging.info('One or both of read_error_rate and max_read_length not known. Estimate from first 10,000 reads...')
            estimated_read_length, estimated_read_error_rate = utils.estimate_max_read_length_and_read_error_rate_from_qual_scores(self.reads_files[0])
            logging.info('Estimated max_read_length=' + str(estimated_read_length) + ' and read_error_rate=' + str(estimated_read_error_rate))

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


        if self.total_splits is not None or self.variants_per_split is not None or os.path.exists(os.path.join(self.split_input_dir, 'data.pickle')):
            self._run_gramtools_with_split_vcf()
        else:
            self._run_gramtools_not_split_vcf()

        logging.info('All done! Thank you for using minos :)')


    def _run_gramtools_not_split_vcf(self):
        self.gramtools_kmer_size = Adjudicator._get_gramtools_kmer_size(self.gramtools_build_dir, self.gramtools_kmer_size)
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
            self.unfiltered_vcf_file,
            self.gramtools_kmer_size,
            sample_name=sample_name,
            max_read_length=self.max_read_length,
            filtered_outfile=self.final_vcf
        )

        if self.clean:
            os.rename(os.path.join(self.gramtools_quasimap_dir, 'report.json'), os.path.join(self.outdir, 'gramtools.quasimap.report.json'))
            shutil.rmtree(self.gramtools_quasimap_dir)

            if not self.user_supplied_gramtools_build_dir:
                os.rename(os.path.join(self.gramtools_build_dir, 'build_report.json'), os.path.join(self.outdir, 'gramtools.build.report.json'))
                shutil.rmtree(self.gramtools_build_dir)


    def _run_gramtools_with_split_vcf(self):
        logging.info('Splitting VCF files into chunks (if not already done)')
        chunker = vcf_chunker.VcfChunker(
            self.split_input_dir,
            vcf_infile=self.clustered_vcf,
            ref_fasta=self.ref_fasta,
            variants_per_split=self.variants_per_split,
            max_read_length=self.max_read_length,
            total_splits=self.total_splits,
            flank_length=self.max_read_length,
            gramtools_kmer_size=self.gramtools_kmer_size,
        )
        chunker.make_split_files()
        self.gramtools_kmer_size = chunker.gramtools_kmer_size

        logging.info('VCF file split into ' + str(chunker.total_split_files) + ' chunks')
        try:
            os.mkdir(self.split_output_dir)
        except:
            raise Error('Error making output split directory ' + self.split_output_dir)

        unmapped_reads_file = os.path.join(self.split_output_dir, 'unmapped_reads.bam')
        bam_read_extract.get_unmapped_reads(self.reads_files[0], unmapped_reads_file)
        split_vcf_outfiles = {}
        split_vcf_outfiles_unfiltered = {}

        for ref_name, split_file_list in chunker.vcf_split_files.items():
            split_vcf_outfiles[ref_name] = []
            split_vcf_outfiles_unfiltered[ref_name] = []
            for split_file in split_file_list:
                logging.info('===== Start analysing variants in VCF split file ' + split_file.filename + ' =====')
                split_reads_file = os.path.join(self.split_output_dir, 'split.' + str(split_file.file_number) + '.reads.bam')
                bam_read_extract.get_region(
                    self.reads_files[0],
                    split_file.chrom,
                    split_file.chrom_start,
                    split_file.chrom_end,
                    split_reads_file,
                )

                gramtools_quasimap_dir = os.path.join(self.split_output_dir, 'split.' + str(split_file.file_number) + '.gramtools.quasimap')

                build_report, quasimap_report = gramtools.run_gramtools(
                    split_file.gramtools_build_dir,
                    gramtools_quasimap_dir,
                    split_file.filename,
                    self.ref_fasta,
                    [unmapped_reads_file, split_reads_file],
                    self.max_read_length,
                    kmer_size=self.gramtools_kmer_size,
                )

                logging.info('Loading split gramtools quasimap output files ' + gramtools_quasimap_dir)
                perl_generated_vcf = os.path.join(split_file.gramtools_build_dir, 'perl_generated_vcf')
                mean_depth, vcf_header, vcf_records, allele_coverage, allele_groups = gramtools.load_gramtools_vcf_and_allele_coverage_files(perl_generated_vcf, gramtools_quasimap_dir)
                logging.info('Finished loading gramtools files')
                if self.sample_name is None:
                    sample_name = vcf_file_read.get_sample_name_from_vcf_header_lines(vcf_header)
                else:
                    sample_name = self.sample_name
                assert sample_name is not None
                split_vcf_out = os.path.join(self.split_output_dir, 'split.' + str(split_file.file_number) + '.out.vcf')
                unfiltered_vcf_out = os.path.join(self.split_output_dir, 'split.' + str(split_file.file_number) + '.out.debug.calls_with_zero_cov_alleles.vcf')
                logging.info('Writing VCf output file ' + split_vcf_out + ' for split VCF file ' + split_file.filename)
                gramtools.write_vcf_annotated_using_coverage_from_gramtools(
                    mean_depth,
                    vcf_records,
                    allele_coverage,
                    allele_groups,
                    self.read_error_rate,
                    unfiltered_vcf_out,
                    self.gramtools_kmer_size,
                    sample_name=sample_name,
                    max_read_length=self.max_read_length,
                    filtered_outfile=split_vcf_out,
                )
                split_vcf_outfiles[ref_name].append(split_vcf_out)
                split_vcf_outfiles_unfiltered[ref_name].append(unfiltered_vcf_out)

                if self.clean:
                    logging.info('Cleaning gramtools files from split VCF file ' + split_file.filename)
                    if not self.user_supplied_gramtools_build_dir:
                        os.rename(os.path.join(split_file.gramtools_build_dir, 'build_report.json'), split_file.gramtools_build_dir + '.report.json')
                        shutil.rmtree(split_file.gramtools_build_dir)
                        os.unlink(split_file.filename)

                    os.rename(os.path.join(gramtools_quasimap_dir, 'report.json'), gramtools_quasimap_dir + '.report.json')
                    shutil.rmtree(gramtools_quasimap_dir)
                    os.unlink(split_reads_file)

                logging.info('===== Finish analysing variants in VCF split file ' + split_file.filename + ' =====')

        logging.info('Merging VCF files into one output file ' + self.final_vcf)
        chunker.merge_files(split_vcf_outfiles, self.final_vcf)
        chunker.merge_files(split_vcf_outfiles_unfiltered, self.unfiltered_vcf_file)

        if self.clean:
            logging.info('Deleting temp split VCF files')
            for d in split_vcf_outfiles, split_vcf_outfiles_unfiltered:
                for file_list in d.values():
                    for filename in file_list:
                        os.unlink(filename)
            os.unlink(unmapped_reads_file)

