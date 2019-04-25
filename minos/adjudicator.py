import json
import logging
import os
import shutil
import statistics
import sys

from cluster_vcf_records import vcf_clusterer, vcf_file_read

from minos import bam_read_extract, dependencies, genotype_confidence_simulator, gramtools, plots, utils, vcf_chunker

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
        gramtools_kmer_size=10,
        sample_name=None,
        variants_per_split=None,
        alleles_per_split=None,
        total_splits=None,
        clean=True,
        genotype_simulation_iterations=10000,
        use_unmapped_reads=False,
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
        self.plots_prefix = os.path.join(self.outdir, 'final.vcf.plots')

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
            if len(self.vcf_files) != 1:
                raise Error('Error! If gramtools_build_dir is used, then There Can Be Only One input VCF file (which is assumed to be clustered')


        self.gramtools_kmer_size = gramtools_kmer_size
        self.gramtools_quasimap_dir = os.path.join(self.outdir, 'gramtools.quasimap')
        self.perl_generated_vcf = os.path.join(self.gramtools_build_dir, 'perl_generated_vcf')
        self.read_error_rate = read_error_rate
        self.max_read_length = max_read_length
        self.variants_per_split = variants_per_split
        self.alleles_per_split = alleles_per_split
        self.total_splits = total_splits

        if (self.total_splits is not None or self.variants_per_split is not None or self.alleles_per_split is not None) and len(self.reads_files) != 1:
            raise Error('Error! If using splitting, must input one reads file (which is assumed to be a sorted indexed BAM file)')

        self.clean = clean
        self.genotype_simulation_iterations = genotype_simulation_iterations
        self.use_unmapped_reads = use_unmapped_reads


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
            return 10
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

        if self.user_supplied_gramtools_build_dir:
            logging.info('User supplied gramtools build dir. Assuming VCF already clustered, so skipping clustering')
            assert len(self.vcf_files) == 1
            self.clustered_vcf = self.vcf_files[0]
        else:
            logging.info('Clustering VCF file(s), to make one VCF input file for gramtools')
            clusterer = vcf_clusterer.VcfClusterer(
                self.vcf_files,
                self.ref_fasta,
                self.clustered_vcf,
                max_distance_between_variants=1,
                max_alleles_per_cluster=self.max_alleles_per_cluster,
            )
            clusterer.run()

            logging.info('Finished clustering VCF file(s)')

        if not vcf_file_read.vcf_file_has_at_least_one_record(self.clustered_vcf):
            error_message = 'No VCF records. Cannot continue. Please check that the input VCF files contained at least one variant'
            logging.error(error_message)
            raise Error(error_message)


        if self.total_splits is not None or self.variants_per_split is not None or self.alleles_per_split is not None or os.path.exists(os.path.join(self.split_input_dir, 'data.pickle')):
            self._run_gramtools_with_split_vcf()
        else:
            self._run_gramtools_not_split_vcf()

        logging.info('Making plots from final.vcf')
        plots.plots_from_minos_vcf(self.final_vcf, self.plots_prefix)

        logging.info('All done! Thank you for using minos :)')


    @classmethod
    def _add_gt_conf_percentile_and_filters_to_vcf_file(cls, vcf_file, mean_depth, depth_variance, error_rate, iterations, min_dp=2, min_gt_conf_percentile=2.5):
        '''Overwrites vcf_file, with new version that has GT_CONF_PERCENTILE added,
        and filter for DP and GT_CONF_PERCENTILE'''
        simulations = genotype_confidence_simulator.GenotypeConfidenceSimulator(mean_depth, depth_variance, error_rate, allele_length=1, iterations=iterations)
        simulations.run_simulations()
        vcf_header, vcf_lines = vcf_file_read.vcf_file_to_list(vcf_file)
        for i, line in enumerate(vcf_header):
            if line.startswith('##FORMAT=<ID=GT_CONF'):
                break
        else:
            raise Exception(f'No GT_CONF description found in header of VCF file {vcf_file}. Cannot continue')

        vcf_header.insert(i+1, r'''##FORMAT=<ID=GT_CONF_PERCENTILE,Number=1,Type=Float,Description="Percentile of GT_CONF">''')
        vcf_header.insert(i+1, f'##FILTER=<ID=MIN_DP,Description="Minimum DP of {min_dp}">')
        vcf_header.insert(i+1, f'##FILTER=<ID=MIN_GCP,Description="Minimum GT_CONF_PERCENTILE of {min_gt_conf_percentile}">')

        with open(vcf_file, 'w') as f:
            print(*vcf_header, sep='\n', file=f)

            for vcf_record in vcf_lines:
                vcf_record.FILTER = set()

                if 'GT_CONF' in vcf_record.FORMAT:
                    conf = int(round(float(vcf_record.FORMAT['GT_CONF'])))
                    if 'GT' in vcf_record.FORMAT and '.' not in vcf_record.FORMAT['GT']:
                        vcf_record.set_format_key_value('GT_CONF_PERCENTILE', str(simulations.get_percentile(conf)))
                        if 'DP' in vcf_record.FORMAT and float(vcf_record.FORMAT['DP']) < min_dp:
                            vcf_record.FILTER.add('MIN_DP')
                        if float(vcf_record.FORMAT['GT_CONF_PERCENTILE']) < min_gt_conf_percentile:
                            vcf_record.FILTER.add('MIN_GCP')
                        if len(vcf_record.FILTER) == 0:
                            vcf_record.FILTER.add('PASS')

                print(vcf_record, file=f)


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
        mean_depth, depth_variance, vcf_header, vcf_records, allele_coverage, allele_groups = gramtools.load_gramtools_vcf_and_allele_coverage_files(self.perl_generated_vcf, self.gramtools_quasimap_dir)
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

        logging.info(f'Adding GT_CONF_PERCENTLE to debug VCF file {self.unfiltered_vcf_file}, using mean depth {mean_depth}, depth variance {depth_variance}, error rate {self.read_error_rate}, and {self.genotype_simulation_iterations} simulation iterations')
        Adjudicator._add_gt_conf_percentile_and_filters_to_vcf_file(self.unfiltered_vcf_file, mean_depth, depth_variance, self.read_error_rate, self.genotype_simulation_iterations)
        logging.info(f'Adding GT_CONF_PERCENTLE to final VCF file {self.final_vcf}, using mean depth {mean_depth}, depth variance {depth_variance}, error rate {self.read_error_rate}, and {self.genotype_simulation_iterations} simulation iterations')
        Adjudicator._add_gt_conf_percentile_and_filters_to_vcf_file(self.final_vcf, mean_depth, depth_variance, self.read_error_rate, self.genotype_simulation_iterations)

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
            alleles_per_split=self.alleles_per_split,
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

        if self.use_unmapped_reads:
            unmapped_reads_file = os.path.join(self.split_output_dir, 'unmapped_reads.bam')
            bam_read_extract.get_unmapped_reads(self.reads_files[0], unmapped_reads_file)

        split_vcf_outfiles = {}
        split_vcf_outfiles_unfiltered = {}
        mean_depths = []
        depth_variances = []

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
                if self.use_unmapped_reads:
                    reads_files = [unmapped_reads_file, split_reads_file]
                else:
                    reads_files = [split_reads_file]

                build_report, quasimap_report = gramtools.run_gramtools(
                    split_file.gramtools_build_dir,
                    gramtools_quasimap_dir,
                    split_file.filename,
                    self.ref_fasta,
                    reads_files,
                    self.max_read_length,
                    kmer_size=self.gramtools_kmer_size,
                )

                logging.info('Loading split gramtools quasimap output files ' + gramtools_quasimap_dir)
                perl_generated_vcf = os.path.join(split_file.gramtools_build_dir, 'perl_generated_vcf')
                mean_depth, depth_variance, vcf_header, vcf_records, allele_coverage, allele_groups = gramtools.load_gramtools_vcf_and_allele_coverage_files(perl_generated_vcf, gramtools_quasimap_dir)
                mean_depths.append(mean_depth)
                depth_variances.append(depth_variance)
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

        mean_depth = statistics.mean(mean_depths)
        depth_variance = statistics.mean(depth_variances)
        logging.info(f'Adding GT_CONF_PERCENTLE to debug VCF file {self.unfiltered_vcf_file}, using mean depth {mean_depth}, depth variance {depth_variance}, error rate {self.read_error_rate}, and {self.genotype_simulation_iterations} simulation iterations')
        Adjudicator._add_gt_conf_percentile_and_filters_to_vcf_file(self.unfiltered_vcf_file, mean_depth, depth_variance, self.read_error_rate, self.genotype_simulation_iterations)
        logging.info(f'Adding GT_CONF_PERCENTLE to final VCF file {self.final_vcf}, using mean depth {mean_depth}, depth variance {depth_variance}, error rate {self.read_error_rate}, and {self.genotype_simulation_iterations} simulation iterations')
        Adjudicator._add_gt_conf_percentile_and_filters_to_vcf_file(self.final_vcf, mean_depth, depth_variance, self.read_error_rate, self.genotype_simulation_iterations)

        if self.clean:
            logging.info('Deleting temp split VCF files')
            for d in split_vcf_outfiles, split_vcf_outfiles_unfiltered:
                for file_list in d.values():
                    for filename in file_list:
                        os.unlink(filename)
            if self.use_unmapped_reads:
                os.unlink(unmapped_reads_file)

