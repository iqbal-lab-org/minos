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
        max_read_length=150,
        read_error_rate=None,
        overwrite_outdir=False,
        max_alleles_per_cluster=5000,
    ):
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.reads_files = [os.path.abspath(x) for x in reads_files]
        self.vcf_files = [os.path.abspath(x) for x in vcf_files]
        self.max_read_length = max_read_length
        self.overwrite_outdir = overwrite_outdir
        self.max_alleles_per_cluster = max_alleles_per_cluster

        self.outdir = os.path.abspath(outdir)
        self.log_file = os.path.join(self.outdir, 'log.txt')
        self.clustered_vcf = os.path.join(self.outdir, 'gramtools.in.vcf')
        self.final_vcf = os.path.join(self.outdir, 'final.vcf')
        self.gramtools_outdir = os.path.join(self.outdir, 'gramtools.out')
        self.perl_generated_vcf = os.path.join(self.gramtools_outdir, 'perl_generated_vcf')

        if read_error_rate is None:
            self.read_error_rate = utils.estimate_read_error_rate_from_qual_scores(self.reads_files[0])
        else:
            self.read_error_rate = read_error_rate


    def run(self):
        if os.path.exists(self.outdir) and self.overwrite_outdir:
            shutil.rmtree(self.outdir)

        try:
            os.mkdir(self.outdir)
        except:
            raise Error('Error making output directory ' + self.outdir)

        dependencies.check_and_report_dependencies(self.log_file, programs=['gramtools'])

        clusterer = vcf_clusterer.VcfClusterer(
            self.vcf_files,
            self.ref_fasta,
            self.clustered_vcf,
            max_distance_between_variants=1,
            max_alleles_per_cluster=self.max_alleles_per_cluster,
        )
        clusterer.run()

        quasimap_dir = gramtools.run_gramtools(
            self.gramtools_outdir,
            self.clustered_vcf,
            self.ref_fasta,
            self.reads_files,
            self.max_read_length,
        )

        mean_depth, vcf_header, vcf_records, allele_coverage, allele_groups = gramtools.load_gramtools_vcf_and_allele_coverage_files(self.perl_generated_vcf, quasimap_dir)
        sample_name = vcf_file_read.get_sample_name_from_vcf_header_lines(vcf_header)
        assert sample_name is not None
        gramtools.write_vcf_annotated_using_coverage_from_gramtools(
            mean_depth,
            vcf_records,
            allele_coverage,
            allele_groups,
            self.read_error_rate,
            self.final_vcf,
            sample_name=sample_name,
        )

