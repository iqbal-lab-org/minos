import json
import logging
import math
import os
import shutil
import statistics
import sys

from cluster_vcf_records import vcf_clusterer, vcf_file_read, variant_tracking
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pyfastaq
import seaborn as sns


from minos import (
    bam_read_extract,
    dependencies,
    genotype_confidence_simulator,
    gramtools,
    utils,
    vcf_chunker,
)


class Adjudicator:
    """
    Runs gramtools build and quasimap, genotyping, and confidence simulations for a set of vcfs and a fasta ref.
    """

    def __init__(
        self,
        outdir,
        ref_fasta,
        reads_files,
        vcf_files,
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
        filter_min_dp=0,
        filter_min_gcp=5,
        filter_max_dp=3.0,
        filter_min_frs=0.9,
        call_hets=False,
        debug=False,
        cluster_input_vcfs=True,
    ):
        self.original_ref_fasta = os.path.abspath(ref_fasta)
        self.reads_files = [os.path.abspath(x) for x in reads_files]
        self.vcf_files = [os.path.abspath(x) for x in vcf_files]
        self.overwrite_outdir = overwrite_outdir
        self.max_alleles_per_cluster = max_alleles_per_cluster
        self.sample_name = sample_name
        self._set_sample_name()
        assert self.sample_name is not None
        self.outdir = os.path.abspath(outdir)
        self.split_output_dir = os.path.join(self.outdir, "split.out")
        self.log_file = os.path.join(self.outdir, "log.txt")
        self.cluster_dir = os.path.join(self.outdir, "cluster_vcfs")
        self.clustered_vcf_prefix = os.path.join(self.outdir, "gramtools.in")
        self.clustered_vcf = f"{self.clustered_vcf_prefix}.vcf"
        self.unfiltered_vcf_file = os.path.join(
            self.outdir, "debug.calls_with_zero_cov_alleles.vcf"
        )
        self.final_vcf = os.path.join(self.outdir, "final.vcf")

        if gramtools_build_dir is None:
            self.split_input_dir = os.path.join(self.outdir, "split.in")
            self.gramtools_build_dir = os.path.join(self.outdir, "gramtools.build")
            self.user_supplied_gramtools_build_dir = False
        else:
            self.gramtools_build_dir = os.path.abspath(gramtools_build_dir)
            self.split_input_dir = self.gramtools_build_dir
            self.user_supplied_gramtools_build_dir = True
            if not os.path.exists(self.gramtools_build_dir):
                raise Exception(
                    "Error! gramtools_build_dir="
                    + self.gramtools_build_dir
                    + " used, but directory not found on disk. Cannot continue"
                )
            if len(self.vcf_files) != 1:
                raise Exception(
                    "Error! If gramtools_build_dir is used, then There Can Be Only One input VCF file (which is assumed to be clustered"
                )

        self.gramtools_kmer_size = gramtools_kmer_size
        self.gramtools_quasimap_dir = os.path.join(self.outdir, "gramtools.quasimap")
        self.gramtools_quasimap_json = self.gramtools_quasimap_dir + ".report.json"
        self.gramtools_build_json = self.gramtools_build_dir + ".report.json"
        self.read_error_rate = read_error_rate
        self.variants_per_split = variants_per_split
        self.alleles_per_split = alleles_per_split
        self.total_splits = total_splits
        self.mean_depth = None
        self.variance_depth = None

        if (
            self.total_splits is not None
            or self.variants_per_split is not None
            or self.alleles_per_split is not None
        ) and len(self.reads_files) != 1:
            raise Exception(
                "Error! If using splitting, must input one reads file (which is assumed to be a sorted indexed BAM file)"
            )

        self.clean = clean
        self.genotype_simulation_iterations = genotype_simulation_iterations
        self.use_unmapped_reads = use_unmapped_reads
        self.filter_min_dp = filter_min_dp
        self.filter_min_gcp = filter_min_gcp
        self.filter_max_dp = filter_max_dp
        self.filter_min_frs = filter_min_frs
        self.call_hets = call_hets
        if self.call_hets:
            raise NotImplementedError("Heterozygous calling is not implemented")
        self.debug = debug
        self.cluster_input_vcfs = cluster_input_vcfs
        self.ref_seq_lengths = {
            x.id.split()[0]: len(x)
            for x in pyfastaq.sequences.file_reader(self.original_ref_fasta)
        }
        if self.debug:
            self.sim_conf_scores_file = os.path.join(
                self.outdir, "debug.genotype_sim_conf_scores.txt"
            )
            self.real_conf_scores_file = os.path.join(
                self.outdir, "debug.genotype_real_conf_scores.txt"
            )
            self.genotype_hist_pdf = os.path.join(
                self.outdir, "debug.genotype_hist.pdf"
            )
        else:
            self.sim_conf_scores_file = None
            self.real_conf_scores_file = None
            self.genotype_hist_pdf = None

    def build_output_dir(self):
        try:
            os.mkdir(self.outdir)
        except FileExistsError:
            if self.overwrite_outdir:
                shutil.rmtree(self.outdir)
                os.mkdir(self.outdir)
            else:
                raise Exception(
                    f"Output directory {self.outdir} already exists. "
                    f"Rerun command with --force flag if you are OK with overwriting it"
                )
        except Exception as e:
            raise Exception(f"Could not make {self.outdir} due to {e}")

    def _set_sample_name(self):
        if self.sample_name is None:
            try:
                self.sample_name = vcf_file_read.get_sample_name_from_vcf_file(
                    self.vcf_files[0]
                )
            except:
                self.sample_name = "sample"

        if self.sample_name is None:
            self.sample_name = "sample"

    @classmethod
    def _get_gramtools_kmer_size(cls, build_dir, input_kmer_size):
        if build_dir is not None and os.path.exists(build_dir):
            json_file = os.path.join(build_dir, "build_report.json")
            with open(json_file) as f:
                json_build_report = json.load(f)

            if "kmer_size" in json_build_report:
                kmer_size_from_json = json_build_report["kmer_size"]
            else:
                raise Exception(
                    "kmer_size not found in graamtools build report" + json_file
                )

            if input_kmer_size is not None:
                logging.warning(
                    "gramtools_kmer_size option used, but going to ignore it because gramtools_build_dir used, which means that gramtools build has already been run"
                )

            return kmer_size_from_json
        elif input_kmer_size is None:
            return 10
        else:
            return input_kmer_size

    def run(self):
        self.build_output_dir()

        fh = logging.FileHandler(self.log_file, mode="w")
        log = logging.getLogger()
        formatter = logging.Formatter(
            "[minos %(asctime)s %(levelname)s] %(message)s", datefmt="%d-%m-%Y %H:%M:%S"
        )
        fh.setFormatter(formatter)
        log.addHandler(fh)
        logging.info("Command run: " + " ".join(sys.argv))
        to_check = [
            "gramtools",
            "vcfbreakmulti",
            "vcfallelicprimitives",
            "vcfuniq",
            "vt",
        ]
        dependencies.check_and_report_dependencies(programs=to_check)
        logging.info("Dependencies look OK")

        self.ref_fasta = os.path.join(self.outdir, "ref.fa")
        utils.fasta_to_upper_and_ACGT_only(self.original_ref_fasta, self.ref_fasta)

        if self.read_error_rate is None:
            logging.info("read_error_rate unknown. Estimate from first 10,000 reads...")
            (
                estimated_read_length,
                estimated_read_error_rate,
            ) = utils.estimate_max_read_length_and_read_error_rate_from_qual_scores(
                self.reads_files[0]
            )
            logging.info(f"Estimated read_error_rate={estimated_read_error_rate}")

            self.read_error_rate = (
                estimated_read_error_rate
                if self.read_error_rate is None
                else self.read_error_rate
            )
            logging.info(f"Using read_error_rate={self.read_error_rate}")

        if self.user_supplied_gramtools_build_dir:
            logging.info(
                "User supplied gramtools build dir. Assuming VCF already clustered, so skipping clustering"
            )
            assert len(self.vcf_files) == 1
            self.clustered_vcf = self.vcf_files[0]
        elif not self.cluster_input_vcfs:
            logging.info("Skipping VCF clustering because user requested to skip")
        else:
            logging.info(
                "Clustering VCF file(s), to make one VCF input file for gramtools"
            )
            tracker = variant_tracking.VariantTracker(self.cluster_dir, self.ref_fasta)
            tracker.merge_vcf_files(self.vcf_files)
            tracker.cluster(self.clustered_vcf_prefix, float("Inf"), max_alleles=5000)
            if not self.debug:
                os.unlink(f"{self.clustered_vcf_prefix}.excluded.tsv")
                utils.rm_rf(self.cluster_dir)
            utils.remove_vars_from_vcf_at_contig_ends(
                self.clustered_vcf, self.clustered_vcf, ref_lengths=self.ref_seq_lengths
            )
            logging.info("Finished clustering VCF file(s)")

        if not vcf_file_read.vcf_file_has_at_least_one_record(self.clustered_vcf):
            error_message = "No VCF records. Cannot continue. Please check that the input VCF files contained at least one variant"
            logging.error(error_message)
            raise Exception(error_message)

        if (
            self.total_splits is not None
            or self.variants_per_split is not None
            or self.alleles_per_split is not None
            or os.path.exists(os.path.join(self.split_input_dir, "data.pickle"))
        ):
            self._run_gramtools_with_split_vcf()
        else:
            self._run_gramtools_not_split_vcf()

        logging.info("All done! Thank you for using minos :)")

    def _run_gramtools_not_split_vcf(self):
        self.gramtools_kmer_size = Adjudicator._get_gramtools_kmer_size(
            self.gramtools_build_dir, self.gramtools_kmer_size
        )

        self.run_adjudicate(
            self.gramtools_build_dir,
            self.gramtools_quasimap_dir,
            self.clustered_vcf,
            self.reads_files,
            self.final_vcf,
            self.unfiltered_vcf_file,
        )
        self.run_gt_conf()

    def _get_read_coverage_one_split(self, split_file, quasimap_dir):
        all_cov = gramtools.coverage_list_from_quasimap_dir(quasimap_dir)
        # We don't use all the records in a split run of gramtools. The
        # records overlap. There was originally a list of variants for the
        # whole ref seq this split came from. The split_file stores
        # the start/end indexes in that list. It also has the start/end
        # indexes of the records that we want to use. Use these to work
        # out the range we actually want.
        assert 1 + split_file.file_end_index - split_file.file_start_index == len(
            all_cov
        )
        start = split_file.use_start_index - split_file.file_start_index
        number_wanted = 1 + split_file.use_end_index - split_file.use_start_index
        # The list of coverages has None for non-snp sites. We only want
        # the SNPs, so take the correct range and remove the non-snps
        return [x for x in all_cov[start : start + number_wanted] if x is not None]

    def _run_quasimap_one_split(self, split_file, unmapped_reads_file=None):
        logging.info(f"Start quasimap on split file {split_file.filename}")
        split_reads_file = os.path.join(
            self.split_output_dir,
            f"split.{split_file.file_number}.reads.bam",
        )
        bam_read_extract.get_region(
            self.reads_files[0],
            split_file.chrom,
            split_file.chrom_start,
            split_file.chrom_end,
            split_reads_file,
        )

        quasimap_dir = os.path.join(
            self.split_output_dir,
            f"split.{split_file.file_number}.gramtools.quasimap",
        )
        if self.use_unmapped_reads:
            reads_files = [unmapped_reads_file, split_reads_file]
        else:
            reads_files = [split_reads_file]

        build_report, quasimap_report = gramtools.run_gramtools(
            split_file.gramtools_build_dir,
            quasimap_dir,
            split_file.filename,
            self.ref_fasta,
            reads_files,
            kmer_size=self.gramtools_kmer_size,
        )
        read_cov = self._get_read_coverage_one_split(split_file, quasimap_dir)

        if self.clean:
            os.unlink(split_reads_file)

        logging.info(f"Finish quasimap on split file {split_file.filename}")
        return read_cov, build_report, quasimap_report

    def _run_gramtools_with_split_vcf(self):
        logging.info("Splitting VCF files into chunks (if not already done)")
        chunker = vcf_chunker.VcfChunker(
            self.split_input_dir,
            vcf_infile=self.clustered_vcf,
            ref_fasta=self.ref_fasta,
            variants_per_split=self.variants_per_split,
            alleles_per_split=self.alleles_per_split,
            total_splits=self.total_splits,
            gramtools_kmer_size=self.gramtools_kmer_size,
        )
        chunker.make_split_files()
        self.gramtools_kmer_size = chunker.gramtools_kmer_size

        logging.info(
            "VCF file split into " + str(chunker.total_split_files) + " chunks"
        )
        try:
            os.mkdir(self.split_output_dir)
        except:
            raise Exception(
                "Error making output split directory " + self.split_output_dir
            )

        if self.use_unmapped_reads:
            unmapped_reads_file = os.path.join(
                self.split_output_dir, "unmapped_reads.bam"
            )
            bam_read_extract.get_unmapped_reads(
                self.reads_files[0], unmapped_reads_file
            )
        else:
            unmapped_reads_file = None

        read_coverage = []
        build_reports = {}
        quasimap_reports = {}

        # Run gramtools quasimap on each split. Get back the read depth
        # from each split, which we need to get the global read depth and
        # variance, to then use for genotyping
        for ref_name, split_file_list in chunker.vcf_split_files.items():
            for split_file in split_file_list:
                read_cov, build_report, quasimap_report = self._run_quasimap_one_split(
                    split_file, unmapped_reads_file
                )
                read_coverage.extend(read_cov)
                build_reports[split_file.file_number] = build_report
                quasimap_reports[split_file.file_number] = quasimap_report

        with open(self.gramtools_quasimap_json, "w") as f:
            json.dump(quasimap_reports, f, indent=2, sort_keys=True)
        if not self.user_supplied_gramtools_build_dir:
            with open(self.gramtools_build_json, "w") as f:
                json.dump(build_reports, f, indent=2, sort_keys=True)

        self.mean_depth = round(statistics.mean(read_coverage), 3)
        self.variance_depth = round(statistics.variance(read_coverage), 3)

        # Can now genotype each split VCF, using the global mean depth and variance
        split_vcf_outfiles = {}
        split_vcf_outfiles_unfiltered = {}
        for ref_name, split_file_list in chunker.vcf_split_files.items():
            split_vcf_outfiles[ref_name] = []
            split_vcf_outfiles_unfiltered[ref_name] = []
            for split_file in split_file_list:
                build_vcf = os.path.join(split_file.gramtools_build_dir, "build.vcf")
                quasimap_dir = os.path.join(
                    self.split_output_dir,
                    f"split.{split_file.file_number}.gramtools.quasimap",
                )
                logging.info(f"Loading gramtools quasimap output files " + quasimap_dir)
                (
                    _,  # mean depth for this split, which we don't want
                    _,  # depth variance for this split, which we don't want
                    vcf_header,
                    vcf_records,
                    allele_coverage,
                    allele_groups,
                ) = gramtools.load_gramtools_vcf_and_allele_coverage_files(
                    build_vcf, quasimap_dir
                )
                logging.info("Finished loading gramtools files")

                if self.clean:
                    shutil.rmtree(quasimap_dir)

                vcf_prefix = os.path.join(
                    self.split_output_dir,
                    f"split.{split_file.file_number}.out",
                )
                split_vcf_out = f"{vcf_prefix}.vcf"
                unfiltered_vcf_out = (
                    f"{vcf_prefix}.debug.calls_with_zero_cov_alleles.vcf"
                )

                gramtools.write_vcf_annotated_using_coverage_from_gramtools(
                    self.mean_depth,
                    self.variance_depth,
                    vcf_records,
                    allele_coverage,
                    allele_groups,
                    self.read_error_rate,
                    unfiltered_vcf_out,
                    sample_name=self.sample_name,
                    filtered_outfile=split_vcf_out,
                    ref_seq_lengths=self.ref_seq_lengths,
                    call_hets=self.call_hets,
                )

                split_vcf_outfiles[ref_name].append(split_vcf_out)
                split_vcf_outfiles_unfiltered[ref_name].append(unfiltered_vcf_out)

        # We now have minos run on each split VCF. Merge into one VCF, then can
        # add gt conf and gcp to the merged VCF.
        logging.info("Merging VCF files into one output file " + self.final_vcf)
        chunker.merge_files(split_vcf_outfiles, self.final_vcf)
        chunker.merge_files(split_vcf_outfiles_unfiltered, self.unfiltered_vcf_file)
        self.run_gt_conf()

        if self.clean:
            logging.info("Deleting temp split VCF files")
            for d in split_vcf_outfiles, split_vcf_outfiles_unfiltered:
                for file_list in d.values():
                    for filename in file_list:
                        os.unlink(filename)
            if self.use_unmapped_reads:
                os.unlink(unmapped_reads_file)

    def run_adjudicate(
        self, build_dir, quasimap_dir, vcf, reads_files, final_vcf, debug_vcf
    ):
        build_report, quasimap_report = gramtools.run_gramtools(
            build_dir,
            quasimap_dir,
            vcf,
            self.ref_fasta,
            reads_files,
            kmer_size=self.gramtools_kmer_size,
        )

        build_vcf = os.path.join(build_dir, "build.vcf")

        logging.info("Loading gramtools quasimap output files " + quasimap_dir)
        (
            self.mean_depth,
            self.variance_depth,
            vcf_header,
            vcf_records,
            allele_coverage,
            allele_groups,
        ) = gramtools.load_gramtools_vcf_and_allele_coverage_files(
            build_vcf, quasimap_dir
        )

        logging.info("Finished loading gramtools files")

        if self.clean:
            os.rename(
                os.path.join(quasimap_dir, "quasimap_outputs", "quasimap_report.json"),
                self.gramtools_quasimap_json,
            )
            shutil.rmtree(quasimap_dir)

            if not self.user_supplied_gramtools_build_dir:
                os.rename(
                    os.path.join(build_dir, "build_report.json"),
                    self.gramtools_build_json,
                )
                shutil.rmtree(build_dir)

        gramtools.write_vcf_annotated_using_coverage_from_gramtools(
            self.mean_depth,
            self.variance_depth,
            vcf_records,
            allele_coverage,
            allele_groups,
            self.read_error_rate,
            debug_vcf,
            sample_name=self.sample_name,
            filtered_outfile=final_vcf,
            ref_seq_lengths=self.ref_seq_lengths,
            call_hets=self.call_hets,
        )

    @classmethod
    def _plot_gt_conf_hists(cls, real_file, sim_file, outfile):
        with open(real_file) as f:
            real_data = [int(x.rstrip()) for x in f]
        with open(sim_file) as f:
            sim_data = [int(x.rstrip()) for x in f]

        real_color = "#7fc97f"
        sim_color = "#beaed4"
        real_patch = mpatches.Patch(color=real_color, label="Real")
        sim_patch = mpatches.Patch(color=sim_color, label="Simulated")
        fig, ax = plt.subplots()
        sns.distplot(real_data, color=real_color, ax=ax)
        sns.distplot(sim_data, color=sim_color, ax=ax)
        plt.legend(handles=[real_patch, sim_patch])
        ax.set_xlabel("GT_CONF")
        plt.savefig(outfile)
        plt.clf()
        plt.close("all")

    def run_gt_conf(self):
        if self.mean_depth > 0:
            simulations = genotype_confidence_simulator.GenotypeConfidenceSimulator(
                self.mean_depth,
                self.variance_depth,
                self.read_error_rate,
                allele_length=1,
                iterations=self.genotype_simulation_iterations,
                call_hets=self.call_hets,
            )
            simulations.run_simulations(conf_scores_file=self.sim_conf_scores_file)
        else:
            simulations = None

        logging.info(
            f"Adding GT_CONF_PERCENTILE to final VCF file {self.final_vcf} & its debug counterpart, "
            f"using mean depth {self.mean_depth}, variance depth {self.variance_depth}, error rate {self.read_error_rate}, "
            f"and {self.genotype_simulation_iterations} simulation iterations"
        )

        for f in [self.unfiltered_vcf_file, self.final_vcf]:
            if self.debug and f == self.final_vcf:
                scores_file = self.real_conf_scores_file
            else:
                scores_file = None
            Adjudicator._add_gt_conf_percentile_and_filters_to_vcf_file(
                f,
                simulations,
                min_dp=self.filter_min_dp,
                min_gcp=self.filter_min_gcp,
                max_dp=self.filter_max_dp,
                min_frs=self.filter_min_frs,
                conf_scores_file=scores_file,
            )
            if scores_file is not None:
                Adjudicator._plot_gt_conf_hists(
                    scores_file,
                    self.sim_conf_scores_file,
                    self.genotype_hist_pdf,
                )

    @classmethod
    def _add_gt_conf_percentile_and_filters_to_vcf_file(
        cls,
        vcf_file,
        geno_simulations,
        min_dp=0,
        min_gcp=5,
        max_dp=3.0,
        min_frs=0.9,
        conf_scores_file=None,
    ):
        """Overwrites vcf_file, with new version that has GT_CONF_PERCENTILE added,
        and filter for DP, GT_CONF_PERCENTILE, and FRS"""
        if conf_scores_file is not None:
            real_conf_scores = []

        vcf_header, vcf_lines = vcf_file_read.vcf_file_to_list(vcf_file)
        found_GT_CONF = False
        mean_depth = None
        variance = None
        for i, line in enumerate(vcf_header):
            if line.startswith("##FORMAT=<ID=GT_CONF"):
                found_GT_CONF = True
            elif line.startswith("##minosMeanReadDepth="):
                mean_depth = float(line.rstrip().split("=")[-1])
            elif line.startswith("##minosReadDepthVariance="):
                variance = float(line.rstrip().split("=")[-1])
        if not found_GT_CONF:
            raise Exception(
                f"No GT_CONF description found in header of VCF file {vcf_file}. Cannot continue"
            )
        if mean_depth is None or variance is None:
            raise Exception(
                f"minosMeanReadDepth and/or minosReadDepthVariance not found in header of VCF file {vcf_file}. Cannot continue"
            )

        max_dp_cutoff = mean_depth + max_dp * math.sqrt(variance)
        vcf_header[-1:-1] = [
            '##FORMAT=<ID=GT_CONF_PERCENTILE,Number=1,Type=Float,Description="Percentile of GT_CONF">',
            f'##FILTER=<ID=MIN_FRS,Description="Minimum FRS of {min_frs}">',
            f'##FILTER=<ID=MIN_DP,Description="Minimum DP of {min_dp}">',
            f'##FILTER=<ID=MIN_GCP,Description="Minimum GT_CONF_PERCENTILE of {min_gcp}">',
            f'##FILTER=<ID=MAX_DP,Description="Maximum DP of {max_dp_cutoff} (= {max_dp} standard deviations from the mean read depth {mean_depth})">',
        ]

        with open(vcf_file, "w") as f:
            print(*vcf_header, sep="\n", file=f)

            for vcf_record in vcf_lines:
                vcf_record.FILTER = set()

                if "GT" in vcf_record.FORMAT and "GT_CONF" in vcf_record.FORMAT:
                    if "." not in vcf_record.FORMAT["GT"]:
                        conf = int(round(float(vcf_record.FORMAT["GT_CONF"])))
                        vcf_record.set_format_key_value(
                            "GT_CONF_PERCENTILE",
                            str(geno_simulations.get_percentile(conf)),
                        )
                        if (
                            "DP" in vcf_record.FORMAT
                            and float(vcf_record.FORMAT["DP"]) < min_dp
                        ):
                            vcf_record.FILTER.add("MIN_DP")
                        if float(vcf_record.FORMAT["GT_CONF_PERCENTILE"]) < min_gcp:
                            vcf_record.FILTER.add("MIN_GCP")
                        if float(vcf_record.FORMAT["DP"]) > max_dp_cutoff:
                            vcf_record.FILTER.add("MAX_DP")
                        if (
                            "FRS" in vcf_record.FORMAT
                            and float(vcf_record.FORMAT["FRS"]) < min_frs
                        ):
                            vcf_record.FILTER.add("MIN_FRS")
                        if len(vcf_record.FILTER) == 0:
                            vcf_record.FILTER.add("PASS")

                        if conf_scores_file is not None:
                            real_conf_scores.append(conf)
                    else:
                        # Add a default null percentile
                        vcf_record.set_format_key_value("GT_CONF_PERCENTILE", "0.0")

                print(vcf_record, file=f)

        if conf_scores_file is not None:
            with open(conf_scores_file, "w") as f:
                print(*real_conf_scores, sep="\n", file=f)
