import json
import logging
import os
import shutil
import statistics
import sys

from cluster_vcf_records import vcf_clusterer, vcf_file_read
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

    # The structures below are designed to cope with chunking
    mean_depths = []
    variance_depths = []

    def __init__(
        self,
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
        filter_min_dp=0,
        filter_min_gcp=5,
        filter_min_frs=0.9,
        call_hets=False,
        debug=False,
    ):
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.reads_files = [os.path.abspath(x) for x in reads_files]
        self.vcf_files = [os.path.abspath(x) for x in vcf_files]
        self.overwrite_outdir = overwrite_outdir
        self.max_alleles_per_cluster = max_alleles_per_cluster
        self.sample_name = sample_name
        self.outdir = os.path.abspath(outdir)
        self.split_output_dir = os.path.join(self.outdir, "split.out")
        self.log_file = os.path.join(self.outdir, "log.txt")
        self.clustered_vcf = os.path.join(self.outdir, "gramtools.in.vcf")
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
        self.read_error_rate = read_error_rate
        self.max_read_length = max_read_length
        self.variants_per_split = variants_per_split
        self.alleles_per_split = alleles_per_split
        self.total_splits = total_splits

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
        self.filter_min_frs = filter_min_frs
        self.call_hets = call_hets
        self.debug = debug
        self.ref_seq_lengths = {
            x.id.split()[0]: len(x)
            for x in pyfastaq.sequences.file_reader(self.ref_fasta)
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
        Adjudicator.mean_depths = []
        Adjudicator.variance_depths = []
        self.build_output_dir()

        fh = logging.FileHandler(self.log_file, mode="w")
        log = logging.getLogger()
        formatter = logging.Formatter(
            "[minos %(asctime)s %(levelname)s] %(message)s", datefmt="%d-%m-%Y %H:%M:%S"
        )
        fh.setFormatter(formatter)
        log.addHandler(fh)
        logging.info("Command run: " + " ".join(sys.argv))
        dependencies.check_and_report_dependencies(programs=["gramtools"])
        logging.info("Dependencies look OK")

        if self.read_error_rate is None or self.max_read_length is None:
            logging.info(
                "One or both of read_error_rate and max_read_length not known. Estimate from first 10,000 reads..."
            )
            (
                estimated_read_length,
                estimated_read_error_rate,
            ) = utils.estimate_max_read_length_and_read_error_rate_from_qual_scores(
                self.reads_files[0]
            )
            logging.info(
                "Estimated max_read_length="
                + str(estimated_read_length)
                + " and read_error_rate="
                + str(estimated_read_error_rate)
            )

            self.read_error_rate = (
                estimated_read_error_rate
                if self.read_error_rate is None
                else self.read_error_rate
            )
            self.max_read_length = (
                estimated_read_length
                if self.max_read_length is None
                else self.max_read_length
            )
        logging.info(
            "Using max_read_length="
            + str(self.max_read_length)
            + " and read_error_rate="
            + str(self.read_error_rate)
        )

        if self.user_supplied_gramtools_build_dir:
            logging.info(
                "User supplied gramtools build dir. Assuming VCF already clustered, so skipping clustering"
            )
            assert len(self.vcf_files) == 1
            self.clustered_vcf = self.vcf_files[0]
        else:
            logging.info(
                "Clustering VCF file(s), to make one VCF input file for gramtools"
            )
            clusterer = vcf_clusterer.VcfClusterer(
                self.vcf_files,
                self.ref_fasta,
                self.clustered_vcf,
                cluster_boundary_size=0,
                max_alleles_per_cluster=self.max_alleles_per_cluster,
            )
            clusterer.run()

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

    def _run_gramtools_with_split_vcf(self):
        logging.info("Splitting VCF files into chunks (if not already done)")
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

        split_vcf_outfiles = {}
        split_vcf_outfiles_unfiltered = {}

        for ref_name, split_file_list in chunker.vcf_split_files.items():
            split_vcf_outfiles[ref_name] = []
            split_vcf_outfiles_unfiltered[ref_name] = []
            for split_file in split_file_list:
                logging.info(
                    "===== Start analysing variants in VCF split file "
                    + split_file.filename
                    + " ====="
                )
                split_reads_file = os.path.join(
                    self.split_output_dir,
                    "split." + str(split_file.file_number) + ".reads.bam",
                )
                bam_read_extract.get_region(
                    self.reads_files[0],
                    split_file.chrom,
                    split_file.chrom_start,
                    split_file.chrom_end,
                    split_reads_file,
                )

                gramtools_quasimap_dir = os.path.join(
                    self.split_output_dir,
                    "split." + str(split_file.file_number) + ".gramtools.quasimap",
                )
                if self.use_unmapped_reads:
                    reads_files = [unmapped_reads_file, split_reads_file]
                else:
                    reads_files = [split_reads_file]

                split_vcf_out = os.path.join(
                    self.split_output_dir,
                    "split." + str(split_file.file_number) + ".out.vcf",
                )
                unfiltered_vcf_out = os.path.join(
                    self.split_output_dir,
                    "split."
                    + str(split_file.file_number)
                    + ".out.debug.calls_with_zero_cov_alleles.vcf",
                )

                self.run_adjudicate(
                    split_file.gramtools_build_dir,
                    gramtools_quasimap_dir,
                    split_file.filename,
                    reads_files,
                    split_vcf_out,
                    unfiltered_vcf_out,
                )

                if self.clean:
                    os.unlink(split_reads_file)
                    if not self.user_supplied_gramtools_build_dir:
                        os.unlink(split_file.filename)

                split_vcf_outfiles[ref_name].append(split_vcf_out)
                split_vcf_outfiles_unfiltered[ref_name].append(unfiltered_vcf_out)

                logging.info(
                    "===== Finish analysing variants in VCF split file "
                    + split_file.filename
                    + " ====="
                )

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
            self.max_read_length,
            kmer_size=self.gramtools_kmer_size,
        )

        build_vcf = os.path.join(build_dir, "build.vcf")

        logging.info("Loading gramtools quasimap output files " + quasimap_dir)
        (
            mean_depth,
            variance_depth,
            vcf_header,
            vcf_records,
            allele_coverage,
            allele_groups,
        ) = gramtools.load_gramtools_vcf_and_allele_coverage_files(
            build_vcf, quasimap_dir
        )
        Adjudicator.mean_depths.append(mean_depth)
        Adjudicator.variance_depths.append(variance_depth)

        logging.info("Finished loading gramtools files")

        if self.clean:
            os.rename(
                os.path.join(quasimap_dir, "quasimap_outputs", "quasimap_report.json"),
                quasimap_dir + ".report.json",
            )
            shutil.rmtree(quasimap_dir)

            if not self.user_supplied_gramtools_build_dir:
                os.rename(
                    os.path.join(build_dir, "build_report.json"),
                    build_dir + ".report.json",
                )
                shutil.rmtree(build_dir)

        if self.sample_name is None:
            sample_name = vcf_file_read.get_sample_name_from_vcf_header_lines(
                vcf_header
            )
        else:
            sample_name = self.sample_name

        gramtools.write_vcf_annotated_using_coverage_from_gramtools(
            mean_depth,
            variance_depth,
            vcf_records,
            allele_coverage,
            allele_groups,
            self.read_error_rate,
            debug_vcf,
            sample_name=sample_name,
            max_read_length=self.max_read_length,
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
        mean_depth = statistics.mean(Adjudicator.mean_depths)
        variance_depth = statistics.mean(Adjudicator.variance_depths)

        if mean_depth > 0:
            simulations = genotype_confidence_simulator.GenotypeConfidenceSimulator(
                mean_depth,
                variance_depth,
                self.read_error_rate,
                allele_length=1,
                iterations=self.genotype_simulation_iterations,
                call_hets=self.call_hets,
            )
            simulations.run_simulations(conf_scores_file=self.sim_conf_scores_file)
        else:
            simulations = None

        logging.info(
            f"Adding GT_CONF_PERCENTLE to final VCF file {self.final_vcf} & its debug counterpart, "
            f"using mean depth {mean_depth}, variance depth {variance_depth}, error rate {self.read_error_rate}, "
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
                min_frs=self.filter_min_frs,
                conf_scores_file=scores_file,
            )
            if scores_file is not None:
                Adjudicator._plot_gt_conf_hists(
                    scores_file, self.sim_conf_scores_file, self.genotype_hist_pdf,
                )

    @classmethod
    def _add_gt_conf_percentile_and_filters_to_vcf_file(
        cls,
        vcf_file,
        geno_simulations,
        min_dp=0,
        min_gcp=5,
        min_frs=0.9,
        conf_scores_file=None,
    ):
        """Overwrites vcf_file, with new version that has GT_CONF_PERCENTILE added,
        and filter for DP, GT_CONF_PERCENTILE, and FRS"""
        if conf_scores_file is not None:
            real_conf_scores = []

        vcf_header, vcf_lines = vcf_file_read.vcf_file_to_list(vcf_file)
        for i, line in enumerate(vcf_header):
            if line.startswith("##FORMAT=<ID=GT_CONF"):
                break
        else:
            raise Exception(
                f"No GT_CONF description found in header of VCF file {vcf_file}. Cannot continue"
            )

        vcf_header[i + 1 : i + 1] = [
            '##FORMAT=<ID=GT_CONF_PERCENTILE,Number=1,Type=Float,Description="Percentile of GT_CONF">',
            f'##FILTER=<ID=MIN_FRS,Description="Minimum FRS of {min_frs}">',
            f'##FILTER=<ID=MIN_DP,Description="Minimum DP of {min_dp}">',
            f'##FILTER=<ID=MIN_GCP,Description="Minimum GT_CONF_PERCENTILE of {min_gcp}">',
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
