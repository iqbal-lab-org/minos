#!/usr/bin/env python3

import argparse
import logging
import minos


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="minos",
        usage="minos <command> <options>",
        description="minos: variant call adjudication",
    )

    parser.add_argument("--version", action="version", version=minos.__version__)
    parser.add_argument(
        "--debug",
        help="More verbose logging, and less file cleaning",
        action="store_true",
    )

    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # ------------------------ adjudicate -----------------------------------------
    subparser_adjudicate = subparsers.add_parser(
        "adjudicate",
        help="Choose correct variants from VCF files",
        usage="minos adjudicate [options] <--reads reads_file> <outdir> <ref_fasta> <vcf_in_1> [vcf_in_2 ...]",
        description="Choose correct variants from VCF files by remapping to graph",
        epilog="IMPORTANT: the --reads option is required",
    )

    subparser_adjudicate.add_argument(
        "--reads",
        action="append",
        required=True,
        help="REQUIRED. Reads file. Can be any format compatible with htslib. Use this option more than once for >1 reads files. If splitting (with one of --total_splits,--variants_per_split,--alleles_per_split), must provide one sorted indexed BAM file of reads",
        metavar="FILENAME",
    )
    subparser_adjudicate.add_argument(
        "--gramtools_build_dir",
        help="Gramtools build directory corresponding to input VCF file. If used, assumes VCF is clustered, skips clustering and gramtools build stages and uses the provided directory instead",
        metavar="DIRNAME",
    )
    subparser_adjudicate.add_argument(
        "--gramtools_kmer_size",
        type=int,
        help="This number is used with gramtools build --kmer-size. Ignored if --gramtools_build_dir used.  [%(default)s]",
        default=10,
        metavar="INT",
    )
    subparser_adjudicate.add_argument(
        "--max_read_length",
        type=int,
        help="Maximum read length, this is used by gramtools. If not given, estimated by taking longest of first 10,000 reads",
        metavar="INT",
    )
    subparser_adjudicate.add_argument(
        "--read_error_rate",
        type=float,
        help="Read error rate between 0 and 1 [%(default)s]",
        metavar="FLOAT",
        default=0.002,
    )
    subparser_adjudicate.add_argument(
        "--max_alleles_per_cluster",
        type=int,
        help="Maximum allowed alleles in one cluster. If there are too many alleles then combinations of SNPs are not generated [%(default)s]",
        metavar="INT",
        default=5000,
    )
    subparser_adjudicate.add_argument(
        "--force", action="store_true", help="Replace outdir, if it already exists"
    )
    subparser_adjudicate.add_argument(
        "--sample_name",
        help="Sample name to put in final VCF output file. Default is to use first sample name found in input VCF file(s)",
        metavar="STRING",
    )
    subparser_adjudicate.add_argument(
        "--total_splits",
        type=int,
        help="Split VCF, aiming for this many chunks with the same number of variants in each chunk. Increases run time, but saves RAM (see also --variants_per_split and --alleles_per_split). If used, then reads must be in one sorted indexed BAM file",
        metavar="INT",
    )
    subparser_adjudicate.add_argument(
        "--variants_per_split",
        type=int,
        help="Split VCF, aiming for this many variants in each split. Takes precedence over --total_splits. Increases run time, but saves RAM. If used, then reads must be in one sorted indexed BAM file",
        metavar="INT",
    )
    subparser_adjudicate.add_argument(
        "--alleles_per_split",
        type=int,
        help="Split VCF, aiming for this many alleles in each split. Takes precedence over --total_splits. Increases run time, but saves RAM. If used, then reads must be in one sorted indexed BAM file",
        metavar="INT",
    )
    subparser_adjudicate.add_argument(
        "--use_unmapped_reads",
        action="store_true",
        help="When using splitting (with --total_splits, --variants_per_split, or --alleles_per_split), use the unmapped reads with each split. Default is to ignore them. Using this option may add huge increase to run time, with little benefit to variant call accuracy",
    )
    subparser_adjudicate.add_argument(
        "--filter_min_dp",
        type=int,
        help="Minimum depth to be used for MIN_DP filter in output VCF file [%(default)s]",
        default=0,
    )
    subparser_adjudicate.add_argument(
        "--filter_min_frs",
        type=int,
        help="Minimum value to be used for MIN_FRS filter in output VCF file [%(default)s]",
        default=0.9,
    )
    subparser_adjudicate.add_argument(
        "--filter_min_gcp",
        type=float,
        help="Minimum genotype confidence percentile to be used for MIN_GCP filter in output VCF file [%(default)s]",
        default=5.0,
    )
    subparser_adjudicate.add_argument(
        "--include_het_calls",
        action="store_true",
        help="Consider heterozygous calls when genotyping each allele, instead of only making homozygous calls.",
    )
    subparser_adjudicate.add_argument("outdir", help="Name of output directory")
    subparser_adjudicate.add_argument(
        "ref_fasta", help="Reference FASTA filename (must match VCF file(s))"
    )
    subparser_adjudicate.add_argument(
        "vcf_files",
        nargs="+",
        help="VCF filename(s) to be merged. Must provide at least one filename.",
    )
    subparser_adjudicate.set_defaults(func=minos.tasks.adjudicate.run)

    # ------------------------ cluster_vcfs ---------------------------------------
    subparser_cluster_vcfs = subparsers.add_parser(
        "cluster_vcfs",
        help="Cluster one or more VCF files",
        usage="minos cluster_vcfs [options] <ref_fasta> <vcf_in_1> [vcf_in_2 ...]",
        description="Clusters one of more VCF files, outputting a single VCF file where overlapping variants are merged into one record",
    )

    subparser_cluster_vcfs.add_argument(
        "--max_alleles_per_cluster",
        type=int,
        help="Maximum allowed alleles in one cluster. If there are too many alleles then combinations of SNPs are not generated [%(default)s]",
        metavar="INT",
        default=5000,
    )
    subparser_cluster_vcfs.add_argument(
        "--max_var_dist",
        type=int,
        help="Variants up to this distance apart are merged [%(default)s]",
        default=1,
        metavar="INT",
    )
    subparser_cluster_vcfs.add_argument(
        "-o",
        "--outfile",
        help="Name of output VCF file (default is stdout)",
        default="-",
        metavar="FILENAME",
    )
    subparser_cluster_vcfs.add_argument(
        "ref_fasta", help="FASTA filename of reference genome that matches the VCF file"
    )
    subparser_cluster_vcfs.add_argument(
        "vcf_files",
        nargs="+",
        help="VCF filename(s) to be merged. Must provide at least one filename.",
    )
    subparser_cluster_vcfs.set_defaults(func=minos.tasks.cluster_vcfs.run)

    # ----------------- make_split_gramtools_build --------------------------------
    subparser_make_split_gramtools_build = subparsers.add_parser(
        "make_split_gramtools_build",
        help="Split VCF/ref into pieces and run gramtools build",
        usage="minos make_split_gramtools_build [options] <outdir> <VCF file> <ref FASTA>",
        description='Split VCF/ref into pieces and run gramtools build. The resulting directory can be used with "minos adjudicate --gramtools_build_dir", so that the same build can be used on multiple samples without generating it during each run of minos',
    )

    subparser_make_split_gramtools_build.add_argument(
        "outdir", help="Name of output directory"
    )
    subparser_make_split_gramtools_build.add_argument(
        "vcf_file", help="Input VCF file to be split"
    )
    subparser_make_split_gramtools_build.add_argument(
        "ref_fasta", help="Reference FASTA file"
    )
    subparser_make_split_gramtools_build.add_argument(
        "--total_splits",
        type=int,
        help="Split VCF, aiming for this many chunks with the same number of variants in each chunk (see also --variants_per_split and --alleles_per_split) [%(default)s]",
        metavar="INT",
        default=100,
    )
    subparser_make_split_gramtools_build.add_argument(
        "--variants_per_split",
        type=int,
        help="Split VCF, aiming for this many variants in each split. If used, --total_splits is ignored",
        metavar="INT",
    )
    subparser_make_split_gramtools_build.add_argument(
        "--alleles_per_split",
        type=int,
        help="Split VCF, aiming for this many alleles in each split. If used, --total_splits is ignored",
        metavar="INT",
    )
    subparser_make_split_gramtools_build.add_argument(
        "--max_read_length",
        type=int,
        help="This number is used with gramtools build --max-read-length [%(default)s]",
        default=200,
    )
    subparser_make_split_gramtools_build.add_argument(
        "--gramtools_kmer_size",
        type=int,
        help="This number is used with gramtools build --kmer-size [%(default)s]",
        default=10,
        metavar="INT",
    )
    subparser_make_split_gramtools_build.add_argument(
        "--threads",
        type=int,
        help="Number of gramtools builds to run in parallel [%(default)s]",
        default=1,
        metavar="INT",
    )
    subparser_make_split_gramtools_build.set_defaults(
        func=minos.tasks.make_split_gramtools_build.run
    )

    # ----------------- multi_sample_pipeline -------------------------------------

    subparser_multi_sample_pipeline = subparsers.add_parser(
        "multi_sample_pipeline",
        help="Re-call across mulitple samples",
        usage="minos multi_sample_pipeline [options] <ref_fasta> <data_tsv> <outdir>",
        description="Re-calls variants across multiple samples. Re-genotypes every sample at every position where there is a variant in the input. Assumes input VCF files made by minos.",
    )

    subparser_multi_sample_pipeline.add_argument(
        "ref_fasta", help="FASTA file of reference, corresponding to all input files"
    )
    subparser_multi_sample_pipeline.add_argument(
        "data_tsv",
        help="TSV file of input files. One sample per line. column1=VCF filename. column2=reads file. Can optionally have more columns for any sample that has more than one reads file.",
    )
    subparser_multi_sample_pipeline.add_argument(
        "outdir", help="Name of output directory"
    )
    subparser_multi_sample_pipeline.add_argument(
        "--max_alleles_per_cluster",
        type=int,
        help="Maximum allowed alleles in one cluster. If there are too many alleles then combinations of SNPs are not generated [%(default)s]",
        metavar="INT",
        default=5000,
    )
    subparser_multi_sample_pipeline.add_argument(
        "--min_large_ref_length",
        type=int,
        help="Variants with REF string at least this long are counted as large. These are currently removed from the analysis [%(default)s]",
        default=50,
        metavar="INT",
    )
    subparser_multi_sample_pipeline.add_argument(
        "--max_read_length",
        type=int,
        help="Maximum read length (--max-read-length option of gramtools build). If not given, determined from header of VCF files (if made by minos). The default of zero means determine from VCF files.",
        default=0,
        metavar="INT",
    )
    subparser_multi_sample_pipeline.add_argument(
        "--gramtools_kmer_size",
        type=int,
        help="This number is used with gramtools build --kmer-size [%(default)s]",
        default=10,
        metavar="INT",
    )
    subparser_multi_sample_pipeline.add_argument(
        "--gramtools_build_threads",
        type=int,
        help="Number of gramtools builds to run in parallel (only used if --total_splits,--variants_per_split, or --alleles_per_split used) [%(default)s]",
        default=1,
        metavar="INT",
    )
    subparser_multi_sample_pipeline.add_argument(
        "--total_splits",
        type=int,
        help="Split VCF, aiming for this many chunks with the same number of variants in each chunk. Increases run time, but saves RAM (see also --variants_per_split and --alleles_per_split). If used, then reads for each sample must be in one sorted indexed BAM file",
        metavar="INT",
    )
    subparser_multi_sample_pipeline.add_argument(
        "--variants_per_split",
        type=int,
        help="Split VCF, aiming for this many variants in each split. Takes precedence over --total_splits. Increases run time, but saves RAM. If used, then reads for each sample must be in one sorted indexed BAM file",
        metavar="INT",
    )
    subparser_multi_sample_pipeline.add_argument(
        "--alleles_per_split",
        type=int,
        help="Split VCF, aiming for this many alleles in each split. Takes precedence over --total_splits. Increases run time, but saves RAM. If used, then reads for each sample must be in one sorted indexed BAM file",
        metavar="INT",
    )
    subparser_multi_sample_pipeline.add_argument(
        "--use_unmapped_reads",
        action="store_true",
        help="When using splitting (with --total_splits, --variants_per_split, or --alleles_per_split), use the unmapped reads with each split. Default is to ignore them. Using this option may add huge increase to run time, with little benefit to variant call accuracy",
    )
    subparser_multi_sample_pipeline.add_argument(
        "--nextflow_work_dir",
        help="Nextflow work directory (-work-dir option of nextflow). If not given subdirectory of outdir is used",
        metavar="DIRECTORY_NAME",
    )
    subparser_multi_sample_pipeline.add_argument(
        "--nextflow_config_file",
        help="Nextflow config file (-c option of nextflow). If not given, nextflow's defaults are used",
        metavar="FILENAME",
    )
    subparser_multi_sample_pipeline.add_argument(
        "--force", action="store_true", help="Overwrite the output directory"
    )
    subparser_multi_sample_pipeline.add_argument(
        "--no_run",
        action="store_true",
        help="Make all input files for nextflow, but do not run nextflow",
    )
    subparser_multi_sample_pipeline.add_argument(
        "--no_clean",
        action="store_true",
        help="Do not clean up temporary nextflow files",
    )
    subparser_multi_sample_pipeline.add_argument(
        "--nf_ram_cluster_small_vars",
        type=float,
        help="Nextflow RAM limit when clustering small variants [%(default)s]",
        metavar="FLOAT",
        default=2,
    )
    subparser_multi_sample_pipeline.add_argument(
        "--nf_ram_gramtools_build_small",
        type=float,
        help="Nextflow RAM limit when running gramtools build on small variants [%(default)s]",
        metavar="FLOAT",
        default=12,
    )
    subparser_multi_sample_pipeline.add_argument(
        "--nf_ram_minos_small_vars",
        type=float,
        help="Nextflow RAM limit when running minos on small variants [%(default)s]",
        metavar="FLOAT",
        default=5,
    )
    subparser_multi_sample_pipeline.add_argument(
        "--nf_ram_merge_small_vars",
        type=float,
        help="Nextflow RAM limit when merging small variant vcf files [%(default)s]",
        metavar="FLOAT",
        default=2,
    )
    subparser_multi_sample_pipeline.add_argument(
        "--testing", action="store_true", help=argparse.SUPPRESS
    )
    subparser_multi_sample_pipeline.set_defaults(
        func=minos.tasks.multi_sample_pipeline.run
    )

    # ------------------------ versions -------------------------------------------
    subparser_versions = subparsers.add_parser(
        "versions",
        help="Report versions of minos and dependencies",
        usage="minos versions",
        description="Checks all dependencies are found. Reports where they are found and their versions",
    )
    subparser_versions.set_defaults(func=minos.tasks.versions.run)

    args = parser.parse_args()

    log = logging.getLogger()
    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
