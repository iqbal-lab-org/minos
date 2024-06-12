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
        default=2,
    )
    subparser_adjudicate.add_argument(
        "--filter_min_frs",
        type=float,
        help="Minimum value to be used for MIN_FRS filter in output VCF file [%(default)s]",
        default=0.9,
    )
    subparser_adjudicate.add_argument(
        "--filter_min_gcp",
        type=float,
        help="Minimum genotype confidence percentile to be used for MIN_GCP filter in output VCF file [%(default)s]",
        default=0.5,
    )
    subparser_adjudicate.add_argument(
        "--filter_max_dp",
        type=float,
        help="Determines maximum value to be used for MAX_DP filter in output VCF file. Using a value of 'x' here sets the cutoff to be more than x standard deviations from the mean read depth [%(default)s]",
        default=3.0,
    )
    subparser_adjudicate.add_argument(
        "--include_het_calls",
        action="store_true",
        #help="NOT IMPLEMENTED! Consider heterozygous calls when genotyping each allele, instead of only making homozygous calls.",
        help=argparse.SUPPRESS,
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

    # ------------------- get_test_data -------------------------------------------
    subparser_get_test_data = subparsers.add_parser(
        "get_test_data",
        help="Get directory of test data files, for testing minos runs ok",
        usage="minos get_test_data <outdir>",
        description="Get directory of test data files, for testing minos runs ok",
    )

    subparser_get_test_data.add_argument(
        "outdir",
        help="Name of output directory (is created, must not already exist)",
    )
    subparser_get_test_data.set_defaults(
        func=minos.tasks.get_test_data.run
    )

    # ------------------ vcf_merge ------------------------------------------------
    subparser_vcf_merge = subparsers.add_parser(
        "vcf_merge",
        help="Merge VCF files; makes input to vcf_cluster",
        usage="minos vcf_merge [options] <vcf_fofn> <ref_fasta> <outdir>",
        description="Gathers variants from VCF files, storing variants and sample info in new directory than can be used as input to vcf_cluster",
    )

    subparser_vcf_merge.add_argument(
        "vcf_fofn", help="File of VCF filenames",
    )
    subparser_vcf_merge.add_argument(
        "ref_fasta",
        help="FASTA file of reference, corresponding to all input VCF files",
    )
    subparser_vcf_merge.add_argument(
        "outdir", help="Output directory",
    )
    subparser_vcf_merge.add_argument(
        "--temp_dir",
        help="Temporary directory used for processing input VCF files. If does not exist, will be created (then destroyed at the end of the run). If this option is not used, then `tmp` is used inside the output directory",
        metavar="FILENAME",
    )
    subparser_vcf_merge.add_argument(
        "--cpus",
        type=int,
        help="Number of CPUs to use. This many input VCF files will be read in parallel [%(default)s]",
        default=1,
        metavar="INT",
    )
    subparser_vcf_merge.add_argument(
        "--sample_limit",
        type=int,
        help=argparse.SUPPRESS,  # hidden option used when running tests
    )
    subparser_vcf_merge.add_argument(
        "--mem_limit",
        type=float,
        help="Sample information is stored in chunks. This limits the RAM in GB per chunk (total ram usage of the program will be higher), writing each chunk to disk when the limit is reached [%(default)s]",
        default=2.0,
        metavar="FLOAT",
    )
    subparser_vcf_merge.add_argument(
        "--force",
        action="store_true",
        help="Replace output directory if it already exists",
    )
    subparser_vcf_merge.set_defaults(func=minos.tasks.vcf_merge.run)

    # ------------------ vcf_cluster ----------------------------------------------
    subparser_vcf_cluster = subparsers.add_parser(
        "vcf_cluster",
        help="Make clustered VCF file, using output of vcf_merge",
        usage="minos vcf_cluster [options] <ref_fasta> <merge_dir> <outprefix>",
        description="Make clustered VCF file, using output of vcf_merge",
    )

    subparser_vcf_cluster.add_argument(
        "ref_fasta",
        help="FASTA file of reference, corresponding to all input VCF files",
    )
    subparser_vcf_cluster.add_argument(
        "merge_dir", help="Input directory, made by vcf_merge",
    )
    subparser_vcf_cluster.add_argument(
        "outprefix", help="Prefix of output files",
    )
    subparser_vcf_cluster.add_argument(
        "--cpus",
        type=int,
        help="Number of CPUs to use. If >1, splits genome into chunks and processes in parallel [%(default)s]",
        default=1,
        metavar="INT",
    )
    subparser_vcf_cluster.add_argument(
        "--max_ref_len",
        type=int,
        help="Maximum allowed length of reference allele. Any longer than this are excluded [%(default)s]",
        default=50,
        metavar="INT",
    )
    subparser_vcf_cluster.add_argument(
        "--max_alleles",
        type=int,
        help="If a variant site has more than this many alleles (after generating all combinations of SNPs/indels), then instead only use the combinations of variants that are seen in the samples [%(default)s]",
        default=500,
        metavar="INT",
    )
    subparser_vcf_cluster.set_defaults(func=minos.tasks.vcf_cluster.run)

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
