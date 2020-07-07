import copy
import datetime
import json
import logging
import os
import statistics

from cluster_vcf_records import vcf_file_read

from minos import dependencies, genotyper, utils
from minos import __version__ as minos_version


def _build_json_file_is_good(json_build_report):
    """Returns true iff looks like gramtools build_report.json
    says that gramtools build ran successfully"""
    if not os.path.exists(json_build_report):
        return False

    with open(json_build_report) as f:
        build_report = json.load(f)
        try:
            success = build_report["success"]
        except:
            return False

        assert type(success) == bool
        return success


def run_gramtools_build(outdir, vcf_file, ref_file, kmer_size=10):
    """Runs gramtools build. Makes new directory called 'outdir' for
    the output"""
    if os.path.exists(outdir):
        raise FileExistsError(
            f"Gramtools build output directory '{outdir}' already exists. Cannot continue"
        )
    os.mkdir(outdir)
    gramtools_exe = dependencies.find_binary("gramtools")
    build_command = " ".join(
        [
            gramtools_exe,
            "build",
            "--gram-directory",
            outdir,
            "--vcf",
            vcf_file,
            "--reference",
            ref_file,
            "--kmer-size",
            str(kmer_size),
            "--no-vcf-clustering",
        ]
    )
    logging.info("Running gramtools build: " + build_command)
    completed_process = utils.syscall(build_command, allow_fail=True)
    logging.info(
        "Finished running gramtools build. Return code: "
        + str(completed_process.returncode)
    )
    build_report = os.path.join(outdir, "build_report.json")
    ran_ok = (
        _build_json_file_is_good(build_report) and completed_process.returncode == 0
    )
    if not ran_ok:
        logging.info(
            "Error running gramtools build. See build report file " + build_report
        )
        raise Exception(
            f"Error running gramtools build: {build_command}\nstdout:{completed_process.stdout}\nstderr:\n{completed_process.stderr}"
        )

    logging.info("Build report file looks good from gramtools build: " + build_report)


def run_gramtools(
    build_dir, quasimap_dir, vcf_file, ref_file, reads, kmer_size=10, seed=42,
):
    """If build_dir does not exist, runs runs gramtools build and quasimap.
    Otherwise, just runs quasimap. quasimap output is in new
    directory called quasimap_dir.
    "reads" can be one filename, or a list of filenames.
    Raises Error if either of the expected json coverage
    files made by quasimap are not found."""
    gramtools_exe = dependencies.find_binary("gramtools")
    if not os.path.exists(build_dir):
        run_gramtools_build(build_dir, vcf_file, ref_file, kmer_size=kmer_size)

    if type(reads) is not list:
        assert type(reads) is str
        reads = [reads]

    quasimap_command = " ".join(
        [
            gramtools_exe,
            "quasimap",
            f"--seed {seed}",
            "--gram-dir",
            build_dir,
            "--run-dir",
            quasimap_dir,
            " ".join(["--reads " + x for x in reads]),
        ]
    )
    logging.info("Running gramtools quasimap: " + quasimap_command)
    utils.syscall(quasimap_command)
    logging.info("Finished running gramtools quasimap")

    build_report = os.path.join(build_dir, "build_report.json")
    quasimap_report = os.path.join(
        quasimap_dir, "quasimap_outputs", "quasimap_report.json"
    )
    allele_base_counts_file = os.path.join(
        quasimap_dir, "quasimap_outputs", "allele_base_coverage.json"
    )
    grouped_allele_counts_file = os.path.join(
        quasimap_dir, "quasimap_outputs", "grouped_allele_counts_coverage.json"
    )
    files_ok = True
    for filename in (
        build_report,
        quasimap_report,
        allele_base_counts_file,
        grouped_allele_counts_file,
    ):
        if not (os.path.exists(filename)):
            files_ok = False
            logging.error("gramtools file not found: " + filename)

    if not files_ok:
        error_message = "Looks like something went wrong during gramtools run. At least one output file not present. Cannot continue."
        logging.error(error_message)
        raise Exception(error_message)

    with open(build_report) as f:
        json_build_report = json.load(f)
    with open(quasimap_report) as f:
        json_quasimap_report = json.load(f)

    return json_build_report, json_quasimap_report


def load_gramtools_vcf_and_allele_coverage_files(vcf_file, quasimap_dir):
    """Loads the perl_generated_vcf file and allele_coverage files.
    Sanity checks that they agree: 1) same number of lines (excluding header
    lines in vcf) and 2) number of alts agree on each line.
    Raises error at the first time somthing wrong is found.
    Returns a list of tuples: (VcfRecord, dict of allele -> coverage)"""
    allele_base_counts_file = os.path.join(
        quasimap_dir, "quasimap_outputs", "allele_base_coverage.json"
    )
    grouped_allele_counts_file = os.path.join(
        quasimap_dir, "quasimap_outputs", "grouped_allele_counts_coverage.json"
    )
    all_allele_coverage, allele_groups = load_allele_files(
        allele_base_counts_file, grouped_allele_counts_file
    )
    vcf_header, vcf_lines = vcf_file_read.vcf_file_to_list(vcf_file)
    coverages = []

    if len(all_allele_coverage) != len(vcf_lines):
        raise Exception(
            "Number of records in VCF ("
            + str(len(vcf_lines))
            + ") does not match number output from gramtools.("
            + str(len(all_allele_coverage))
            + "). Cannot continue"
        )

    for i, (allele_combi_coverage, allele_per_base_coverage) in enumerate(
        all_allele_coverage
    ):
        if len(allele_per_base_coverage) != 1 + len(vcf_lines[i].ALT):
            raise Exception(
                "Mismatch in number of alleles for this VCF record:\n"
                + str(vcf_lines[i])
                + "\nLine number is "
                + str(i + 1)
            )

        coverages.append(sum(allele_combi_coverage.values()))

    assert len(coverages) > 0
    # Unlikely to happen edge case on real data is when coverages has length 1.
    # It happens when running test_run in adjudicator_test, with a split VCf.
    # One of the splits only has 1 record.
    if len(coverages) == 1:
        variance = 1.000
    else:
        variance = round(statistics.variance(coverages), 3)

    return (
        round(statistics.mean(coverages), 3),
        variance,
        vcf_header,
        vcf_lines,
        all_allele_coverage,
        allele_groups,
    )


def update_vcf_record_using_gramtools_allele_depths(
    vcf_record, gtyper, allele_combination_cov, allele_per_base_cov, allele_groups_dict,
):
    """allele_depths should be a dict of allele -> coverage.
    The REF allele must also be in the dict.
    So keys of dict must be equal to REF + ALTs sequences.
    This also changes all columns from QUAL onwards.
    Returns a VcfRecord the same as vcf_record, but with all zero
    coverage alleles removed, and GT and COV fixed accordingly"""
    logging.debug(f"Genotyping. allele_combination_cov={[str(allele_groups_dict[x]) + ': ' + str(allele_combination_cov[x]) for x in allele_combination_cov]}")
    logging.debug(f"Genotyping. allele_per_base_cov={allele_per_base_cov}")
    gtyper.run(allele_combination_cov, allele_per_base_cov, allele_groups_dict)
    genotype_indexes = set()

    if "." in gtyper.genotype:
        genotype = "./."
    else:
        if 0 in gtyper.genotype:
            genotype_indexes.add(0)
        for i in range(len(vcf_record.ALT)):
            if i + 1 in gtyper.genotype:
                genotype_indexes.add(i + 1)

        if len(genotype_indexes) == 1:
            genotype_index = genotype_indexes.pop()
            genotype = str(genotype_index) + "/" + str(genotype_index)
            genotype_indexes.add(genotype_index)
        else:
            genotype = "/".join([str(x) for x in sorted(list(genotype_indexes))])

    cov_values = [
        gtyper.singleton_allele_coverages.get(x, 0)
        for x in range(1 + len(vcf_record.ALT))
    ]
    cov_string = ",".join([str(x) for x in cov_values])
    vcf_record.QUAL = None
    vcf_record.INFO.clear()
    vcf_record.FILTER = set()
    vcf_record.FORMAT.clear()
    dp = sum(allele_combination_cov.values())
    vcf_record.set_format_key_value("DP", str(dp))
    dpf = 0 if dp == 0 else str(round(dp / gtyper.mean_depth, 4))
    vcf_record.set_format_key_value("DPF", str(dpf))
    vcf_record.set_format_key_value("GT", genotype)
    vcf_record.set_format_key_value("COV", cov_string)
    vcf_record.set_format_key_value("FRS", str(gtyper.genotype_frs))
    vcf_record.set_format_key_value("GT_CONF", str(gtyper.genotype_confidence))

    # Make new record where all zero coverage alleles are removed
    filtered_record = copy.deepcopy(vcf_record)
    if genotype in ["./.", "0/0"]:
        return filtered_record

    indexes_to_keep = set(
        [i for i in range(len(cov_values)) if i == 0 or cov_values[i] > 0]
    )
    indexes_to_keep.update(genotype_indexes)
    indexes_to_keep = list(indexes_to_keep)
    indexes_to_keep.sort()
    filtered_record.set_format_key_value(
        "COV", ",".join([str(cov_values[i]) for i in indexes_to_keep])
    )
    assert indexes_to_keep[0] == 0
    filtered_record.ALT = [filtered_record.ALT[i - 1] for i in indexes_to_keep[1:]]

    # The indexes of the genotype string 'n/m' are shifted because
    # we probably removed some alleles
    genotype_strings = {
        vcf_record.REF if i == 0 else vcf_record.ALT[i - 1] for i in genotype_indexes
    }
    new_genotype_indexes = set()
    if 0 in genotype_indexes:
        new_genotype_indexes.add(0)
    for i, genotype_string in enumerate(filtered_record.ALT):
        if genotype_string in genotype_strings:
            new_genotype_indexes.add(i + 1)
            if len(genotype_strings) == len(new_genotype_indexes):
                break

    new_genotype_indexes = list(new_genotype_indexes)
    if len(new_genotype_indexes) == 1:
        new_genotype_indexes.append(new_genotype_indexes[0])
    assert len(new_genotype_indexes) == 2
    filtered_record.set_format_key_value(
        "GT", "/".join([str(x) for x in new_genotype_indexes])
    )
    return filtered_record


def write_vcf_annotated_using_coverage_from_gramtools(
    mean_depth,
    depth_variance,
    vcf_records,
    all_allele_coverage,
    allele_groups,
    read_error_rate,
    outfile,
    sample_name="SAMPLE",
    filtered_outfile=None,
    ref_seq_lengths=None,
    call_hets=False,
):
    """mean_depth, vcf_records, all_allele_coverage, allele_groups should be those
    returned by load_gramtools_vcf_and_allele_coverage_files().
    Writes a new VCF that has allele counts for all the ALTs"""
    assert len(vcf_records) == len(all_allele_coverage)

    header_lines = [
        "##fileformat=VCFv4.2",
        "##source=minos, version " + minos_version,
        "##fileDate=" + str(datetime.date.today()),
        '##FORMAT=<ID=COV,Number=R,Type=Integer,Description="Number of reads on ref and alt alleles">',
        '##FORMAT=<ID=FRS,Number=1,Type=Float,Description="Fraction of reads that support the genotype call">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="total read depth from gramtools">',
        '##FORMAT=<ID=DPF,Number=1,Type=Float,Description="Depth Fraction, defined as DP divided by mean depth">',
        '##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description="Genotype confidence. Difference in log likelihood of most likely and next most likely genotype">',
        f"##minosMeanReadDepth={mean_depth}",
    ]

    if ref_seq_lengths is not None:
        for name, length in sorted(ref_seq_lengths.items()):
            header_lines.append(f"##contig=<ID={name},length={length}>")

    header_lines.append(
        "\t".join(
            [
                "#CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                sample_name,
            ]
        )
    )

    gtyper = genotyper.Genotyper(
        mean_depth, depth_variance, read_error_rate, call_hets=call_hets,
    )

    if filtered_outfile is not None:
        f_filter = open(filtered_outfile, "w")
        print(*header_lines, sep="\n", file=f_filter)

    with open(outfile, "w") as f:
        print(*header_lines, sep="\n", file=f)

        for i in range(len(vcf_records)):
            logging.debug("Genotyping: " + str(vcf_records[i]))
            filtered_record = update_vcf_record_using_gramtools_allele_depths(
                vcf_records[i],
                gtyper,
                all_allele_coverage[i][0],
                all_allele_coverage[i][1],
                allele_groups,
            )
            print(vcf_records[i], file=f)
            if filtered_outfile is not None:
                print(filtered_record, file=f_filter)

    if filtered_outfile is not None:
        f_filter.close()


def load_allele_files(allele_base_counts_file, grouped_allele_counts_file):
    """Loads the allele base counts and groupeed allele counts files
    made by gramtools qausimap."""
    with open(allele_base_counts_file) as f:
        json_base_counts_data = json.load(f)
    with open(grouped_allele_counts_file) as f:
        json_allele_counts_data = json.load(f)

    try:
        allele_base_counts = json_base_counts_data["allele_base_counts"]
    except:
        raise Exception(
            "Error in json file "
            + allele_base_counts_file
            + ". allele_base_counts not found."
        )

    try:
        site_counts = json_allele_counts_data["grouped_allele_counts"]["site_counts"]
    except:
        raise Exception(
            "Error in json file "
            + grouped_allele_counts_file
            + ". site_counts not found."
        )

    try:
        allele_groups = json_allele_counts_data["grouped_allele_counts"][
            "allele_groups"
        ]
    except:
        raise Exception(
            "Error in json file "
            + grouped_allele_counts_file
            + ". allele_groups not found."
        )

    if len(allele_base_counts) != len(site_counts):
        raise Exception(
            "Mismatch between number of records in json files "
            + allele_base_counts_file
            + " and "
            + grouped_allele_counts_file
        )

    for key, value in allele_groups.items():
        allele_groups[key] = set(value)

    return list(zip(site_counts, allele_base_counts)), allele_groups
