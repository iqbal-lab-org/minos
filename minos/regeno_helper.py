import csv
import itertools
import json
import logging
import multiprocessing
import os
import re

import pyfastaq
from cluster_vcf_records import vcf_file_read, vcf_record

from minos import utils


def vcf_has_too_many_variants(
    vcf_file, max_number_of_records, max_ref_proportion, ref_length
):
    number_of_records = 0
    total_length = 0

    with vcf_file_read.open_vcf_file_for_reading(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            record = vcf_record.VcfRecord(line)
            if "GT" in record.FORMAT:
                gt = set(re.split("[/|]", record.FORMAT["GT"]))
                gt.discard("0")
                gt.discard(".")
                if len(gt) > 0:
                    number_of_records += 1
                    total_length += len(record.REF)

            if (
                number_of_records > max_number_of_records
                or (total_length / ref_length) > max_ref_proportion
            ):
                return True

    return False


def compress_file(filenames):
    infile, outfile = filenames
    zipper = "bgzip" if infile.endswith(".vcf") else "gzip -9"
    utils.syscall(f"{zipper} -c {infile} > {outfile}")


def _sample_is_ok(sample_tuple, max_number_of_records, max_ref_proportion, ref_length):
    return not vcf_has_too_many_variants(
        sample_tuple[1], max_number_of_records, max_ref_proportion, ref_length
    )


def parse_manifest_file(
    infile,
    merge_fofn,
    adjudicate_tsv,
    ref_fasta,
    cpus=1,
    max_number_of_records=15000,
    max_ref_proportion=0.1,
):
    all_samples = []
    seq_reader = pyfastaq.sequences.file_reader(ref_fasta)
    ref_length = sum([len(x) for x in seq_reader])

    logging.info(f"Load manifest file {infile}")
    with open(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for c in ("name", "vcf", "reads"):
            if c not in reader.fieldnames:
                raise ValueError(f"Column called '{c}' not found in file {infile}")

        for d in reader:
            reads = "--reads " + " --reads ".join(d["reads"].split())
            all_samples.append((d["name"], d["vcf"], reads))

    logging.info(f"Checking VCFs for each sample in parallel using {cpus} cpu(s)")
    with multiprocessing.Pool(cpus) as pool:
        samples_ok = pool.starmap(
            _sample_is_ok,
            zip(
                all_samples,
                itertools.repeat(max_number_of_records),
                itertools.repeat(max_ref_proportion),
                itertools.repeat(ref_length),
            ),
        )
    logging.info(
        f"Finished checking VCFs for each sample. Keeping {sum(samples_ok)} of {len(samples_ok)} samples."
    )
    logging.info("Writing output files")

    with open(merge_fofn, "w") as f_merge, open(adjudicate_tsv, "w") as f_adj:
        print("name", "reads", sep="\t", file=f_adj)
        for sample_ok, sample_data in zip(samples_ok, all_samples):
            if not sample_ok:
                continue

            print(sample_data[1], file=f_merge)
            print(sample_data[0], sample_data[2], sep="\t", file=f_adj)

    logging.info("Finished writing output files")


def manifest_to_set_of_sample_names(infile):
    with open(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return {d["name"] for d in reader}


def make_per_sample_vcfs_dir(
    sample_data_tsv, root_outdir, original_manifest=None, samples_per_dir=1000, cpus=1
):
    vcf_root_out = os.path.join("VCFs")
    logs_root_out = os.path.join("Logs")
    if not os.path.exists(root_outdir):
        os.mkdir(root_outdir)
    # utils.rm_rf(f"{root_outdir}/*")
    utils.rm_rf(vcf_root_out)
    utils.rm_rf(logs_root_out)
    os.mkdir(os.path.join(root_outdir, vcf_root_out))
    os.mkdir(os.path.join(root_outdir, logs_root_out))
    sample_number = 0
    tsv_out = os.path.join(root_outdir, "manifest.tsv")
    utils.rm_rf(tsv_out)
    json_out = os.path.join(root_outdir, "manifest.json")
    utils.rm_rf(json_out)
    data = {}
    parallel_jobs_data = []

    with open(sample_data_tsv) as f_in, open(tsv_out, "w") as f_out:
        print("sample", "vcf_file", "log_file", sep="\t", file=f_out)
        for line in f_in:
            if sample_number % samples_per_dir == 0:
                outdir = str(sample_number // samples_per_dir)
                vcf_dir = os.path.join(vcf_root_out, outdir)
                vcf_dir_full = os.path.join(root_outdir, vcf_dir)
                os.mkdir(vcf_dir_full)
                log_dir = os.path.join(logs_root_out, outdir)
                log_dir_full = os.path.join(root_outdir, log_dir)
                os.mkdir(log_dir_full)

            sample_name, minos_indir = line.rstrip().split()
            vcf_in = os.path.join(minos_indir, "debug.calls_with_zero_cov_alleles.vcf")
            log_in = os.path.join(minos_indir, "log.txt")
            vcf_out = os.path.join(vcf_dir, f"{sample_number}.vcf.gz")
            vcf_out_full = os.path.join(root_outdir, vcf_out)
            log_out = os.path.join(log_dir, f"{sample_number}.log.gz")
            log_out_full = os.path.join(root_outdir, log_out)
            parallel_jobs_data.append((vcf_in, vcf_out_full))
            parallel_jobs_data.append((log_in, log_out_full))
            sample_number += 1
            print(sample_name, vcf_out, log_out, sep="\t", file=f_out)
            data[sample_name] = {"vcf_file": vcf_out, "log_file": log_out}

    with multiprocessing.Pool(cpus) as pool:
        pool.map(compress_file, parallel_jobs_data)

    with open(json_out, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)

    if original_manifest is None:
        return

    expect_samples = manifest_to_set_of_sample_names(original_manifest)
    failed_samples = expect_samples.difference(data)
    if len(failed_samples) > 0:
        failed_samples = sorted(list(failed_samples))
        with open(os.path.join(root_outdir, "failed_samples.txt"), "w") as f:
            print(*failed_samples, sep="\n", file=f)
