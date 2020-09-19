import csv
import json
import multiprocessing
import os

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
                gt = set(record.FORMAT["GT"].split("/"))
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


def parse_manifest_file(infile, merge_fofn, adjudicate_tsv):
    with open(infile) as f_in, open(merge_fofn, "w") as f_out_merge, open(
        adjudicate_tsv, "w"
    ) as f_out_adj:
        print("name", "reads", sep="\t", file=f_out_adj)
        reader = csv.DictReader(f_in, delimiter="\t")
        for d in reader:
            if "vcf" in d:
                print(d["vcf"], file=f_out_merge)

            reads = "--reads " + " --reads ".join(d["reads"].split())
            print(d["name"], reads, sep="\t", file=f_out_adj)


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
    utils.rm_rf(f"{root_outdir}/*")
    os.mkdir(os.path.join(root_outdir, vcf_root_out))
    os.mkdir(os.path.join(root_outdir, logs_root_out))
    sample_number = 0
    tsv_out = os.path.join(root_outdir, "manifest.tsv")
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

    with open(os.path.join(root_outdir, "manifest.json"), "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)

    if original_manifest is None:
        return

    expect_samples = manifest_to_set_of_sample_names(original_manifest)
    failed_samples = expect_samples.difference(data)
    if len(failed_samples) > 0:
        failed_samples = sorted(list(failed_samples))
        with open(os.path.join(root_outdir, "failed_samples.txt"), "w") as f:
            print(*failed_samples, sep="\n", file=f)
