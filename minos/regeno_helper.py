import json
import multiprocessing
import os

from minos import utils


def compress_file(filenames):
    infile, outfile = filenames
    zipper = "bgzip" if infile.endswith(".vcf") else "gzip -9"
    utils.syscall(f"{zipper} -c {infile} > {outfile}")


def make_per_sample_vcfs_dir(
    sample_data_tsv, root_outdir, samples_per_dir=1000, cpus=1
):
    utils.rm_rf(root_outdir)
    vcf_root_out = os.path.join("VCFs")
    logs_root_out = os.path.join("Logs")
    os.mkdir(root_outdir)
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
            print(sample_number, sample_name, vcf_out, log_out, sep="\t", file=f_out)
            data[sample_name] = {"vcf_file": vcf_out, "log_file": log_out}

    with multiprocessing.Pool(cpus) as pool:
        pool.map(compress_file, parallel_jobs_data)

    with open(os.path.join(root_outdir, "manifest.json"), "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)
