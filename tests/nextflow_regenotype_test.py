import csv
import filecmp
import json
import os
import pytest
import subprocess

from cluster_vcf_records import vcf_file_read

from minos import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "nextflow_regenotype")
minos_dir = os.path.abspath(os.path.join(this_dir, os.pardir))


def _write_manifest(filename):
    with open(filename, "w") as f:
        print("name", "vcf", "reads", sep="\t", file=f)
        for i in range(5):
            vcf = os.path.join(data_dir, f"data.sample.{i}.vcf")
            if i < 4:
                bam = os.path.join(data_dir, f"data.sample.{i}.bam")
            else:
                bam = "does_not_exist"

            if i == 0:
                bam += "\t" + bam
            print(f"sample.{i}", vcf, bam, sep="\t", file=f)


def test_regenotype_pipeline():
    outdir = "tmp.nextflow_regeno_test.out"
    utils.rm_rf(outdir)
    os.mkdir(outdir)
    manifest = "tmp.nextflow_regeno_test.tsv"
    _write_manifest(os.path.join(outdir, manifest))

    regeno_nf = os.path.join(minos_dir, "nextflow", "regenotype.nf")
    regeno_config = os.path.join(minos_dir, "nextflow", "regenotype.config")
    dag = "tmp.nextflow_regeno_test.dag.pdf"
    ref_fasta = os.path.join(data_dir, "data.ref.fa")
    mask_bed = os.path.join(data_dir, "mask.bed")
    command = f"nextflow run -c {regeno_config} -profile tiny -with-dag {dag} {regeno_nf} --make_distance_matrix --mask_bed_file {mask_bed} --max_variants_per_sample 10 --ref_fasta {ref_fasta} --manifest {manifest} --outdir OUT"
    utils.syscall(command, cwd=outdir)

    expect_failed_samples = os.path.join(data_dir, "expect.failed_samples.txt")
    got_failed_samples = os.path.join(outdir, "OUT", "failed_samples.txt")
    assert filecmp.cmp(got_failed_samples, expect_failed_samples, shallow=False)

    expect_dist_matrix = os.path.join(data_dir, "expect.distance_matrix.txt")
    got_dist_matrix = os.path.join(outdir, "OUT", "distance_matrix.txt")
    assert filecmp.cmp(got_dist_matrix, expect_dist_matrix, shallow=False)

    # Don't know order of lines in the manifest tsv, or the filename that will
    # be given to each sample. We'll load in each VCF and check it matches the
    # sample name from the manifest. Also check info in json and tsv files
    # match
    manifest_json = os.path.join(outdir, "OUT", "manifest.json")
    assert os.path.exists(manifest_json)
    manifest_tsv = os.path.join(outdir, "OUT", "manifest.tsv")

    with open(manifest_json) as f:
        manifest_data = json.load(f)

    with open(manifest_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for d in reader:
            vcf = os.path.join(outdir, "OUT", d["vcf_file"])
            assert d["sample"] == vcf_file_read.get_sample_name_from_vcf_file(vcf)
            assert os.path.exists(os.path.join(outdir, "OUT", d["log_file"]))
            assert manifest_data[d["sample"]]["log_file"] == d["log_file"]
            assert manifest_data[d["sample"]]["vcf_file"] == d["vcf_file"]

    utils.rm_rf(outdir)
