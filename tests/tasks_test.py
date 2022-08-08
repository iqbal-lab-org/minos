import os
import pytest
from unittest import mock

from minos import tasks, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "tasks")


def test_vcf_merge_and_cluster():
    # Merge and cluster are tested properly in `cluster_vcf_records` repo.
    # Here, we'll check they run and make the expected output files (but
    # don't check the contents of the files).
    options = mock.Mock()
    options.ref_fasta = os.path.join(data_dir, "merge_and_cluster.ref.fa")
    options.outdir = "tmp.merge_and_cluster"
    utils.rm_rf(options.outdir)
    options.vcf_fofn = "tmp.merge_and_cluster.vcf.fofn"
    with open(options.vcf_fofn, "w") as f:
        for i in (1, 2, 3):
            print(os.path.join(data_dir, f"merge_and_cluster.{i}.vcf"), file=f)

    options.temp_dir = None
    options.cpus = 2
    options.mem_limit = 1
    options.force = False
    options.sample_limit = 2
    tasks.vcf_merge.run(options)
    expect = [
        "block.0.tsv.gz",
        "block.0.tsv.gz.tbi",
        "block.1.tsv.gz",
        "block.1.tsv.gz.tbi",
        "metadata.json",
        "variants.tsv.gz",
    ]
    for fname in expect:
        assert os.path.exists(os.path.join(options.outdir, fname))
    os.unlink(options.vcf_fofn)

    options = mock.Mock()
    options.ref_fasta = os.path.join(data_dir, "merge_and_cluster.ref.fa")
    options.merge_dir = "tmp.merge_and_cluster"
    options.outprefix = "tmp.merge_and_cluster.out"
    options.max_ref_len = 6
    options.max_alleles = 50
    options.cpus = 1
    expect = [f"{options.outprefix}.excluded.tsv", f"{options.outprefix}.vcf"]
    utils.rm_rf(*expect)
    tasks.vcf_cluster.run(options)
    for fname in expect:
        assert os.path.exists(fname)
        os.unlink(fname)
    utils.rm_rf(options.merge_dir)


def test_get_test_data():
    options = mock.Mock()
    options.outdir = "tmp.get_test_data"
    utils.syscall(f"rm -rf {options.outdir}")
    tasks.get_test_data.run(options)
    outdir_abs = os.path.abspath(options.outdir)
    assert os.path.exists(options.outdir)
    with open(os.path.join(options.outdir, "manifest.tsv")) as f:
        got_lines = [x.rstrip().split("\t") for x in f]
    assert got_lines == [
        ["name", "reads", "vcf"],
        ["sample1", os.path.join(outdir_abs, "sample1.bam"), os.path.join(outdir_abs, "in.1.vcf")],
        ["sample2", os.path.join(outdir_abs, "sample2.bam"), os.path.join(outdir_abs, "in.2.vcf")],
    ]
    utils.syscall(f"rm -r {options.outdir}")
