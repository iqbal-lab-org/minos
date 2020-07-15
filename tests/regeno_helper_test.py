import filecmp
import json
import os
import pytest

import pyfastaq

from minos import regeno_helper, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "regeno_helper")


def test_compress_file():
    """test compress_file"""
    vcf_in = os.path.join(data_dir, "compress_file.vcf")
    vcf_out = "tmp.compress_file.vcf.gz"
    txt_in = os.path.join(data_dir, "compress_file.txt")
    txt_out = "tmp.compress_file.txt.gz"
    utils.rm_rf(vcf_out, txt_out)
    regeno_helper.compress_file((vcf_in, vcf_out))
    regeno_helper.compress_file((txt_in, txt_out))
    assert os.path.exists(vcf_out)
    assert os.path.exists(txt_out)
    os.unlink(vcf_out)
    os.unlink(txt_out)


def _file_to_lines(infile):
    f = pyfastaq.utils.open_file_read(infile)
    lines = [x.rstrip() for x in f]
    f.close()
    return lines


def _file_contents_the_same(file1, file2):
    return _file_to_lines(file1) == _file_to_lines(file2)


def test_make_per_sample_vcfs_dir():
    """test make_per_sample_vcfs_dir"""
    manifest_file = "tmp.make_per_sample_vcfs_dir.tsv"
    indir = os.path.join(data_dir, "make_per_sample_vcfs_dir")
    minos_indirs = {}
    with open(manifest_file, "w") as f:
        for i in range(1, 6):
            minos_dir = os.path.join(indir, f"minos.{i}")
            print(f"sample.{i}", minos_dir, sep="\t", file=f)
            minos_indirs[f"sample.{i}"] = minos_dir

    root_out = "tmp.make_per_sample_vcfs_dir.out"
    utils.rm_rf(root_out)
    regeno_helper.make_per_sample_vcfs_dir(
        manifest_file, root_out, samples_per_dir=2, cpus=2
    )
    os.unlink(manifest_file)
    expect_tsv = os.path.join(data_dir, "make_per_sample_vcfs_dir.expect.tsv")
    got_tsv = os.path.join(root_out, "manifest.tsv")
    assert filecmp.cmp(expect_tsv, got_tsv, shallow=False)
    expect_json = os.path.join(data_dir, "make_per_sample_vcfs_dir.expect.json")
    got_json = os.path.join(root_out, "manifest.json")
    assert filecmp.cmp(expect_json, got_json, shallow=False)
    with open(got_json) as f:
        json_data = json.load(f)

    for sample, minos_dir in minos_indirs.items():
        original_vcf = os.path.join(minos_dir, "debug.calls_with_zero_cov_alleles.vcf")
        original_log = os.path.join(minos_dir, "log.txt")
        new_vcf = os.path.join(root_out, json_data[sample]["vcf_file"])
        new_log = os.path.join(root_out, json_data[sample]["log_file"])
        assert _file_contents_the_same(original_vcf, new_vcf)
        assert _file_contents_the_same(original_log, new_log)

    utils.rm_rf(root_out)
