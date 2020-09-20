import filecmp
import json
import os
import pytest

import pyfastaq

from minos import regeno_helper, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "regeno_helper")


def test_vcf_has_too_many_variants():
    vcf_file = os.path.join(data_dir, "vcf_has_too_many_variants.vcf")
    assert not regeno_helper.vcf_has_too_many_variants(vcf_file, 5, 0.5, 10)
    assert regeno_helper.vcf_has_too_many_variants(vcf_file, 5, 0.49, 10)
    assert regeno_helper.vcf_has_too_many_variants(vcf_file, 4, 0.5, 10)
    assert regeno_helper.vcf_has_too_many_variants(vcf_file, 5, 0.5, 9)


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


def test_parse_manifest_file():
    """test parse_manifest_file"""
    manifest_tsv = "tmp.parse_manifest_file.in.tsv"
    with open(manifest_tsv, "w") as f:
        vcf_prefix =  os.path.join(data_dir, "parse_manifest_file")
        vcf1 = f"{vcf_prefix}.1.vcf"
        vcf2 = f"{vcf_prefix}.2.vcf"
        vcf3 = f"{vcf_prefix}.3.vcf"
        print("name", "vcf", "reads", sep="\t", file=f)
        print("sample1", vcf1, "1.reads.fq", sep="\t", file=f)
        print("sample2", vcf2, "2.reads.1.fq 2.reads.2.fq", sep="\t", file=f)
        print("sample3", vcf3, "3.reads.1.fq 3.reads.2.fq", sep="\t", file=f)
    merge_out = "tmp.parse_manifest_file.merge.fofn"
    adjudicate_out = "tmp.parse_manifest_file.adjudicate.tsv"
    ref_fasta = os.path.join(data_dir, "parse_manifest_file.ref.fa")
    utils.rm_rf(merge_out, adjudicate_out)
    regeno_helper.parse_manifest_file(manifest_tsv, merge_out, adjudicate_out, ref_fasta)
    os.unlink(manifest_tsv)
    expect_adj = os.path.join(data_dir, "parse_manifest_file.out.tsv")
    assert filecmp.cmp(expect_adj, adjudicate_out, shallow=False)
    os.unlink(adjudicate_out)
    with open(merge_out) as f:
        got_lines = [x.rstrip() for x in f]
    assert got_lines == [vcf1, vcf3]
    os.unlink(merge_out)


def test_manifest_to_set_of_sample_names():
    """test manifest_to_set_of_sample_names"""
    infile = os.path.join(data_dir, "manifest_to_set_of_sample_names.tsv")
    expect = {"sample1", "sample2", "sample3"}
    assert regeno_helper.manifest_to_set_of_sample_names(infile) == expect


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
