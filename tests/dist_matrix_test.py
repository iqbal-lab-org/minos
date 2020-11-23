import filecmp
import os
import pytest

import numpy as np

from minos import dist_matrix, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "dist_matrix")


def test_bed_mask_file_to_dict():
    infile = os.path.join(data_dir, "bed_mask_file_to_dict.bed")
    got = dist_matrix.bed_mask_file_to_dict(infile)
    expect = {
        "ref1": {0, 1, 10},
        "ref2": {42},
    }
    assert got == expect


def test_chrom_pos_in_mask():
    assert not dist_matrix.chrom_pos_in_mask("c1", 42, "A", {})
    assert not dist_matrix.chrom_pos_in_mask("c1", 42, "A", {"c2": {42}})
    assert dist_matrix.chrom_pos_in_mask("c1", 42, "A", {"c2": {42}, "c1": {42}})
    assert dist_matrix.chrom_pos_in_mask("c1", 41, "AT", {"c1": {42}})
    assert not dist_matrix.chrom_pos_in_mask("c1", 41, "AT", {"c1": {43}})


def test_gt_str_to_number():
    assert dist_matrix.gt_str_to_number("") == 0
    assert dist_matrix.gt_str_to_number(".") == 0
    assert dist_matrix.gt_str_to_number("0") == 1
    assert dist_matrix.gt_str_to_number("0/1") == 0
    assert dist_matrix.gt_str_to_number("0/0") == 1
    assert dist_matrix.gt_str_to_number("1/1") == 2


def test_load_genotypes_from_vcf_file():
    infile = os.path.join(data_dir, "load_genotypes_from_vcf_file.vcf")
    got = dist_matrix.load_genotypes_from_vcf_file(infile)
    expect = {
        "sample1": np.array([1, 2, 1, 2, 2], dtype=np.uint16),
        "sample2": np.array([0, 2, 1, 2, 2], dtype=np.uint16),
    }
    np.testing.assert_equal(got, expect)

    mask = {"ref": {0, 3, 4}}
    got = dist_matrix.load_genotypes_from_vcf_file(infile, mask)
    expect = {
        "sample1": np.array([2, 1], dtype=np.uint16),
        "sample2": np.array([2, 1], dtype=np.uint16),
    }
    np.testing.assert_equal(got, expect)


def test_distance_matrix_from_vcf_file():
    vcf_file = os.path.join(data_dir, "distance_matrix_from_vcf_file.vcf")
    outfile = "tmp.distance_matrix_from_vcf_file.out"
    utils.rm_rf(outfile)
    dist_matrix.distance_matrix_from_vcf_file(vcf_file, outfile)
    expect = os.path.join(data_dir, "distance_matrix_from_vcf_file.expect.no_mask")
    assert filecmp.cmp(outfile, expect, shallow=False)
    os.unlink(outfile)

    mask_bed = os.path.join(data_dir, "distance_matrix_from_vcf_file.mask.bed")
    dist_matrix.distance_matrix_from_vcf_file(vcf_file, outfile, mask_bed_file=mask_bed)
    expect = os.path.join(data_dir, "distance_matrix_from_vcf_file.expect.mask")
    assert filecmp.cmp(outfile, expect, shallow=False)
    os.unlink(outfile)
