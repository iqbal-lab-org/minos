import os
import pytest

from minos import genotyper

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "genotyper")


def test_init():
    """test init"""
    gtyper = genotyper.Genotyper(0, 20, 0.0001)
    assert gtyper.min_cov_more_than_error == 0
    assert gtyper.no_of_successes == 0
    assert gtyper.prob_of_success == 0

    gtyper = genotyper.Genotyper(10, 20, 0.0001)
    assert gtyper.no_of_successes == 10
    assert gtyper.prob_of_success == 0.5
    assert gtyper.min_cov_more_than_error == 1

    gtyper = genotyper.Genotyper(10, 20, 0.01)
    assert gtyper.no_of_successes == 10
    assert gtyper.prob_of_success == 0.5
    assert gtyper.min_cov_more_than_error == 2

    gtyper = genotyper.Genotyper(100, 200, 0.001)
    assert gtyper.no_of_successes == 100
    assert gtyper.prob_of_success == 0.5
    assert gtyper.min_cov_more_than_error == 8

    # variance < mean, so will hit the code where it forces
    # variance = 2 * mean = 20
    gtyper = genotyper.Genotyper(10, 5, 0.01)
    assert gtyper.no_of_successes == 10
    assert gtyper.prob_of_success == 0.5
    assert gtyper.min_cov_more_than_error == 2


def test_singleton_alleles_and_coverage():
    """test _singleton_alleles_and_coverage"""
    allele_combination_cov = {"1": 20, "3": 1}
    allele_groups_dict = {"1": {0}, "2": {1}, "3": {1, 2}, "4": {5, 6}}
    got = genotyper.Genotyper._singleton_alleles_and_coverage(
        allele_combination_cov, allele_groups_dict
    )
    assert got == {0: 20}

    allele_combination_cov["2"] = 42
    got = genotyper.Genotyper._singleton_alleles_and_coverage(
        allele_combination_cov, allele_groups_dict
    )
    assert got == {0: 20, 1: 42}


def test_total_coverage():
    """test _total_coverage"""
    f = genotyper.Genotyper._total_coverage
    assert f({}) == 0
    assert f({"x": 1}) == 1
    assert f({"x": 1, "y": 41}) == 42


def test_haploid_allele_coverages():
    """test  _haploid_allele_coverages"""
    allele_combination_cov = {"1": 20, "2": 1}
    allele_groups_dict = {"0": {0}, "1": {1}, "2": {1, 2}, "3": {5, 6}}
    num_distinct_alleles = 7  # 1 + the max allele index
    got = genotyper.Genotyper._haploid_allele_coverages(
        num_distinct_alleles, allele_combination_cov, allele_groups_dict
    )
    assert got == [0, 21, 1, 0, 0, 0, 0]


def test_log_likelihood_homozygous():
    """test _log_likelihood_homozygous"""
    gtyper = genotyper.Genotyper(100, 200, 0.01)
    allele_depth = 90
    total_depth = 95
    allele_length = 5
    non_zeros = allele_length
    got = gtyper._log_likelihood_homozygous(
        allele_depth, total_depth, allele_length, non_zeros
    )
    assert round(got, 2) == -26.78

    gtyper = genotyper.Genotyper(10, 200, 0.01)
    allele_depth = 1
    total_depth = 9
    got = gtyper._log_likelihood_homozygous(
        allele_depth, total_depth, allele_length, non_zeros
    )
    assert round(got, 2) == -39.34


def test_calculate_log_likelihoods():
    """test _calculate_log_likelihoods"""
    gtyper = genotyper.Genotyper(20, 40, 0.01)
    allele_combination_cov = {"1": 2, "2": 20, "3": 1}
    allele_groups_dict = {"1": {0}, "2": {1}, "3": {0, 1}, "4": {2}}
    allele_per_base_cov = [[0, 1], [20, 19]]
    depth0 = round(3 / 23, 4)
    depth01 = 1
    depth1 = round(21 / 23, 4)
    gtyper._init_alleles_and_genotypes(
        allele_combination_cov=allele_combination_cov,
        allele_per_base_cov=allele_per_base_cov,
        allele_groups_dict=allele_groups_dict,
    )
    gtyper._calculate_log_likelihoods()
    assert len(gtyper.likelihoods) == 2
    expected = [
        ({1}, -12.03, depth1),
        ({0}, -114.57, depth0),
    ]
    gtyper.likelihoods = [(x[0], round(x[1], 2), x[2]) for x in gtyper.likelihoods]
    assert gtyper.likelihoods == expected


def test_run():
    """test run"""
    gtyper = genotyper.Genotyper(20, 40, 0.01)
    allele_combination_cov = {"1": 2, "2": 20, "3": 1}
    allele_groups_dict = {"1": {0}, "2": {1}, "3": {0, 1}, "4": {2}}
    allele_per_base_cov = [[0, 1], [20, 19]]
    gtyper.run(allele_combination_cov, allele_per_base_cov, allele_groups_dict)
    depth0 = round(3 / 23, 4)
    depth1 = round(21 / 23, 4)
    expected = [({1}, -12.03, depth1), ({0}, -114.57, depth0)]
    assert len(gtyper.likelihoods) == len(expected)
    for i in range(len(expected)):
        assert gtyper.likelihoods[i][0] == expected[i][0]
        assert round(gtyper.likelihoods[i][1], 2) == round(expected[i][1], 2)
        assert gtyper.likelihoods[i][2] == expected[i][2]


def test_run_zero_coverage():
    """test run when all alleles have zero coverage"""
    gtyper = genotyper.Genotyper(20, 40, 0.01)
    allele_combination_cov = {}
    allele_groups_dict = {"1": {0}, "2": {1}, "3": {0, 1}, "4": {2}}
    allele_per_base_cov = [[0], [0, 0]]
    gtyper.run(allele_combination_cov, allele_per_base_cov, allele_groups_dict)
    assert gtyper.genotype == {"."}
    assert gtyper.genotype_confidence == 0.0
    assert gtyper.genotype_frs == "."


def test_nomatherror_mean_depth0():
    """
    Can get a mean_depth of zero but try to genotype a non-zero coverage site due to rounding imprecision.
    In which case we need to avoid trying to do log(0) in likelihood calculation and should return no call.
    """
    gtyper = genotyper.Genotyper(0, 0, 0.01)
    allele_combination_cov = {"1": 1}
    allele_groups_dict = {"1": {0}, "2": {1}}
    allele_per_base_cov = [[1], [0, 0]]
    gtyper.run(allele_combination_cov, allele_per_base_cov, allele_groups_dict)
    assert gtyper.genotype == {"."}
    assert gtyper.genotype_confidence == 0.0
    assert gtyper.genotype_frs == "."
