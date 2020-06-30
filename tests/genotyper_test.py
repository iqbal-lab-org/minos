import os
import pytest

from minos import genotyper

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "genotyper")


def test_init():
    """test init"""
    gtyper = genotyper.Genotyper(0, 20, 0.0001)
    assert gtyper.use_nbinom
    assert gtyper.min_cov_more_than_error == 0
    assert gtyper.no_of_successes == 0
    assert gtyper.prob_of_success == 0

    gtyper = genotyper.Genotyper(10, 20, 0.0001)
    assert gtyper.use_nbinom
    assert gtyper.no_of_successes == 10
    assert gtyper.prob_of_success == 0.5
    assert gtyper.min_cov_more_than_error == 1

    gtyper = genotyper.Genotyper(10, 20, 0.01)
    assert gtyper.use_nbinom
    assert gtyper.no_of_successes == 10
    assert gtyper.prob_of_success == 0.5
    assert gtyper.min_cov_more_than_error == 2

    gtyper = genotyper.Genotyper(100, 200, 0.001)
    assert gtyper.use_nbinom
    assert gtyper.no_of_successes == 100
    assert gtyper.prob_of_success == 0.5
    assert gtyper.min_cov_more_than_error == 8

    gtyper = genotyper.Genotyper(10, 5, 0.01)
    assert not gtyper.use_nbinom
    assert gtyper.no_of_successes == None
    assert gtyper.prob_of_success == None
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


def test_coverage_of_diploid_alleles_equal_dispatching():
    """
    test _coverage_of_diploid_alleles function
    # Equiv class -> coverage:
    #   {0,1}: 9
    We collect coverages of alleles 0 and 1.
    This is an edge case where the only coverage is on both: then we dispatch equally.
    """
    allele_combination_cov = {"1": 9}
    allele_groups_dict = {"1": {0, 1}}
    got = genotyper.Genotyper._coverage_of_diploid_alleles(
        0, 1, allele_combination_cov, allele_groups_dict,
    )
    assert got == (4.5, 4.5)


def test_coverage_of_diploid_alleles_correct_dispatching():
    """
    test _coverage_of_diploid_alleles function
    # Equiv class -> coverage:
    #   {0}: 17
    #   {1}: 80
    #   {0,1}: 10
    #   {0,2): 3
    We collect coverage for alleles 0 and 1: 20 units specific to 0, 80 specific to 1.
    Thus we expect a 1:4 dispatching ratio of the 10 reads which agree with both 0 and 1.
    """
    allele_combination_cov = {"1": 17, "2": 80, "3": 10, "4": 3}
    allele_groups_dict = {"1": {0}, "2": {1}, "3": {0, 1}, "4": {0, 2}}
    got = genotyper.Genotyper._coverage_of_diploid_alleles(
        0, 1, allele_combination_cov, allele_groups_dict,
    )
    assert got == (22, 88)


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
    assert round(got, 2) == -26.71

    gtyper = genotyper.Genotyper(10, 200, 0.01)
    allele_depth = 1
    total_depth = 9
    got = gtyper._log_likelihood_homozygous(
        allele_depth, total_depth, allele_length, non_zeros
    )
    assert round(got, 2) == -44.77

    gtyper = genotyper.Genotyper(10, 200, 0.01, force_poisson=True)
    allele_depth = 1
    total_depth = 9
    got = gtyper._log_likelihood_homozygous(
        allele_depth, total_depth, allele_length, non_zeros
    )
    assert round(got, 2) == -44.54


def test_log_likelihood_heterozygous():
    """test _log_likelihood_heterozygous"""
    gtyper = genotyper.Genotyper(100, 200, 0.01)
    allele_depth1 = 45
    allele_depth2 = 40
    total_depth = 95
    allele_length1 = 3
    allele_length2 = 3
    non_zeros1 = 3
    non_zeros2 = 3
    got = gtyper._log_likelihood_heterozygous(
        allele_depth1,
        allele_depth2,
        total_depth,
        allele_length1,
        allele_length2,
        non_zeros1,
        non_zeros2,
    )
    assert round(got, 2) == -52.97

    non_zeros1 = 2
    non_zeros2 = 2
    got = gtyper._log_likelihood_heterozygous(
        allele_depth1,
        allele_depth2,
        total_depth,
        allele_length1,
        allele_length2,
        non_zeros1,
        non_zeros2,
    )
    assert round(got, 2) == -86.31


def test_calculate_log_likelihoods():
    """test _calculate_log_likelihoods"""
    gtyper = genotyper.Genotyper(20, 40, 0.01, call_hets=True)
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
    assert len(gtyper.likelihoods) == 3
    expected = [
        ({1}, -11.68, depth1),
        ({0, 1}, -22.93, depth01),
        ({0}, -124.91, depth0),
    ]
    gtyper.likelihoods = [(x[0], round(x[1], 2), x[2]) for x in gtyper.likelihoods]
    assert gtyper.likelihoods == expected


def test_run_with_call_hets_false():
    """test run call_hets False"""
    gtyper = genotyper.Genotyper(20, 40, 0.01, call_hets=False)
    allele_combination_cov = {"1": 2, "2": 20, "3": 1}
    allele_groups_dict = {"1": {0}, "2": {1}, "3": {0, 1}, "4": {2}}
    allele_per_base_cov = [[0, 1], [20, 19]]
    gtyper.run(allele_combination_cov, allele_per_base_cov, allele_groups_dict)
    depth0 = round(3 / 23, 4)
    depth1 = round(21 / 23, 4)
    expected = [({1}, -11.68, depth1), ({0}, -124.91, depth0)]
    assert len(gtyper.likelihoods) == len(expected)
    for i in range(len(expected)):
        assert gtyper.likelihoods[i][0] == expected[i][0]
        assert round(gtyper.likelihoods[i][1], 2) == round(expected[i][1], 2)
        assert gtyper.likelihoods[i][2] == expected[i][2]


def test_run():
    """test run call_hets True"""
    gtyper = genotyper.Genotyper(20, 40, 0.01, call_hets=True)
    allele_combination_cov = {"1": 2, "2": 20, "3": 1}
    allele_groups_dict = {"1": {0}, "2": {1}, "3": {0, 1}, "4": {2}}
    allele_per_base_cov = [[0, 1], [20, 19]]
    gtyper.run(allele_combination_cov, allele_per_base_cov, allele_groups_dict)
    depth0 = round(3 / 23, 4)
    depth01 = 1
    depth1 = round(21 / 23, 4)
    expected = [
        ({1}, -11.68, depth1),
        ({0, 1}, -22.93, depth01),
        ({0}, -124.91, depth0),
    ]
    assert len(gtyper.likelihoods) == len(expected)
    for i in range(len(expected)):
        assert gtyper.likelihoods[i][0] == expected[i][0]
        assert round(gtyper.likelihoods[i][1], 2) == round(expected[i][1], 2)
        assert gtyper.likelihoods[i][2] == expected[i][2]


def test_run_zero_coverage():
    """test run when all alleles have zero coverage"""
    gtyper = genotyper.Genotyper(20, 40, 0.01, call_hets=True)
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
    gtyper = genotyper.Genotyper(0, 0, 0.01, call_hets=True)
    allele_combination_cov = {"1": 1}
    allele_groups_dict = {"1": {0}, "2": {1}}
    allele_per_base_cov = [[1], [0, 0]]
    gtyper.run(allele_combination_cov, allele_per_base_cov, allele_groups_dict)
    assert gtyper.genotype == {"."}
    assert gtyper.genotype_confidence == 0.0
    assert gtyper.genotype_frs == "."
