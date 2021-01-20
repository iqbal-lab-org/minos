import os
import pytest

from minos import genotype_confidence_simulator

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "genotype_confidence_simulator")


def test_simulate_confidence_scores():
    """test _simulate_confidence_scores"""
    got = genotype_confidence_simulator.GenotypeConfidenceSimulator._simulate_confidence_scores(
        50, 300, 0.1, iterations=5
    )
    expected = [76, 90, 119, 140, 156]
    assert got == expected

    #  Since the genotype confidence normalises by length, we shpuld get the same
    # results with allele lengths 1 and 2.
    got = genotype_confidence_simulator.GenotypeConfidenceSimulator._simulate_confidence_scores(
        50, 300, 0.1, iterations=5, allele_length=2
    )
    assert got == expected


def test_make_conf_to_percentile_dict():
    """test _make_conf_to_percentile_dict"""
    confidence_scores = [1, 1, 2, 3, 4, 4, 4, 5, 6, 8]
    got = genotype_confidence_simulator.GenotypeConfidenceSimulator._make_conf_to_percentile_dict(
        confidence_scores
    )
    expected = {1: 15, 2: 30, 3: 40, 4: 60, 5: 80, 6: 90, 8: 100}
    assert got == expected


def test_run_simulations_and_get_percentile_allele_length_1():
    """test run_simulations and get_percentile"""
    simulator = genotype_confidence_simulator.GenotypeConfidenceSimulator(
        50, 300, 0.01, iterations=5
    )
    simulator.run_simulations()
    expected_confidence_scores_percentiles = {
        193: 20.0,
        221: 40.0,
        271: 60.0,
        278: 80.0,
        303: 100.0
    }
    assert (
        simulator.confidence_scores_percentiles
        == expected_confidence_scores_percentiles
    )
    assert simulator.get_percentile(193) == 20.00
    assert simulator.get_percentile(221) == 40.00
    # Try getting number that is not in the dict and will have to be inferred
    assert simulator.get_percentile(207) == 30.0
    # Try values outside the range of what we already have
    simulator.get_percentile(192) == 0.00
    simulator.get_percentile(191) == 0.00
    simulator.get_percentile(304) == 100.00
    simulator.get_percentile(305) == 100.00


def test_run_simulations_and_get_percentile_allele_length_2():
    """test run_simulations and get_percentile"""
    simulator = genotype_confidence_simulator.GenotypeConfidenceSimulator(
        50, 300, 0.01, allele_length=2, iterations=5
    )
    simulator.run_simulations()
    expected_confidence_scores_percentiles = {
        193: 20.0,
        221: 40.0,
        271: 60.0,
        278: 80.0,
        303: 100.0
    }
    assert (
        simulator.confidence_scores_percentiles
        == expected_confidence_scores_percentiles
    )
    assert simulator.get_percentile(193) == 20.00
    assert simulator.get_percentile(221) == 40.00
    # Try getting number that is not in the dict and will have to be inferred
    assert simulator.get_percentile(207) == 30.0
    # Try values outside the range of what we already have
    simulator.get_percentile(192) == 0.00
    simulator.get_percentile(191) == 0.00
    simulator.get_percentile(304) == 100.00
    simulator.get_percentile(305) == 100.00


def test_simulations():
    """test simulations"""
    mean_depth = 50
    error_rate = 0.1
    depth_variance = 300
    allele_lengths = [1, 2, 10]
    simulations = genotype_confidence_simulator.Simulations(
        mean_depth, depth_variance, error_rate, allele_lengths=allele_lengths
    )
    assert sorted(list(simulations.sims.keys())) == allele_lengths
    assert simulations.get_percentile(1, 40) == simulations.sims[1].get_percentile(40)
    assert simulations.get_percentile(2, 40) == simulations.sims[2].get_percentile(40)
    # Check we get nearest allele length, when asking for an allele length that wasn't simulated
    assert simulations.get_percentile(3, 40) == simulations.sims[2].get_percentile(40)
