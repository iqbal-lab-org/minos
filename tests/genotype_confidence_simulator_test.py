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
    expected = [99, 114, 143, 200, 215]
    assert got == expected

    #  Since the genotype confidence normalises by length, we shpuld get the same
    # results with allele lengths 1 and 2.
    got = genotype_confidence_simulator.GenotypeConfidenceSimulator._simulate_confidence_scores(
        50, 300, 0.1, iterations=5, allele_length=2
    )
    expected = [99, 114, 143, 200, 215]
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
        50, 300, 0.01, iterations=5, call_hets=True
    )
    simulator.run_simulations()
    expected_confidence_scores_percentiles = {
        55: 20.0,
        256: 40.0,
        285: 60.0,
        337: 80.0,
        368: 100.0,
    }
    assert (
        simulator.confidence_scores_percentiles
        == expected_confidence_scores_percentiles
    )
    assert simulator.get_percentile(55) == 20.00
    assert simulator.get_percentile(256) == 40.00
    # Try getting number that is not in the dict and will have to be inferred
    assert simulator.get_percentile(150) == 29.45
    # Try values outside the range of what we already have
    simulator.get_percentile(53) == 0.00
    simulator.get_percentile(52) == 0.00
    simulator.get_percentile(369) == 100.00
    simulator.get_percentile(370) == 100.00


def test_run_simulations_and_get_percentile_allele_length_1_no_het_calls():
    """test run_simulations and get_percentile with no het calls"""
    simulator = genotype_confidence_simulator.GenotypeConfidenceSimulator(
        50, 300, 0.1, iterations=5, call_hets=False,
    )
    simulator.run_simulations()
    expected_confidence_scores_percentiles = {
        99: 20.0,
        114: 40.0,
        143: 60.0,
        200: 80.0,
        215: 100.0,
    }
    assert (
        simulator.confidence_scores_percentiles
        == expected_confidence_scores_percentiles
    )
    assert simulator.get_percentile(200) == 80.00


def test_run_simulations_and_get_percentile_allele_length_2():
    """test run_simulations and get_percentile"""
    simulator = genotype_confidence_simulator.GenotypeConfidenceSimulator(
        50, 300, 0.01, allele_length=2, iterations=5, call_hets=True,
    )
    simulator.run_simulations()
    expected_confidence_scores_percentiles = {
        55: 20.0,
        256: 40.0,
        285: 60.0,
        337: 80.0,
        368: 100.0,
    }
    assert (
        simulator.confidence_scores_percentiles
        == expected_confidence_scores_percentiles
    )
    assert simulator.get_percentile(55) == 20.00
    assert simulator.get_percentile(256) == 40.00
    # Try getting number that is not in the dict and will have to be inferred
    assert simulator.get_percentile(150) == 29.45
    # Try values outside the range of what we already have
    simulator.get_percentile(53) == 0.00
    simulator.get_percentile(52) == 0.00
    simulator.get_percentile(369) == 100.00
    simulator.get_percentile(370) == 100.00


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
