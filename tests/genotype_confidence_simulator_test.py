import os
import unittest

from minos import genotype_confidence_simulator

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "genotype_confidence_simulator")


class TestGenotypeConfidenceSimulator(unittest.TestCase):
    def test_simulate_confidence_scores(self):
        """test _simulate_confidence_scores"""
        got = genotype_confidence_simulator.GenotypeConfidenceSimulator._simulate_confidence_scores(
            50, 300, 0.1, iterations=5
        )
        expected = [26, 31, 37, 46, 51]
        self.assertEqual(expected, got)

        #  Since the genotype confidence normalises by length, we shpuld get the same
        # results with allele lengths 1 and 2.
        got = genotype_confidence_simulator.GenotypeConfidenceSimulator._simulate_confidence_scores(
            50, 300, 0.1, iterations=5, allele_length=2
        )
        expected = [26, 31, 37, 46, 51]
        self.assertEqual(expected, got)

    def test_make_conf_to_percentile_dict(self):
        """test _make_conf_to_percentile_dict"""
        confidence_scores = [1, 1, 2, 3, 4, 4, 4, 5, 6, 8]
        got = genotype_confidence_simulator.GenotypeConfidenceSimulator._make_conf_to_percentile_dict(
            confidence_scores
        )
        expected = {1: 15, 2: 30, 3: 40, 4: 50, 5: 80, 6: 90, 8: 100}
        self.assertTrue(expected, got)

    def test_run_simulations_and_get_percentile_allele_length_1(self):
        """test run_simulations and get_percentile"""
        simulator = genotype_confidence_simulator.GenotypeConfidenceSimulator(
            50, 300, 0.1, iterations=5
        )
        simulator.run_simulations()
        expected_confidence_scores_percentiles = {
            26: 20.0,
            31: 40.0,
            37: 60.0,
            46: 80.0,
            51: 100.0,
        }
        self.assertEqual(
            expected_confidence_scores_percentiles,
            simulator.confidence_scores_percentiles,
        )
        self.assertEqual(20.00, simulator.get_percentile(26))
        self.assertEqual(40.00, simulator.get_percentile(31))
        # Try getting numbers that are not in the dict and will have to be inferred
        self.assertEqual(28.00, simulator.get_percentile(28))
        self.assertEqual(84.00, simulator.get_percentile(47))
        self.assertEqual(88.00, simulator.get_percentile(48))
        # Try values outside the range of what we already have
        self.assertEqual(0.00, simulator.get_percentile(25))
        self.assertEqual(0.00, simulator.get_percentile(24))
        self.assertEqual(100.00, simulator.get_percentile(51))
        self.assertEqual(100.00, simulator.get_percentile(52))

    def test_run_simulations_and_get_percentile_allele_length_1_no_het_calls(self):
        """test run_simulations and get_percentile with no het calls"""
        simulator = genotype_confidence_simulator.GenotypeConfidenceSimulator(
            50, 300, 0.1, iterations=5, call_hets=False,
        )
        simulator.run_simulations()
        expected_confidence_scores_percentiles = {
            149: 20.0, 164: 40.0, 193: 60.0, 200: 80.0, 215: 100.0,
        }
        self.assertEqual(
            expected_confidence_scores_percentiles,
            simulator.confidence_scores_percentiles,
        )
        self.assertEqual(80.00, simulator.get_percentile(200))

    def test_run_simulations_and_get_percentile_allele_length_2(self):
        """test run_simulations and get_percentile"""
        simulator = genotype_confidence_simulator.GenotypeConfidenceSimulator(
            50, 300, 0.1, allele_length=2, iterations=5
        )
        simulator.run_simulations()
        expected_confidence_scores_percentiles = {
            26: 20.0,
            31: 40.0,
            37: 60.0,
            46: 80.0,
            51: 100.0,
        }
        self.assertEqual(
            expected_confidence_scores_percentiles,
            simulator.confidence_scores_percentiles,
        )
        self.assertEqual(20.00, simulator.get_percentile(26))
        self.assertEqual(40.00, simulator.get_percentile(31))
        # Try getting numbers that are not in the dict and will have to be inferred
        self.assertEqual(28.00, simulator.get_percentile(28))
        self.assertEqual(84.00, simulator.get_percentile(47))
        self.assertEqual(88.00, simulator.get_percentile(48))
        # Try values outside the range of what we already have
        self.assertEqual(0.00, simulator.get_percentile(25))
        self.assertEqual(0.00, simulator.get_percentile(24))
        self.assertEqual(100.00, simulator.get_percentile(51))
        self.assertEqual(100.00, simulator.get_percentile(52))

    def test_simulations(self):
        """test simulations"""
        mean_depth = 50
        error_rate = 0.1
        depth_variance = 300
        allele_lengths = [1, 2, 10]
        simulations = genotype_confidence_simulator.Simulations(
            mean_depth, depth_variance, error_rate, allele_lengths=allele_lengths
        )
        self.assertEqual(allele_lengths, sorted(list(simulations.sims.keys())))
        self.assertEqual(
            simulations.sims[1].get_percentile(40), simulations.get_percentile(1, 40)
        )
        self.assertEqual(
            simulations.sims[2].get_percentile(40), simulations.get_percentile(2, 40)
        )
        # Check we get nearest allele length, when asking for an allele length that wasn't simulated
        self.assertEqual(
            simulations.sims[2].get_percentile(40), simulations.get_percentile(3, 40)
        )
