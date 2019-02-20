import shutil
import os
import unittest

from minos import genotype_confidence_simulator

modules_dir = os.path.dirname(os.path.abspath(genotype_confidence_simulator.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'genotype_confidence_simulator')

class TestGenotypeConfidenceSimulator(unittest.TestCase):
    def test_simulate_confidence_scores(self):
        '''test _simulate_confidence_scores'''
        got = genotype_confidence_simulator.GenotypeConfidenceSimulator._simulate_confidence_scores(50, 0.1, 5)
        expected = [31, 34, 40, 46, 56]
        self.assertEqual(expected, got)

        got = genotype_confidence_simulator.GenotypeConfidenceSimulator._simulate_confidence_scores(50, 0.1, 5, allele_length=2)
        expected = [56, 59, 65, 71, 81]
        self.assertEqual(expected, got)


    def test_make_conf_to_percentile_dict(self):
        '''test _make_conf_to_percentile_dict'''
        confidence_scores = [1, 1, 2, 3, 4, 4, 4, 5, 6, 8]
        got = genotype_confidence_simulator.GenotypeConfidenceSimulator._make_conf_to_percentile_dict(confidence_scores)
        expected = {1: 15, 2: 30, 3: 40, 4: 50, 5: 80, 6: 90, 8: 100}
        self.assertTrue(expected, got)


    def test_run_simulations_and_get_percentile_allele_length_1(self):
        '''test run_simulations and get_percentile'''
        simulator = genotype_confidence_simulator.GenotypeConfidenceSimulator(50, 0.1, iterations=5)
        simulator.run_simulations()
        expected_confidence_scores_percentiles = {31: 20.0, 34: 40.0, 40: 60.0, 46: 80.0, 56: 100.0}
        self.assertEqual(expected_confidence_scores_percentiles, simulator.confidence_scores_percentiles)
        self.assertEqual(20.00, simulator.get_percentile(31))
        self.assertEqual(40.00, simulator.get_percentile(34))
        # Try getting numbers that are not in the dict and will have to be inferred
        self.assertEqual(26.67, simulator.get_percentile(32))
        self.assertEqual(90.00, simulator.get_percentile(51))
        self.assertEqual(92.00, simulator.get_percentile(52))
        # Try values outside the range of what we already have
        self.assertEqual(0.00, simulator.get_percentile(30))
        self.assertEqual(0.00, simulator.get_percentile(29))
        self.assertEqual(100.00, simulator.get_percentile(57))
        self.assertEqual(100.00, simulator.get_percentile(58))


    def test_run_simulations_and_get_percentile_allele_length_2(self):
        '''test run_simulations and get_percentile'''
        simulator = genotype_confidence_simulator.GenotypeConfidenceSimulator(50, 0.1, allele_length=2, iterations=5)
        simulator.run_simulations()
        expected_confidence_scores_percentiles = {56: 20.0, 59: 40.0, 65: 60.0, 71: 80.0, 81: 100.0}
        self.assertEqual(expected_confidence_scores_percentiles, simulator.confidence_scores_percentiles)
        self.assertEqual(20.00, simulator.get_percentile(56))
        self.assertEqual(40.00, simulator.get_percentile(59))
        # Try getting numbers that are not in the dict and will have to be inferred
        self.assertEqual(26.67, simulator.get_percentile(57))
        self.assertEqual(90.00, simulator.get_percentile(76))
        self.assertEqual(92.00, simulator.get_percentile(77))
        # Try values outside the range of what we already have
        self.assertEqual(0.00, simulator.get_percentile(55))
        self.assertEqual(0.00, simulator.get_percentile(54))
        self.assertEqual(100.00, simulator.get_percentile(82))
        self.assertEqual(100.00, simulator.get_percentile(83))


    def test_simulations(self):
        '''test simulations'''
        mean_depth = 50
        error_rate = 0.1
        allele_lengths = [1, 2, 10]
        simulations = genotype_confidence_simulator.Simulations(mean_depth, error_rate, allele_lengths=allele_lengths)
        self.assertEqual(allele_lengths, sorted(list(simulations.sims.keys())))
        self.assertEqual(simulations.sims[1].get_percentile(40), simulations.get_percentile(1, 40))
        self.assertEqual(simulations.sims[2].get_percentile(40), simulations.get_percentile(2, 40))
        # Check we get nearest allele length, when asking for an allele length that wasn't simulated
        self.assertEqual(simulations.sims[2].get_percentile(40), simulations.get_percentile(3, 40))

