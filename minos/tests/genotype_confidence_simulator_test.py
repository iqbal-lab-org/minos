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
        expected = [6, 9, 15, 21, 31]
        self.assertEqual(expected, got)


    def test_make_conf_to_percentile_dict(self):
        '''test _make_conf_to_percentile_dict'''
        confidence_scores = [1, 1, 2, 3, 4, 4, 4, 5, 6, 8]
        got = genotype_confidence_simulator.GenotypeConfidenceSimulator._make_conf_to_percentile_dict(confidence_scores)
        expected = {1: 15, 2: 30, 3: 40, 4: 50, 5: 80, 6: 90, 8: 100}
        self.assertTrue(expected, got)


    def test_run_simulations_and_get_percentile(self):
        '''test run_simulations and get_percentile'''
        simulator = genotype_confidence_simulator.GenotypeConfidenceSimulator(50, 0.1, 5)
        simulator.run_simulations()
        expected_confidence_scores_percentiles = {6: 20.00, 9: 40.00, 15: 60.00, 21: 80.00, 31: 100.00}
        self.assertEqual(expected_confidence_scores_percentiles, simulator.confidence_scores_percentiles)
        self.assertEqual(20.00, simulator.get_percentile(6))
        self.assertEqual(40.00, simulator.get_percentile(9))
        # Try getting numbers that are not in the dict and will have to be inferred
        self.assertEqual(26.67, simulator.get_percentile(7))
        self.assertEqual(90.00, simulator.get_percentile(26))
        self.assertEqual(92.00, simulator.get_percentile(27))
        # Try values outside the range of what we already have
        self.assertEqual(0.00, simulator.get_percentile(4))
        self.assertEqual(0.00, simulator.get_percentile(5))
        self.assertEqual(100.00, simulator.get_percentile(32))
        self.assertEqual(100.00, simulator.get_percentile(33))

