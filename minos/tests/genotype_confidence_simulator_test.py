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


    def test_run_simulations(self):
        '''test run_simulations'''
        simulator = genotype_confidence_simulator.GenotypeConfidenceSimulator(50, 0.1, 5)
        simulator.run_simulations()
        expected_confidence_scores_percentiles = {6: 20.0, 9: 40.0, 15: 60.0, 21: 80.0, 31: 100.0}
        self.assertEqual(expected_confidence_scores_percentiles, simulator.confidence_scores_percentiles)
