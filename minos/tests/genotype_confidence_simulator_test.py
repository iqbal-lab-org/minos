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



