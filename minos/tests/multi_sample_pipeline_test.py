import filecmp
import os
import shutil
import unittest

from minos import multi_sample_pipeline

modules_dir = os.path.dirname(os.path.abspath(multi_sample_pipeline.py.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'multi_sample_pipeline.py')

class TestMultiSamplePipeline(unittest.TestCase):
    def test_run(self):
        '''test run'''
        self.assertTrue(False)

