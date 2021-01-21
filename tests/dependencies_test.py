import os
import re
import unittest

from minos import dependencies


class TestDependencies(unittest.TestCase):
    def test_find_binary_gramtools(self):
        """test find_binary gramtools"""
        self.assertIsNotNone(dependencies.find_binary("gramtools"))

    def test_get_version_of_program_gramtools(self):
        """test get_version_of_program gramtools"""
        # We don't know what the version might be, so just
        # check that we got something that looks like X.Y.Z.
        got = dependencies.get_version_of_program("gramtools")
        self.assertIsNotNone(got)
        version_regex = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+$")
        self.assertIsNotNone(version_regex.search(got))

    def test_get_version_of_program_nextflow(self):
        """test get_version_of_program nextflow"""
        # We don't know what the version might be, so just
        # check that we got something that has numbers and dots
        got = dependencies.get_version_of_program("nextflow")
        self.assertIsNotNone(got)
        version_regex = re.compile(r"^[0-9\.]+ build [0-9]+$")
        self.assertIsNotNone(version_regex.search(got))

    def test_find_python_packages(self):
        """test find_python_packages"""
        # We don't know what the versions will be.
        # Just check that none of them are None
        got = dependencies.find_python_packages()
        for version, path in got.items():
            self.assertIsNotNone(version)
            self.assertIsNotNone(path)

    def test_find_binaries_and_versions(self):
        """test find_binaries_and_versions"""
        #  We don't know what the vetsions will be.
        # Just check we don't get None
        got = dependencies.find_binaries_and_versions()
        expect_keys = [
            "gramtools",
            "nextflow",
            "vcfallelicprimitives",
            "vcfbreakmulti",
            "vcfuniq",
            "vt",
        ]
        self.assertEqual(expect_keys, sorted(list(got.keys())))
        for version, path in got.items():
            self.assertIsNotNone(version)
            self.assertIsNotNone(path)

    def test_dependencies_report(self):
        """test dependencies_report"""
        got_ok, got_lines = dependencies.dependencies_report()
        self.assertTrue(got_ok)

        os.environ["MINOS_GRAMTOOLS"] = "oops_this_is_wrong"
        got_ok, got_lines = dependencies.dependencies_report()
        del os.environ["MINOS_GRAMTOOLS"]
        self.assertFalse(got_ok)

    def test_check_and_report_dependencies(self):
        """test check_and_report_dependencies"""
        tmpfile = "tmp.check_and_report_dependencies.out"
        if os.path.exists(tmpfile):
            os.unlink(tmpfile)
        dependencies.check_and_report_dependencies(outfile=tmpfile)
        # We don't inow what the versions etc will be, so just
        # check file got written
        self.assertTrue(os.path.exists(tmpfile))
        os.unlink(tmpfile)

        os.environ["MINOS_GRAMTOOLS"] = "oops_this_is_wrong"
        with self.assertRaises(Exception):
            dependencies.check_and_report_dependencies(outfile=tmpfile)
        del os.environ["MINOS_GRAMTOOLS"]
        self.assertTrue(os.path.exists(tmpfile))
        os.unlink(tmpfile)
