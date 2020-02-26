import os
import unittest

from minos import genotyper

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "genotyper")


class TestGenotyper(unittest.TestCase):
    def test_get_min_cov_to_be_more_likely_than_error(self):
        """test get_min_cov_to_be_more_likely_than_error"""
        self.assertEqual(
            1, genotyper.Genotyper.get_min_cov_to_be_more_likely_than_error(10, 0.0001)
        )
        self.assertEqual(
            2, genotyper.Genotyper.get_min_cov_to_be_more_likely_than_error(10, 0.001)
        )
        self.assertEqual(
            10, genotyper.Genotyper.get_min_cov_to_be_more_likely_than_error(100, 0.001)
        )
        self.assertEqual(
            0, genotyper.Genotyper.get_min_cov_to_be_more_likely_than_error(0, 0.001)
        )
        with self.assertRaises(RuntimeError):
            genotyper.Genotyper.get_min_cov_to_be_more_likely_than_error(10000, 0.5)

    def test_singleton_alleles_and_coverage(self):
        """test _singleton_alleles_and_coverage"""
        allele_combination_cov = {"1": 20, "3": 1}
        allele_groups_dict = {"1": {0}, "2": {1}, "3": {1, 2}, "4": {5, 6}}
        self.assertEqual(
            {0: 20},
            genotyper.Genotyper._singleton_alleles_and_coverage(
                allele_combination_cov, allele_groups_dict
            ),
        )
        allele_combination_cov["2"] = 42
        self.assertEqual(
            {0: 20, 1: 42},
            genotyper.Genotyper._singleton_alleles_and_coverage(
                allele_combination_cov, allele_groups_dict
            ),
        )

    def test_total_coverage(self):
        """test _total_coverage"""
        self.assertEqual(0, genotyper.Genotyper._total_coverage({}))
        self.assertEqual(1, genotyper.Genotyper._total_coverage({"x": 1}))
        self.assertEqual(42, genotyper.Genotyper._total_coverage({"x": 1, "y": 41}))


    def test_haploid_allele_coverages(self):
        """test  _haploid_allele_coverages"""
        allele_combination_cov = {"1": 20, "2": 1}
        allele_groups_dict = {"0": {0}, "1": {1}, "2": {1, 2}, "3": {5, 6}}
        num_distinct_alleles = 7 # 1 + the max allele index
        self.assertEqual(
            [0, 21, 1, 0, 0, 0, 0],
            genotyper.Genotyper._haploid_allele_coverages(
                7, allele_combination_cov, allele_groups_dict
            ),
        )

    def test_coverage_of_diploid_alleles_equal_dispatching(self):
        """
        test _coverage_of_diploid_alleles function
        # Equiv class -> coverage:
        #   {0,1}: 9
        We collect coverages of alleles 0 and 1.
        This is an edge case where the only coverage is on both: then we dispatch equally.
        """
        allele_combination_cov = {"1": 9}
        allele_groups_dict = {"1": {0, 1}}
        self.assertEqual(
            (4.5, 4.5),
            genotyper.Genotyper._coverage_of_diploid_alleles(
                0, 1, allele_combination_cov, allele_groups_dict,
            ),
        )

    def test_coverage_of_diploid_alleles_correct_dispatching(self):
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
        self.assertEqual(
            (22, 88),
            genotyper.Genotyper._coverage_of_diploid_alleles(
                0, 1, allele_combination_cov, allele_groups_dict,
            ),
        )

    def test_log_likelihood_homozygous(self):
        """test _log_likelihood_homozygous"""
        self.assertEqual(
            -26.71,
            round(
                genotyper.Genotyper._log_likelihood_homozygous(100, 90, 95, 0.01, 5, 5),
                2,
            ),
        )
        self.assertEqual(
            -44.54,
            round(
                genotyper.Genotyper._log_likelihood_homozygous(10, 1, 9, 0.01, 5, 5), 2
            ),
        )

    def test_log_likelihood_heterozygous(self):
        """test _log_likelihood_heterozygous"""
        # _log_likelihood_heterozygous(cls, mean_depth, allele_depth1, allele_depth2, total_depth,
        #    error_rate, allele_length1, allele_length2, non_zeros1, non_zeros2):
        self.assertEqual(
            -52.97,
            round(
                genotyper.Genotyper._log_likelihood_heterozygous(
                    100, 45, 40, 95, 0.01, 3, 3, 3, 3
                ),
                2,
            ),
        )
        self.assertEqual(
            -86.31,
            round(
                genotyper.Genotyper._log_likelihood_heterozygous(
                    100, 45, 40, 95, 0.01, 3, 3, 2, 2
                ),
                2,
            ),
        )

    def test_calculate_log_likelihoods(self):
        """test _calculate_log_likelihoods"""
        mean_depth = 20
        error_rate = 0.01
        allele_combination_cov = {"1": 2, "2": 20, "3": 1}
        allele_groups_dict = {"1": {0}, "2": {1}, "3": {0, 1}, "4": {2}}
        allele_per_base_cov = [[0, 1], [20, 19]]
        gtyper = genotyper.Genotyper(
            mean_depth,
            error_rate,
            allele_combination_cov,
            allele_per_base_cov,
            allele_groups_dict,
        )
        gtyper._calculate_log_likelihoods()
        expected = [({1}, -11.68), ({0, 1}, -22.92), ({0}, -124.91)]
        self.assertEqual(3, len(gtyper.likelihoods))
        gtyper.likelihoods = [(x[0], round(x[1], 2)) for x in gtyper.likelihoods]
        self.assertEqual(expected, gtyper.likelihoods)

    def test_run_with_call_hets_false(self):
        """test run"""
        mean_depth = 20
        error_rate = 0.01
        allele_combination_cov = {"1": 2, "2": 20, "3": 1}
        allele_groups_dict = {"1": {0}, "2": {1}, "3": {0, 1}, "4": {2}}
        allele_per_base_cov = [[0, 1], [20, 19]]
        gtyper = genotyper.Genotyper(
            mean_depth,
            error_rate,
            allele_combination_cov,
            allele_per_base_cov,
            allele_groups_dict,
            call_hets=False,
        )
        expected = [({1}, -11.68), ({0}, -124.91)]
        gtyper.run()
        self.assertEqual(len(expected), len(gtyper.likelihoods))
        for i in range(len(expected)):
            self.assertEqual(expected[i][0], gtyper.likelihoods[i][0])
            self.assertAlmostEqual(expected[i][1], gtyper.likelihoods[i][1], places=2)

    def test_run(self):
        """test run"""
        mean_depth = 20
        error_rate = 0.01
        allele_combination_cov = {"1": 2, "2": 20, "3": 1}
        allele_groups_dict = {"1": {0}, "2": {1}, "3": {0, 1}, "4": {2}}
        allele_per_base_cov = [[0, 1], [20, 19]]
        gtyper = genotyper.Genotyper(
            mean_depth,
            error_rate,
            allele_combination_cov,
            allele_per_base_cov,
            allele_groups_dict,
        )
        expected = [({1}, -11.68), ({0, 1}, -22.92), ({0}, -124.91)]
        gtyper.run()
        self.assertEqual(len(expected), len(gtyper.likelihoods))
        for i in range(len(expected)):
            self.assertEqual(expected[i][0], gtyper.likelihoods[i][0])
            self.assertAlmostEqual(expected[i][1], gtyper.likelihoods[i][1], places=2)

    def test_run_zero_coverage(self):
        """test run when all alleles have zero coverage"""
        mean_depth = 20
        error_rate = 0.01
        allele_combination_cov = {}
        allele_groups_dict = {"1": {0}, "2": {1}, "3": {0, 1}, "4": {2}}
        allele_per_base_cov = [[0], [0, 0]]
        gtyper = genotyper.Genotyper(
            mean_depth,
            error_rate,
            allele_combination_cov,
            allele_per_base_cov,
            allele_groups_dict,
        )
        gtyper.run()
        self.assertEqual({"."}, gtyper.genotype)
        self.assertEqual(0.0, gtyper.genotype_confidence)

    def test_nomatherror_mean_depth0(self):
        """
        Can get a mean_depth of zero but try to genotype a non-zero coverage site due to rounding imprecision.
        In which case we need to avoid trying to do log(0) in likelihood calculation and should return no call.
        """
        mean_depth = 0
        error_rate = 0.01
        allele_combination_cov = {"1": 1}
        allele_groups_dict = {"1": {0}, "2": {1}}
        allele_per_base_cov = [[1], [0, 0]]
        gtyper = genotyper.Genotyper(
            mean_depth,
            error_rate,
            allele_combination_cov,
            allele_per_base_cov,
            allele_groups_dict,
        )
        gtyper.run()
        self.assertEqual({"."}, gtyper.genotype)
        self.assertEqual(0.0, gtyper.genotype_confidence)
