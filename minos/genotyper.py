import itertools
import logging
import math
import operator

from scipy.stats import nbinom


class Genotyper:
    def __init__(
        self,
        mean_depth,
        depth_variance,
        error_rate,
        call_hets=False,
    ):
        self.mean_depth = mean_depth
        self.depth_variance = depth_variance
        self.error_rate = error_rate
        self.call_hets = call_hets
        if self.call_hets:
            raise NotImplementedError("Heterozygous calling is not implemented")
        self._set_nbinom_parameters()
        self._set_min_cov_more_than_error()
        self._init_alleles_and_genotypes()
        logging.info(
            f"Genotyper parameters: mean_depth={mean_depth}, depth_variance={depth_variance}, error_rate={error_rate}, nbinom_no_of_successes={self.no_of_successes}, nbinom_prob_of_success={self.prob_of_success}, min_cov_more_than_error={self.min_cov_more_than_error}"
        )

    def _set_nbinom_parameters(self):
        # We can only use nbinom if variance > depth. Otherwise, force it by
        # setting variance = 2 * depth. This is very rare - only seen it on
        # simulated data, never on real data.
        if self.depth_variance <= self.mean_depth:
            logging.warning(f"Mean read depth ({self.mean_depth}) > depth variance ({self.depth_variance}) . Setting variance = 2 * mean for genotype model")
            self.depth_variance = 2 * self.mean_depth

        if self.mean_depth == 0:
            self.no_of_successes = 0
            self.prob_of_success = 0
        else:
            self.no_of_successes = (self.mean_depth ** 2) / (
                self.depth_variance - self.mean_depth
            )
            self.prob_of_success = (
                1 - (self.depth_variance - self.mean_depth) / self.depth_variance
            )

    def _nbinom_pmf(self, value, log=False):
        if log:
            return nbinom.logpmf(value, self.no_of_successes, self.prob_of_success)
        else:
            return nbinom.pmf(value, self.no_of_successes, self.prob_of_success)

    def _set_min_cov_more_than_error(self):
        self.min_cov_more_than_error = 0
        if self.mean_depth == 0:
            return
        while self._nbinom_pmf(self.min_cov_more_than_error) <= pow(
            self.error_rate, self.min_cov_more_than_error
        ):
            self.min_cov_more_than_error += 1
            if self.min_cov_more_than_error >= 1000:
                raise RuntimeError(
                    f"Something has likely gone wrong calculating minimum read depth that is more likely to happen than by chance due to read errors. error_rate={self.error_rate}, mean_depth={self.mean_depth}, depth_variance={self.depth_variance}. Cannot continue"
                )

    def _init_alleles_and_genotypes(
        self,
        allele_combination_cov=None,
        allele_per_base_cov=None,
        allele_groups_dict=None,
    ):
        self.allele_combination_cov = allele_combination_cov
        self.allele_per_base_cov = allele_per_base_cov
        self.allele_groups_dict = allele_groups_dict
        self.likelihoods = None
        self.genotype = None
        self.genotype_frs = None
        self.genotype_confidence = None
        self.singleton_allele_coverages = {}
        self.haploid_allele_coverages = []

    @classmethod
    def _singleton_alleles_and_coverage(
        cls, allele_combination_cov, allele_groups_dict
    ):
        singleton_alleles = {}
        for allele_key in allele_combination_cov:
            if len(allele_groups_dict[allele_key]) == 1:
                allele = next(iter(allele_groups_dict[allele_key]))
                singleton_alleles[allele] = allele_combination_cov[allele_key]
        return singleton_alleles

    @classmethod
    def _total_coverage(cls, allele_combination_cov):
        return sum(allele_combination_cov.values())

    @classmethod
    def _haploid_allele_coverages(
        cls, num_distinct_alleles, allele_combination_cov, allele_groups_dict
    ):
        haploid_allele_coverages = [0 for i in range(num_distinct_alleles)]
        for allele_key, coverage in allele_combination_cov.items():
            alleles = list(iter(allele_groups_dict[allele_key]))
            for allele in alleles:
                haploid_allele_coverages[allele] += coverage
        return haploid_allele_coverages

    def _log_likelihood_homozygous(
        self, allele_depth, total_depth, allele_length, non_zeros
    ):
        return sum([
                self._nbinom_pmf(allele_depth, log=True),
                (total_depth - allele_depth) * math.log(self.error_rate),
                math.log(1 - self._nbinom_pmf(0)) * non_zeros / allele_length,
                self._nbinom_pmf(0, log=True) * (allele_length - non_zeros) / allele_length,
            ]
        )

    @classmethod
    def _non_zeros_from_allele_per_base_cov(cls, allele_per_base_cov, min_cov):
        non_zeros = []
        for per_base_cov in allele_per_base_cov:
            non_zero_count = len([c for c in per_base_cov if c >= min_cov])
            non_zeros.append(non_zero_count)
        return non_zeros

    def _calculate_log_likelihoods(self):
        """Makes a list of tuples: ( (allele(s) tuple), log likelihood).
        List is sorted from most to least likely"""
        self.likelihoods = []
        total_depth = sum(self.allele_combination_cov.values())

        non_zeros_per_allele = Genotyper._non_zeros_from_allele_per_base_cov(
            self.allele_per_base_cov, self.min_cov_more_than_error
        )

        self.haploid_allele_coverages = Genotyper._haploid_allele_coverages(
            len(self.allele_per_base_cov),
            self.allele_combination_cov,
            self.allele_groups_dict,
        )

        for allele_number, per_base_cov in enumerate(self.allele_per_base_cov):
            allele_length = len(per_base_cov)
            non_zeros = non_zeros_per_allele[allele_number]
            allele_depth = self.haploid_allele_coverages[allele_number]

            log_likelihood = self._log_likelihood_homozygous(
                allele_depth, total_depth, allele_length, non_zeros,
            )
            frac_support = (
                0 if total_depth == 0 else round(allele_depth / total_depth, 4)
            )
            self.likelihoods.append(
                ({allele_number, allele_number}, log_likelihood, frac_support)
            )

        self.singleton_allele_coverages = Genotyper._singleton_alleles_and_coverage(
            self.allele_combination_cov, self.allele_groups_dict
        )

        self.likelihoods.sort(key=operator.itemgetter(1), reverse=True)
        logging.debug(f"likelihoods: {self.likelihoods}")


    def run(self, allele_combination_cov, allele_per_base_cov, allele_groups_dict):
        self._init_alleles_and_genotypes(
            allele_combination_cov=allele_combination_cov,
            allele_per_base_cov=allele_per_base_cov,
            allele_groups_dict=allele_groups_dict,
        )

        if (
            len(self.allele_combination_cov) == 0
            or Genotyper._total_coverage(self.allele_combination_cov) == 0
            or self.mean_depth == 0
        ):
            self.genotype = {"."}
            self.genotype_confidence = 0.0
            self.genotype_frs = "."
        else:
            self._calculate_log_likelihoods()
            assert self.likelihoods is not None and len(self.likelihoods) > 1
            self.genotype, best_log_likelihood, self.genotype_frs = self.likelihoods[0]

            for allele in self.genotype:
                if self.singleton_allele_coverages.get(allele, 0) == 0:
                    self.genotype = {"."}
                    self.genotype_confidence = 0.0
                    self.genotype_frs = 0
                    break
            else:
                self.genotype_confidence = round(
                    best_log_likelihood - self.likelihoods[1][1], 2
                )
