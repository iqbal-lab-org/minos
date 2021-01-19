import itertools
import logging
import math
import operator

from scipy.stats import nbinom, poisson


class Genotyper:
    def __init__(
        self,
        mean_depth,
        depth_variance,
        error_rate,
        call_hets=False,
        force_poisson=False,
    ):
        self.mean_depth = mean_depth
        self.depth_variance = depth_variance
        self.error_rate = error_rate
        self.call_hets = call_hets
        self.force_poisson = force_poisson
        self._determine_poisson_or_nbinom()
        self._set_min_cov_more_than_error()
        self._init_alleles_and_genotypes()
        logging.info(
            f"Genotyper. use_het_model={call_hets}, mean_depth={mean_depth}, depth_variance={depth_variance}, error_rate={error_rate}, using_nbinom={self.use_nbinom}, nbinom_no_of_successes={self.no_of_successes}, nbinom_prob_of_success={self.prob_of_success}, min_cov_more_than_error={self.min_cov_more_than_error}"
        )

    def _determine_poisson_or_nbinom(self):
        # We can only use nbinom if variance > depth. Otherwise fall back on
        # poisson model, which is likely to be worse. Should only happen when
        # read depth is very low
        if self.force_poisson:
            self.use_nbinom = False
        else:
            self.use_nbinom = self.depth_variance > self.mean_depth

        if self.use_nbinom:
            self.no_of_successes = (self.mean_depth ** 2) / (
                self.depth_variance - self.mean_depth
            )
            self.prob_of_success = (
                1 - (self.depth_variance - self.mean_depth) / self.depth_variance
            )
            self.half_no_of_successes = ((0.5 * self.mean_depth) ** 2) / (
                self.depth_variance - (0.5 * self.mean_depth)
            )
            self.half_prob_of_success = (
                1
                - (self.depth_variance - (0.5 * self.mean_depth)) / self.depth_variance
            )
        else:
            self.no_of_successes = None
            self.prob_of_success = None
            self.half_no_of_successes = None
            self.half_prob_of_success = None

    def _nbinom_or_poisson_pmf(self, value, half_mean=False):
        if self.use_nbinom:
            if half_mean:
                return nbinom.pmf(
                    value, self.half_no_of_successes, self.half_prob_of_success
                )
            else:
                return nbinom.pmf(value, self.no_of_successes, self.prob_of_success)
        else:
            if half_mean:
                return poisson.pmf(value, 0.5 * self.mean_depth)
            else:
                return poisson.pmf(value, self.mean_depth)

    def _set_min_cov_more_than_error(self):
        self.min_cov_more_than_error = 0
        if self.mean_depth == 0:
            return
        while self._nbinom_or_poisson_pmf(self.min_cov_more_than_error) <= pow(
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

    @classmethod
    def _coverage_of_diploid_alleles(
        cls, allele1, allele2, allele_combination_cov, allele_groups_dict,
    ):
        """
        Collects coverage on each of two alleles.
        Define specific coverage for one allele as the sum of coverages for equivalence classes involving that
        allele but not the other one.

        When an equivalence class contains both alleles, the coverage on it is dispatched
        to each allele proportionally to how much specific coverage each has.
        """
        shared_cov = 0
        allele1_total_cov, allele2_total_cov = 0, 0

        for allele_key, coverage in allele_combination_cov.items():
            assert coverage >= 0
            allele_combination = allele_groups_dict[allele_key]
            has_allele_1 = allele1 in allele_combination
            has_allele_2 = allele2 in allele_combination
            if has_allele_1 and has_allele_2:
                shared_cov += coverage
            elif has_allele_1:
                allele1_total_cov += coverage
            elif has_allele_2:
                allele2_total_cov += coverage

        ## Perform the dispatching in ambiguous equiv classes ##
        if (
            allele1_total_cov == 0 and allele2_total_cov == 0
        ):  # If both are zero, do equal dispatching
            allele1_belonging = 0.5
        else:
            allele1_belonging = allele1_total_cov / (
                allele1_total_cov + allele2_total_cov
            )

        allele1_total_cov += allele1_belonging * shared_cov
        allele2_total_cov += (1 - allele1_belonging) * shared_cov
        return allele1_total_cov, allele2_total_cov

    def _log_likelihood_homozygous(
        self, allele_depth, total_depth, allele_length, non_zeros
    ):
        p_nonzero = 1 - self._nbinom_or_poisson_pmf(0)
        return sum([
                math.log(self._nbinom_or_poisson_pmf(allele_depth)),
                (total_depth - allele_depth) * math.log(self.error_rate),
                math.log(p_nonzero) * non_zeros / allele_length,
                math.log(1-p_nonzero) * (allele_length - non_zeros) / allele_length,
            ]
        )

    def _log_likelihood_heterozygous(
        self,
        allele_depth1,
        allele_depth2,
        total_depth,
        allele_length1,
        allele_length2,
        non_zeros1,
        non_zeros2,
    ):
        return sum(
            [
                -self.mean_depth
                * (
                    1
                    + 0.5
                    * (
                        (1 - (non_zeros1 / allele_length1))
                        + (1 - (non_zeros2 / allele_length2))
                    )
                ),
                (allele_depth1 + allele_depth2) * math.log(0.5 * self.mean_depth),
                -math.lgamma(allele_depth1 + 1),
                -math.lgamma(allele_depth2 + 1),
                (total_depth - allele_depth1 - allele_depth2)
                * math.log(self.error_rate),
                ((non_zeros1 / allele_length1) + (non_zeros2 / allele_length2))
                * math.log(1 - self._nbinom_or_poisson_pmf(0, half_mean=True)),
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

        if not self.call_hets:
            self.likelihoods.sort(key=operator.itemgetter(1), reverse=True)
            logging.debug(f"likelihoods: {self.likelihoods}")
            return

        for (allele_number1, allele_number2) in itertools.combinations(
            self.singleton_allele_coverages.keys(), 2
        ):
            allele1_depth, allele2_depth = Genotyper._coverage_of_diploid_alleles(
                allele_number1,
                allele_number2,
                self.allele_combination_cov,
                self.allele_groups_dict,
            )
            allele1_length = len(self.allele_per_base_cov[allele_number1])
            allele2_length = len(self.allele_per_base_cov[allele_number2])
            non_zeros1 = non_zeros_per_allele[allele_number1]
            non_zeros2 = non_zeros_per_allele[allele_number2]

            log_likelihood = self._log_likelihood_heterozygous(
                allele1_depth,
                allele2_depth,
                total_depth,
                allele1_length,
                allele2_length,
                non_zeros1,
                non_zeros2,
            )
            frac_support = (
                0
                if total_depth == 0
                else round((allele1_depth + allele2_depth) / total_depth, 4)
            )
            self.likelihoods.append(
                (set([allele_number1, allele_number2]), log_likelihood, frac_support)
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
