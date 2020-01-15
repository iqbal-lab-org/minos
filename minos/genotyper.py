import itertools
import math
import operator

from scipy.stats import poisson


class Genotyper:
    def __init__(
        self,
        mean_depth,
        error_rate,
        allele_combination_cov,
        allele_per_base_cov,
        allele_groups_dict,
        min_cov_more_than_error=None,
    ):
        self.mean_depth = mean_depth
        self.error_rate = error_rate
        self.allele_combination_cov = allele_combination_cov
        self.allele_per_base_cov = allele_per_base_cov
        self.allele_groups_dict = allele_groups_dict
        self.likelihoods = None
        self.genotype = None
        self.genotype_confidence = None
        self.singleton_allele_coverages = {}
        self.haploid_allele_coverages = []
        if min_cov_more_than_error is None:
            self.min_cov_more_than_error = Genotyper.get_min_cov_to_be_more_likely_than_error(
                self.mean_depth, self.error_rate
            )
        else:
            self.min_cov_more_than_error = min_cov_more_than_error

    @classmethod
    def get_min_cov_to_be_more_likely_than_error(cls, mean_depth, error_rate):
        # For each allele, We want to get the number of positions that have read
        # depth higher than we would expect by chance. Work out the minimum
        # depth needed once here (code was previously repeatedly calling
        # poisson.pmf in _non_zeros_from_allele_per_base_cov()). Doing it once
        # here is about 80 to 90 times faster.

        # If the mean depth is zero, we'll be returning null genotypes with zero
        # confidence. Doesn't matter what we calculate here.
        if mean_depth == 0:
            return 0
        min_cov = 0
        while poisson.pmf(min_cov, mean_depth) <= pow(error_rate, min_cov):
            min_cov += 1
            if min_cov >= 1000:
                raise RuntimeError(
                    f"Something has likely gone wrong calculating minimum read depth that is more likely to happen than by chance due to read errors. error_rate={error_rate}, mean_depth={mean_depth}. Cannot continue"
                )
        return min_cov

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
        cls,
        allele1,
        allele2,
        allele_combination_cov,
        allele_groups_dict,
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
        if allele1_total_cov != 0 or allele2_total_cov != 0:  # If both are zero, there is no dispatching to do.
            allele1_belonging = allele1_total_cov / (allele1_total_cov + allele2_total_cov)
            allele1_total_cov += allele1_belonging * shared_cov
            allele2_total_cov += (1 - allele1_belonging) * shared_cov
        return allele1_total_cov, allele2_total_cov

    @classmethod
    def _log_likelihood_homozygous(
        cls, mean_depth, allele_depth, total_depth, error_rate, allele_length, non_zeros
    ):
        return sum(
            [
                -mean_depth * (1 + (allele_length - non_zeros) / allele_length),
                allele_depth * math.log(mean_depth),
                -math.lgamma(allele_depth + 1),
                (total_depth - allele_depth) * math.log(error_rate),
                non_zeros * math.log(1 - poisson.pmf(0, mean_depth)) / allele_length,
            ]
        )

    @classmethod
    def _log_likelihood_heterozygous(
        cls,
        mean_depth,
        allele_depth1,
        allele_depth2,
        total_depth,
        error_rate,
        allele_length1,
        allele_length2,
        non_zeros1,
        non_zeros2,
    ):
        return sum(
            [
                -mean_depth
                * (
                    1
                    + 0.5
                    * (
                        (1 - (non_zeros1 / allele_length1))
                        + (1 - (non_zeros2 / allele_length2))
                    )
                ),
                (allele_depth1 + allele_depth2) * math.log(0.5 * mean_depth),
                -math.lgamma(allele_depth1 + 1),
                -math.lgamma(allele_depth2 + 1),
                (total_depth - allele_depth1 - allele_depth2) * math.log(error_rate),
                ((non_zeros1 / allele_length1) + (non_zeros2 / allele_length2))
                * math.log(1 - poisson.pmf(0, 0.5 * mean_depth)),
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
            len(self.allele_per_base_cov), self.allele_combination_cov, self.allele_groups_dict
        )

        for allele_number, per_base_cov in enumerate(self.allele_per_base_cov):
            allele_length = len(per_base_cov)
            non_zeros = non_zeros_per_allele[allele_number]
            allele_depth = self.haploid_allele_coverages[allele_number]

            log_likelihood = Genotyper._log_likelihood_homozygous(
                self.mean_depth,
                allele_depth,
                total_depth,
                self.error_rate,
                allele_length,
                non_zeros,
            )
            self.likelihoods.append(({allele_number, allele_number}, log_likelihood))

        self.singleton_allele_coverages = Genotyper._singleton_alleles_and_coverage(
            self.allele_combination_cov, self.allele_groups_dict
        )

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

            log_likelihood = Genotyper._log_likelihood_heterozygous(
                self.mean_depth,
                allele1_depth,
                allele2_depth,
                total_depth,
                self.error_rate,
                allele1_length,
                allele2_length,
                non_zeros1,
                non_zeros2,
            )
            self.likelihoods.append(
                (set([allele_number1, allele_number2]), log_likelihood)
            )

        self.likelihoods.sort(key=operator.itemgetter(1), reverse=True)

    def run(self):
        if (
            len(self.allele_combination_cov) == 0
            or Genotyper._total_coverage(self.allele_combination_cov) == 0
            or self.mean_depth == 0
        ):
            self.genotype = {"."}
            self.genotype_confidence = 0.0
        else:
            self._calculate_log_likelihoods()
            assert self.likelihoods is not None and len(self.likelihoods) > 1
            self.genotype, best_log_likelihood = self.likelihoods[0]

            for allele in self.genotype:
                if self.singleton_allele_coverages.get(allele, 0) == 0:
                    self.genotype = {"."}
                    self.genotype_confidence = 0.0
                    break
            else:
                self.genotype_confidence = round(
                    best_log_likelihood - self.likelihoods[1][1], 2
                )
