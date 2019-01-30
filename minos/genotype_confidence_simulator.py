import numpy as np
from scipy import stats
import random

from minos import genotyper


class GenotypeConfidenceSimulator:
    def __init__(self, mean_depth, error_rate, iterations=10000):
        self.mean_depth = mean_depth
        self.error_rate = error_rate
        self.iterations = iterations
        self.confidence_scores_percentiles = {}


    @classmethod
    def _simulate_confidence_scores(cls, mean_depth, error_rate, iterations, seed=42):
        np.random.seed(seed)
        allele_groups_dict = {'1': {0}, '2': {1}}
        i = 0
        confidences = []

        while i < iterations:
            correct_coverage = np.random.poisson(mean_depth)
            incorrect_coverage = np.random.binomial(mean_depth, error_rate)
            if correct_coverage + incorrect_coverage == 0:
                continue

            allele_combination_cov = {'1': correct_coverage, '2': incorrect_coverage}
            allele_per_base_cov = [[correct_coverage], [incorrect_coverage]]
            gtyper = genotyper.Genotyper(mean_depth, error_rate, allele_combination_cov, allele_per_base_cov, allele_groups_dict)
            gtyper.run()
            confidences.append(round(gtyper.genotype_confidence))
            i += 1

        assert len(confidences) == iterations
        confidences.sort()
        return confidences


    @classmethod
    def _make_conf_to_percentile_dict(cls, confidence_scores):
        assert len(confidence_scores) > 0
        confidence_scores_percentiles = 100 * stats.rankdata(confidence_scores) / len(confidence_scores)
        conf_to_percentile = {}
        for conf, percentile in zip(confidence_scores, confidence_scores_percentiles):
            conf_to_percentile[conf] = round(percentile)
        return conf_to_percentile


    def run_simulations(self):
        confidence_scores = GenotypeConfidenceSimulator._simulate_confidence_scores(self.mean_depth, self.error_rate, self.iterations)
        self.confidence_scores_percentiles = GenotypeConfidenceSimulator._make_conf_to_percentile_dict(confidence_scores)

