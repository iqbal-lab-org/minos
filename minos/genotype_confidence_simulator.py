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
        self.min_conf_score = None
        self.max_conf_score = None


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

            allele_combination_cov = {}
            if incorrect_coverage > 0:
                allele_combination_cov['1'] = incorrect_coverage
            if correct_coverage > 0:
                allele_combination_cov['2'] = correct_coverage
            allele_per_base_cov = [[incorrect_coverage], [correct_coverage]]
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
            conf_to_percentile[conf] = round(percentile, 2)
        return conf_to_percentile


    def get_percentile(self, confidence):
        percentile = self.confidence_scores_percentiles.get(confidence, None)
        if percentile is not None:
            return percentile

        # If we don't know this percentile, then infer it by linearly interpolating
        # between the two nearest values, and add
        # it to the dict, so we don't have to work it out again
        if self.min_conf_score is None or self.max_conf_score is None:
            self.min_conf_score = min(self.confidence_scores_percentiles.keys())
            self.max_conf_score = max(self.confidence_scores_percentiles.keys())

        if confidence < self.min_conf_score:
            self.confidence_scores_percentiles[confidence] = 0.00
        elif confidence > self.max_conf_score:
            self.confidence_scores_percentiles[confidence] = 100.00
        else:
            i = 0
            left_conf = None
            right_conf = None

            for i in range(1, self.max_conf_score, 1):
                if left_conf is None:
                    left_conf = confidence - i if confidence - i in self.confidence_scores_percentiles else None
                if right_conf is None:
                    right_conf = self.confidence_scores_percentiles.get(confidence + 1, None)
                    right_conf = confidence + i if confidence + i in self.confidence_scores_percentiles else None
                if left_conf is not None and right_conf is not None:
                    break
            else:
                raise Exception(f'Error inferring new confidence score {confidence}. Cannot continue')

            assert self.min_conf_score <=  left_conf < confidence < right_conf <= self.max_conf_score
            left_pc = self.confidence_scores_percentiles[left_conf]
            right_pc = self.confidence_scores_percentiles[right_conf]
            if left_pc == right_pc:
                self.confidence_scores_percentiles[confidence] = left_pc
            else:
                self.confidence_scores_percentiles[confidence] = round(left_pc + ((confidence - left_conf) / (right_conf - left_conf)) * (right_pc - left_pc), 2)

        return self.confidence_scores_percentiles[confidence]


    def run_simulations(self):
        confidence_scores = GenotypeConfidenceSimulator._simulate_confidence_scores(self.mean_depth, self.error_rate, self.iterations)
        self.confidence_scores_percentiles = GenotypeConfidenceSimulator._make_conf_to_percentile_dict(confidence_scores)

