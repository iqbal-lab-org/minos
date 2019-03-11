import logging
import numpy as np
from scipy import stats

from minos import genotyper


class GenotypeConfidenceSimulator:
    def __init__(self, mean_depth, depth_variance, error_rate, allele_length=1, iterations=10000):
        self.mean_depth = mean_depth
        self.depth_variance = depth_variance
        self.error_rate = error_rate
        self.iterations = iterations
        self.allele_length = allele_length
        self.confidence_scores_percentiles = {}
        self.min_conf_score = None
        self.max_conf_score = None


    @classmethod
    def _simulate_confidence_scores(cls, mean_depth, depth_variance, error_rate, iterations, allele_length=1, seed=42):
        np.random.seed(seed)
        allele_groups_dict = {'1': {0}, '2': {1}}
        i = 0
        confidences = []
        # We can't use the negative binomial unless depth_variance > mean_depth.
        # So force it to be so.
        if depth_variance < mean_depth:
            depth_variance = 2 * mean_depth
            logging.warn('Variance in read depth is smaller than mean read depth. Setting variance = 2 * mean, so that variant simulations can run. GT_CONF_PERCENTILE in the output VCF file may not be very useful as a result of this.')
        no_of_successes = (mean_depth ** 2) / (depth_variance - mean_depth)
        prob_of_success = 1 - (depth_variance - mean_depth) / depth_variance

        while i < iterations:
            correct_coverage = np.random.negative_binomial(no_of_successes, prob_of_success)
            incorrect_coverage = np.random.binomial(mean_depth, error_rate)
            if correct_coverage + incorrect_coverage == 0:
                continue

            allele_combination_cov = {}
            if incorrect_coverage > 0:
                allele_combination_cov['1'] = incorrect_coverage
            if correct_coverage > 0:
                allele_combination_cov['2'] = correct_coverage
            allele_per_base_cov = [[incorrect_coverage] * allele_length, [correct_coverage] * allele_length]
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
        confidence_scores = GenotypeConfidenceSimulator._simulate_confidence_scores(self.mean_depth, self.depth_variance, self.error_rate, self.iterations, allele_length=self.allele_length)
        self.confidence_scores_percentiles = GenotypeConfidenceSimulator._make_conf_to_percentile_dict(confidence_scores)



# This class is here for when we were simulating different allele lengths.
# It was used before the genotyper was changed to normalise non-zero coverage positions
# by allele length. It's no longer used, but leave it here just in case.
# Note that using it is currently pointless because the genotyper normalises by length,
# so would only make sense to use it if the genotyper was changed.
class Simulations:
    def __init__(self, mean_depth, depth_variance, error_rate, allele_lengths=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50], iterations=1000):
        self.sims = {l: GenotypeConfidenceSimulator(mean_depth, depth_variance, error_rate, allele_length=l, iterations=iterations) for l in allele_lengths}

        for allele_length in sorted(self.sims):
            logging.info(f'Run simulation, mean_depth={mean_depth}, depth_variance={depth_variance}, error_rate={error_rate}, allele_length={allele_length}, iterations={iterations}')
            self.sims[allele_length].run_simulations()
        logging.info('Finished running simulations')
        self.max_allele_length = max(self.sims.keys())
        self.nearest_allele_length = {}


    def get_percentile(self, allele_length, confidence):
        if allele_length in self.sims:
            return self.sims[allele_length].get_percentile(confidence)
        elif allele_length > self.max_allele_length:
            return self.sims[self.max_allele_length].get_percentile(confidence)
        elif allele_length in self.nearest_allele_length:
            return self.sims[self.nearest_allele_length[allele_length]].get_percentile(confidence)
        else:
            # Find the nearest allele length to the one we want, because
            # we haven't had this allele length before
            nearest_length = None

            for i in range(1, self.max_allele_length, 1):
                if (allele_length - i) in self.sims:
                    nearest_length = allele_length - i
                    break
                elif (allele_length + i) in self.sims:
                    nearest_length = allele_length + i
                    break

            assert nearest_length is not None

            self.nearest_allele_length[allele_length] = nearest_length
            return self.sims[nearest_length].get_percentile(confidence)


