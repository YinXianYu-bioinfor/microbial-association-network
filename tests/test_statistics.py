import unittest
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from network_core import rank_correlation, bh_adjust, exact_permutation_p, exact_pairwise_p


class StatisticsTests(unittest.TestCase):
    def test_exact_pairwise_matches_enumeration_with_ties(self):
        values = np.array([[1, 1, 2, 3, 4], [4, 2, 2, 3, 1], [1, 3, 2, 5, 4]], dtype=float)
        pairs = np.triu_indices(3, 1)
        actual = exact_pairwise_p(values, *pairs)
        expected = [exact_permutation_p(values[i], values[j]) for i, j in zip(*pairs)]
        np.testing.assert_allclose(actual, expected)

    def test_tied_spearman_and_constant_feature(self):
        r = rank_correlation(np.array([[1, 1, 3, 3], [1, 2, 3, 4], [1, 1, 1, 1]], dtype=float))
        self.assertAlmostEqual(r[0, 1], 0.8944271909999159)
        self.assertTrue(np.isnan(r[0, 2]))

    def test_exact_permutation_and_bh(self):
        self.assertAlmostEqual(exact_permutation_p([1, 2, 3, 4], [4, 3, 2, 1]), 2/24)
        np.testing.assert_allclose(bh_adjust([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])
