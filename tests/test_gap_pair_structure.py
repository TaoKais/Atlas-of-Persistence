import unittest

import numpy as np

from src.gap_geometry import load_inputs
from src.gap_pair_structure import (
    RATIO_METRICS,
    pair_structure_table,
    ratio_summary_table,
    selected_gap_pairs,
)


class GapPairStructureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.gaps, cls.phases = load_inputs()
        cls.selected = selected_gap_pairs(cls.gaps)
        cls.pairs = pair_structure_table(cls.gaps, cls.phases)
        cls.summary = ratio_summary_table(cls.pairs)

    def test_exactly_five_gap_pairs_are_selected(self) -> None:
        self.assertEqual(len(self.selected), 5)
        self.assertEqual(len(self.pairs), 5)

    def test_selected_pairs_are_the_largest_base_two_gaps(self) -> None:
        expected = (
            self.gaps[self.gaps["base"].astype(str) == "2"]
            .nlargest(5, "gap_deg")["gap_deg"]
            .to_numpy()
        )
        np.testing.assert_allclose(self.selected["gap_deg"].to_numpy(), expected)

    def test_ratios_and_geometric_means_are_positive(self) -> None:
        for metric in RATIO_METRICS:
            self.assertTrue(self.pairs[f"ratio_{metric}"].gt(0.0).all())
            self.assertTrue(self.pairs[f"geometric_mean_{metric}"].gt(0.0).all())

    def test_logarithms_match_cycle_counts(self) -> None:
        np.testing.assert_allclose(
            self.pairs["from_log2_N"], np.log2(self.pairs["from_compton_cycles"])
        )
        np.testing.assert_allclose(
            self.pairs["to_log10_N"], np.log10(self.pairs["to_compton_cycles"])
        )

    def test_ratio_summary_is_finite(self) -> None:
        numeric = self.summary.drop(columns="metric").drop(columns="shared_scale")
        self.assertTrue(np.isfinite(numeric.to_numpy(dtype=float)).all())
        self.assertTrue(self.summary["multiplicative_spread"].ge(1.0).all())


if __name__ == "__main__":
    unittest.main()
