import math
import unittest

import numpy as np

from src.complex_phase import BASES
from src.gap_geometry import (
    TOP_K_VALUES,
    geometry_records,
    load_inputs,
    midpoint_angle_deg,
    polygon_metrics,
    select_largest_gaps,
)


class GapGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.gaps, _ = load_inputs()
        cls.summary, cls.selections = geometry_records(cls.gaps)

    def test_midpoint_angles_are_in_expected_interval(self) -> None:
        self.assertEqual(midpoint_angle_deg(350.0, 30.0), 5.0)
        for selected in self.selections.values():
            self.assertTrue(selected["mid_theta_deg"].ge(0.0).all())
            self.assertTrue(selected["mid_theta_deg"].lt(360.0).all())

    def test_selected_gap_count_matches_requested_k_when_available(self) -> None:
        for base in BASES:
            available = len(self.gaps[self.gaps["base"].astype(str) == base])
            for top_k in TOP_K_VALUES:
                selected = select_largest_gaps(self.gaps, base, top_k)
                self.assertEqual(len(selected), min(top_k, available))

    def test_polygon_area_is_non_negative(self) -> None:
        self.assertTrue(self.summary["polygon_area"].ge(0.0).all())

    def test_symmetry_score_is_finite_and_positive(self) -> None:
        scores = self.summary["symmetry_score"].to_numpy()
        self.assertTrue(np.isfinite(scores).all())
        self.assertTrue((scores > 0.0).all())

    def test_angular_gaps_sum_to_full_circle(self) -> None:
        for selected in self.selections.values():
            metrics = polygon_metrics(selected["mid_theta_deg"].to_numpy())
            self.assertTrue(
                math.isclose(
                    sum(metrics["angle_gaps_between_vertices"]),
                    360.0,
                    rel_tol=1e-12,
                    abs_tol=1e-12,
                )
            )


if __name__ == "__main__":
    unittest.main()
