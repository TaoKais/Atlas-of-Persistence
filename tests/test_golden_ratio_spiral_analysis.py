import math
import unittest

import numpy as np

from src.golden_ratio_spiral_analysis import (
    GOLDEN_K,
    PHI,
    PHI_NEIGHBORS,
    SPECIAL_BASES,
    base_metrics,
    control_metrics,
    load_entities,
    radial_coordinates,
    spiral_fit,
)


class GoldenRatioSpiralAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.entities = load_entities()
        cls.frame = radial_coordinates(cls.entities, PHI)

    def test_finite_N(self) -> None:
        self.assertTrue(np.isfinite(self.entities["compton_cycles"]).all())
        self.assertTrue(self.entities["compton_cycles"].gt(0.0).all())

    def test_theta_is_wrapped(self) -> None:
        self.assertTrue(self.frame["theta"].ge(0.0).all())
        self.assertTrue(self.frame["theta"].lt(2.0 * math.pi).all())

    def test_log_spiral_fit_uses_positive_radius(self) -> None:
        fit = spiral_fit(self.frame)
        self.assertTrue(np.isfinite(fit["log_r"]).all())
        self.assertGreater(len(fit["log_r"]), 2)

    def test_metrics_have_no_nan(self) -> None:
        metrics = base_metrics(self.entities, PHI)
        self.assertTrue(all(np.isfinite(value) for value in metrics.values()))

    def test_controls_run_successfully(self) -> None:
        controls = control_metrics(self.entities, replicates=2)
        self.assertEqual(set(controls["control"]), {
            "observed", "shuffle_lifetimes", "shuffle_masses", "shuffle_N_values",
            "random_logN_distribution", "random_angular_phase",
        })
        self.assertFalse(controls["geometry_coherence_score"].isna().any())

    def test_phi_is_included_exactly(self) -> None:
        self.assertEqual(SPECIAL_BASES["phi"], PHI)
        self.assertTrue(any(value == PHI for value in PHI_NEIGHBORS))

    def test_golden_k_formula(self) -> None:
        self.assertEqual(GOLDEN_K, 2.0 * math.log(PHI) / math.pi)

    def test_infinity_marker_mode_preserves_finite_fit_rows(self) -> None:
        entities = load_entities("stable_as_infinity_marker")
        self.assertIn("electron", entities.attrs["infinity_markers"])
        self.assertTrue(np.isfinite(entities["compton_cycles"]).all())


if __name__ == "__main__":
    unittest.main()
