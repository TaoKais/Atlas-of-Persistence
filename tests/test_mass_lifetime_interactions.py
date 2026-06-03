import math
import unittest

import numpy as np

from src.mass_lifetime_interactions import (
    H,
    MEV_J,
    expectation_comparison,
    finite_unstable_particles,
    global_correlation,
    interaction_regression,
    interaction_significance,
    load_tolerances,
    report,
    safe_log10,
    tautology_decomposition,
)


class MassLifetimeInteractionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tolerances = load_tolerances()
        cls.frame = finite_unstable_particles()

    def test_stable_particles_excluded_by_default(self) -> None:
        self.assertNotIn("stable", set(self.frame["stability"]))
        self.assertNotIn("electron", set(self.frame["name"]))
        self.assertNotIn("proton", set(self.frame["name"]))

    def test_log_transforms_reject_non_positive_values(self) -> None:
        with self.assertRaises(ValueError):
            safe_log10(0.0)
        with self.assertRaises(ValueError):
            safe_log10(-1.0)

    def test_gamma_ev_is_positive_for_finite_lifetime(self) -> None:
        self.assertTrue((self.frame["Gamma_eV"] > 0.0).all())
        self.assertTrue(np.isfinite(self.frame["Gamma_eV"]).all())

    def test_log_n_decomposition_identity_holds_within_tolerance(self) -> None:
        constant = math.log10(MEV_J / H)
        reconstructed = self.frame["log_mass"] + self.frame["log_lifetime"] + constant
        np.testing.assert_allclose(
            reconstructed,
            self.frame["logN"],
            rtol=self.tolerances["relative_numeric_tolerance"],
            atol=self.tolerances["absolute_numeric_tolerance"],
        )

    def test_regression_outputs_contain_n_and_confidence_intervals(self) -> None:
        regression = interaction_regression(self.frame, self.tolerances)
        self.assertIn("n", regression.columns)
        self.assertIn("alpha_ci_low", regression.columns)
        self.assertIn("alpha_ci_high", regression.columns)
        fitted = regression[regression["status"] == "fit"]
        self.assertTrue(fitted["alpha_ci_low"].notna().all())
        self.assertTrue(fitted["alpha_ci_high"].notna().all())

    def test_permutation_p_values_are_in_unit_interval(self) -> None:
        significance = interaction_significance(self.frame, self.tolerances)
        value = float(significance.loc[0, "permutation_p_value"])
        self.assertGreaterEqual(value, 0.0)
        self.assertLessEqual(value, 1.0)

    def test_reported_conclusions_include_tolerance_and_n(self) -> None:
        global_corr = global_correlation(self.frame, self.tolerances)
        regression = interaction_regression(self.frame, self.tolerances)
        decomposition = tautology_decomposition(self.frame)
        significance = interaction_significance(self.frame, self.tolerances)
        text = report(self.frame, global_corr, regression, decomposition, significance, self.tolerances)
        self.assertIn("n=", text)
        self.assertIn("tolerance", text)
        self.assertIn(str(self.tolerances["p_value_threshold"]), text)

    def test_weak_expected_exponent_is_approximate_not_universal(self) -> None:
        regression = interaction_regression(self.frame, self.tolerances)
        comparison = expectation_comparison(regression, self.tolerances)
        weak_note = comparison.loc[comparison["interaction"] == "weak", "caution_note"].iloc[0]
        self.assertIn("approximate", weak_note)
        self.assertIn("not universal", weak_note)
        non_weak = comparison[comparison["interaction"] != "weak"]
        self.assertTrue(non_weak["expected_alpha"].isna().all())


if __name__ == "__main__":
    unittest.main()
