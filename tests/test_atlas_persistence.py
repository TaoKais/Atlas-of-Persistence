import math
import unittest

from src.atlas_persistence import (
    C,
    G,
    Particle,
    compactness,
    exergy_destruction_fraction,
    leave_one_out_metrics,
    logarithmic_gaps,
    ranks,
    rho_energy,
    ridge_predict,
    spearman,
)


class AtlasPersistenceTests(unittest.TestCase):
    def test_compton_frequency_and_mass_have_identical_ratios(self) -> None:
        electron = Particle("electron", "lepton", 0.51099895069, None, "stable", 0, "stable", "")
        muon = Particle("muon", "lepton", 105.6583755, 2.1969811e-6, "unstable", 1, "weak", "")
        self.assertAlmostEqual(
            muon.compton_frequency_hz / electron.compton_frequency_hz,
            muon.mass_mev / electron.mass_mev,
        )

    def test_logarithmic_gap_preserves_mass_ratio(self) -> None:
        low = Particle("low", "test", 1.0, None, "stable", 0, "stable", "")
        high = Particle("high", "test", 10.0, 1.0, "unstable", 1, "weak", "")
        gap = logarithmic_gaps([high, low])[0]
        self.assertAlmostEqual(float(gap["frequency_ratio"]), 10.0)
        self.assertAlmostEqual(float(gap["mass_ratio"]), 10.0)
        self.assertAlmostEqual(float(gap["log10_gap"]), 1.0)

    def test_schwarzschild_radius_has_half_compactness(self) -> None:
        mass = 1.9885e30
        schwarzschild_radius = 2.0 * G * mass / C**2
        self.assertAlmostEqual(compactness(mass, schwarzschild_radius), 0.5)

    def test_exergy_destruction_fraction_uses_gouy_stodola(self) -> None:
        self.assertAlmostEqual(exergy_destruction_fraction(1000.0, 2.5, 298.15), 0.745375)

    def test_lhc_proton_rho_is_total_energy_over_rest_energy(self) -> None:
        rho = rho_energy(6.8e6, 938.27208943)
        self.assertTrue(math.isclose(rho, 7247.36, rel_tol=1e-4))

    def test_rank_correlation_handles_ties(self) -> None:
        self.assertEqual(ranks([10.0, 20.0, 20.0, 40.0]), [1.0, 2.5, 2.5, 4.0])
        self.assertAlmostEqual(spearman([1.0, 2.0, 3.0], [2.0, 4.0, 6.0]), 1.0)

    def test_ridge_predict_fits_simple_trend(self) -> None:
        prediction = ridge_predict([[0.0], [1.0], [2.0]], [1.0, 3.0, 5.0], [3.0], alpha=0.0)
        self.assertAlmostEqual(prediction, 7.0)

    def test_leave_one_out_supports_intercept_baseline(self) -> None:
        particles = [
            Particle("a", "test", 1.0, 1.0, "unstable", 1, "weak", ""),
            Particle("b", "test", 2.0, 10.0, "unstable", 1, "weak", ""),
            Particle("c", "test", 3.0, 100.0, "unstable", 1, "weak", ""),
        ]
        mae, rmse = leave_one_out_metrics(particles, [])
        self.assertAlmostEqual(mae, 1.0)
        self.assertGreater(rmse, mae)


if __name__ == "__main__":
    unittest.main()
