import unittest

import numpy as np
import pandas as pd

from src.dynamic_phase_structure import load_entities
from src.harmonic_base_geometry import (
    constant_matches,
    detect_peaks,
    fft_spectrum,
    gap_geometry_records,
    harmonic_constants,
)


class HarmonicBaseGeometryTests(unittest.TestCase):
    def test_peak_detection_ranks_prominent_maximum_first(self) -> None:
        scan = pd.DataFrame(
            {
                "base": np.arange(9.0),
                "top4_symmetry_score": [0.0, 1.0, 0.0, 0.5, 0.0, 3.0, 0.0, 0.2, 0.0],
            }
        )
        peaks, minima = detect_peaks(scan)
        self.assertEqual(peaks.iloc[0]["base"], 5.0)
        self.assertEqual(peaks.iloc[0]["rank"], 1)
        self.assertGreater(len(minima), 0)

    def test_harmonic_candidates_cover_requested_grid(self) -> None:
        constants = harmonic_constants()
        self.assertEqual(len(constants), 20 * 7)
        self.assertEqual(set(constants["expression"]), {
            "n*pi", "n*e", "n*phi", "pi^n", "e^n", "phi^n", "sqrt(n*pi)"
        })

    def test_constant_matches_rank_closest_candidate(self) -> None:
        peaks = pd.DataFrame({"rank": [1], "base": [np.pi], "score": [1.0]})
        matches = constant_matches(peaks)
        closest = matches.iloc[0]
        self.assertEqual(closest["expression"], "n*pi")
        self.assertEqual(closest["n"], 1)
        self.assertEqual(closest["absolute_error"], 0.0)

    def test_fft_identifies_sinusoid_frequency(self) -> None:
        bases = np.linspace(0.0, 10.0, 1001)
        scan = pd.DataFrame({"base": bases, "top4_symmetry_score": np.sin(2 * np.pi * 2 * bases)})
        spectrum = fft_spectrum(scan)
        dominant = spectrum.iloc[1:].nlargest(1, "power").iloc[0]
        self.assertAlmostEqual(dominant["frequency_cycles_per_base"], 2.0, places=2)

    def test_selected_geometry_contains_top_3_through_top_6(self) -> None:
        entities = load_entities()
        selected = pd.DataFrame({"rank": [1], "base": [2.0]})
        geometry, _ = gap_geometry_records(entities, selected)
        self.assertEqual(set(geometry["top_k"]), {3, 4, 5, 6})
        self.assertTrue(geometry["compactness"].ge(0.0).all())


if __name__ == "__main__":
    unittest.main()
