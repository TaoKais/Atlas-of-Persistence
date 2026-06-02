import math
import unittest

import numpy as np

from src.dynamic_phase_structure import (
    BASE_MAX,
    BASE_MIN,
    BASE_SAMPLES,
    circular_clusters,
    circular_gaps,
    load_entities,
    phase_matrix,
    scan_metrics,
)


class DynamicPhaseStructureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.entities = load_entities()
        cls.bases = np.linspace(BASE_MIN, BASE_MAX, BASE_SAMPLES)
        cls.phases = phase_matrix(cls.entities, cls.bases)

    def test_documented_stable_bounds_are_included(self) -> None:
        stable = self.entities.set_index("name")
        self.assertEqual(stable.loc["electron", "lifetime_handling"], "documented_lower_bound")
        self.assertEqual(stable.loc["proton", "lifetime_handling"], "documented_lower_bound")

    def test_phase_matrix_has_requested_scan_shape_and_range(self) -> None:
        self.assertEqual(self.phases.shape, (BASE_SAMPLES, len(self.entities)))
        self.assertTrue(np.all(self.phases >= 0.0))
        self.assertTrue(np.all(self.phases < 2.0 * math.pi))

    def test_circular_gaps_sum_to_full_circle(self) -> None:
        _, _, gaps = circular_gaps(self.phases[0])
        self.assertTrue(math.isclose(float(gaps.sum()), 2.0 * math.pi, rel_tol=1e-12))

    def test_cluster_split_handles_wraparound(self) -> None:
        theta = np.radians([355.0, 2.0, 8.0, 180.0])
        clusters = circular_clusters(theta, threshold_deg=20.0)
        self.assertEqual(sorted(map(len, clusters)), [1, 3])

    def test_scan_contains_requested_metrics(self) -> None:
        scan, geometry = scan_metrics(self.entities, self.bases[:3], self.phases[:3])
        self.assertEqual(len(scan), 3)
        self.assertEqual(len(geometry), 12)
        self.assertIn("largest_gap_deg", scan)
        self.assertIn("top4_symmetry_score", scan)


if __name__ == "__main__":
    unittest.main()
