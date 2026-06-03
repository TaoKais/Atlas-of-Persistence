import unittest

import numpy as np

from src.invariant_graph_validation import (
    angular_distance_deg,
    angular_graph,
    apply_control,
    base_grid,
    graph_metrics,
    load_entities,
    validation_summary,
)


class InvariantGraphValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.entities = load_entities()
        _, cls.bases = base_grid(20)

    def test_angular_distance_wraps_at_360_degrees(self) -> None:
        self.assertAlmostEqual(float(angular_distance_deg(358.0, 2.0)), 4.0)

    def test_graph_edges_respect_tolerance(self) -> None:
        strict = angular_graph(self.entities, self.bases, 5.0)
        loose = angular_graph(self.entities, self.bases, 30.0)
        self.assertLessEqual(strict.connected_counts.sum(), loose.connected_counts.sum())

    def test_persistence_fraction_in_unit_interval(self) -> None:
        build = angular_graph(self.entities, self.bases, 10.0)
        self.assertTrue(np.all(build.weights >= 0.0))
        self.assertTrue(np.all(build.weights <= 1.0))

    def test_empirical_p_values_in_unit_interval(self) -> None:
        summary, _, _ = validation_summary(
            self.entities, self.bases, self.bases[:5], "angular_closeness", "angular_close_deg", 10.0,
            0.7, 2, 0.05, np.random.default_rng(1)
        )
        self.assertTrue(np.all(summary["empirical_p_value"].between(0.0, 1.0)))

    def test_randomized_controls_preserve_intended_distributions(self) -> None:
        rng = np.random.default_rng(2)
        controlled = apply_control(self.entities, "shuffle_lifetimes", rng)
        self.assertCountEqual(controlled["lifetime_s"].round(20), self.entities["lifetime_s"].round(20))
        self.assertCountEqual(controlled["mass_mev"].round(12), self.entities["mass_mev"].round(12))

    def test_monte_carlo_summary_has_no_nans(self) -> None:
        summary, _, _ = validation_summary(
            self.entities, self.bases, self.bases[:5], "angular_closeness", "angular_close_deg", 10.0,
            0.7, 2, 0.05, np.random.default_rng(3)
        )
        numeric = summary.select_dtypes(include=[float, int])
        self.assertFalse(numeric.isna().any().any())

    def test_every_output_row_includes_tolerance_value(self) -> None:
        summary, _, _ = validation_summary(
            self.entities, self.bases, self.bases[:5], "angular_closeness", "angular_close_deg", 10.0,
            0.7, 2, 0.05, np.random.default_rng(4)
        )
        self.assertIn("tolerance_value", summary.columns)
        self.assertFalse(summary["tolerance_value"].isna().any())

    def test_every_graph_has_node_count_equal_to_valid_entities(self) -> None:
        build = angular_graph(self.entities, self.bases, 10.0)
        metrics = graph_metrics(build.adjacency, build.weights, self.entities["family"])
        self.assertEqual(build.adjacency.shape[0], len(self.entities))
        self.assertGreaterEqual(metrics["largest_component_size"], 1)


if __name__ == "__main__":
    unittest.main()
