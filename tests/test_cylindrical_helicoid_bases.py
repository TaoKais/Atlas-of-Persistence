import unittest

import numpy as np

from src.cylindrical_helicoid_bases import (
    BASES,
    coordinate_table,
    dbscan,
    distance_matrix,
    metric_table,
    persistent_neighbors,
    silhouette_score,
)
from src.stable_particle_sensitivity import load_variant


class CylindricalHelicoidBasesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.entities = load_variant()
        cls.coordinates = coordinate_table(cls.entities)

    def test_coordinate_table_covers_all_bases_and_finite_entities(self) -> None:
        self.assertEqual(len(self.coordinates), len(BASES) * len(self.entities))
        self.assertEqual(set(self.coordinates["base_name"]), set(BASES))
        self.assertTrue(np.allclose(self.coordinates["x"] ** 2 + self.coordinates["y"] ** 2, 1.0))

    def test_distance_matrix_is_symmetric(self) -> None:
        distances = distance_matrix(np.array([[0.0, 0.0], [3.0, 4.0]]))
        np.testing.assert_allclose(distances, [[0.0, 5.0], [5.0, 0.0]])

    def test_dbscan_and_silhouette_find_two_simple_clusters(self) -> None:
        points = np.array([[0.0, 0.0], [0.1, 0.0], [5.0, 5.0], [5.1, 5.0]])
        labels = dbscan(points, eps=0.2, min_samples=2)
        self.assertEqual(set(labels), {0, 1})
        self.assertGreater(silhouette_score(points, labels), 0.9)

    def test_metrics_and_persistent_neighbors_cover_requested_fields(self) -> None:
        metrics = metric_table(self.coordinates)
        neighbors = persistent_neighbors(self.coordinates)
        self.assertEqual(len(metrics), len(BASES))
        self.assertIn("largest_angular_gap_deg", metrics)
        self.assertIn("fraction_of_bases_close", neighbors)
        self.assertTrue(neighbors["fraction_of_bases_close"].between(0.0, 1.0).all())


if __name__ == "__main__":
    unittest.main()
