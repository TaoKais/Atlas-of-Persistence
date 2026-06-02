import unittest

import numpy as np

from src.exploratory_3d_geometry import (
    REFERENCE_BASE,
    geometry_coordinates,
    normalized,
    randomized_entities,
)
from src.stable_particle_sensitivity import load_variant


class Exploratory3DGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.entities = load_variant()

    def test_normalized_values_span_unit_interval(self) -> None:
        values = normalized(np.array([3.0, 5.0, 7.0]))
        np.testing.assert_allclose(values, [0.0, 0.5, 1.0])

    def test_reference_base_cylindrical_phase_is_wrapped_log10(self) -> None:
        frame = geometry_coordinates(self.entities, REFERENCE_BASE)
        expected = 2.0 * np.pi * np.mod(frame["log10_N"], 1.0)
        np.testing.assert_allclose(frame["theta_rad"], expected)

    def test_spherical_coordinates_lie_on_unit_sphere(self) -> None:
        frame = geometry_coordinates(self.entities)
        radius = np.sqrt(
            frame["spherical_x"] ** 2 + frame["spherical_y"] ** 2 + frame["spherical_z"] ** 2
        )
        np.testing.assert_allclose(radius, 1.0)

    def test_randomized_control_preserves_positive_cycles(self) -> None:
        frame = randomized_entities(self.entities, "lifetime")
        self.assertTrue(frame["compton_cycles"].gt(0.0).all())
        self.assertEqual(set(frame["lifetime_s"]), set(self.entities["lifetime_s"]))


if __name__ == "__main__":
    unittest.main()
