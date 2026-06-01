import math
import unittest

import numpy as np

from src.complex_phase import (
    BASES,
    TWO_PI,
    gap_table,
    load_entities,
    phase_table,
    summary_table,
    width_ratio_from_cycles,
)


class ComplexPhaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.entities = load_entities()
        cls.phases = phase_table(cls.entities)
        cls.gaps = gap_table(cls.phases)
        cls.summary = summary_table(cls.phases, cls.gaps)

    def test_stable_entities_without_finite_bounds_are_excluded(self) -> None:
        self.assertNotIn("electron", set(self.entities["name"]))
        self.assertNotIn("proton", set(self.entities["name"]))

    def test_phase_fraction_is_in_unit_interval(self) -> None:
        self.assertTrue(self.phases["phase_fraction"].ge(0.0).all())
        self.assertTrue(self.phases["phase_fraction"].lt(1.0).all())

    def test_theta_is_in_expected_interval(self) -> None:
        self.assertTrue(self.phases["theta_rad"].ge(0.0).all())
        self.assertTrue(self.phases["theta_rad"].lt(TWO_PI).all())

    def test_gap_fractions_sum_to_one_for_each_base(self) -> None:
        sums = self.gaps.groupby("base")["gap_fraction"].sum()
        for base in BASES:
            self.assertTrue(math.isclose(sums[base], 1.0, rel_tol=1e-12, abs_tol=1e-12))

    def test_resultant_length_is_in_unit_interval(self) -> None:
        self.assertTrue(self.summary["R"].ge(0.0).all())
        self.assertTrue(self.summary["R"].le(1.0).all())

    def test_width_relation_matches_established_lifetime_identity(self) -> None:
        cycles = np.array([1.0, 3.5, 12.0])
        np.testing.assert_allclose(width_ratio_from_cycles(cycles), TWO_PI * cycles)


if __name__ == "__main__":
    unittest.main()
