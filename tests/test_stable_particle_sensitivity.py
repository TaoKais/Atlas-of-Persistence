import unittest

from src.stable_particle_sensitivity import (
    BASES,
    DEFAULT_STABLE_TAU,
    STABLE_TAU_VALUES,
    load_variant,
    sensitivity_tables,
)


class StableParticleSensitivityTests(unittest.TestCase):
    def test_finite_variant_excludes_stable_particles(self) -> None:
        entities = load_variant()
        self.assertNotIn("electron", set(entities["name"]))
        self.assertNotIn("proton", set(entities["name"]))
        self.assertTrue(entities["lifetime_s"].notna().all())

    def test_truncated_variant_includes_stable_particles_with_declared_tau(self) -> None:
        entities = load_variant(DEFAULT_STABLE_TAU).set_index("name")
        self.assertEqual(entities.loc["electron", "lifetime_s"], DEFAULT_STABLE_TAU)
        self.assertEqual(entities.loc["proton", "lifetime_s"], DEFAULT_STABLE_TAU)
        self.assertEqual(entities.loc["electron", "lifetime_handling"], "stable_truncation")

    def test_tables_cover_all_variants_and_bases(self) -> None:
        summary, centroids, phases = sensitivity_tables()
        self.assertEqual(len(summary), (1 + len(STABLE_TAU_VALUES)) * len(BASES))
        self.assertIn("top_6_symmetry_score", summary)
        self.assertIn("dominant_interaction", set(centroids["category_type"]))
        self.assertIn("electron", set(phases["name"]))


if __name__ == "__main__":
    unittest.main()
