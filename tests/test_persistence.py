import math
import unittest

from src.persistence import (
    FrequencyType,
    characteristic_frequency,
    cluster_rows,
    compton_frequency,
    gap_rows,
    landscape_metric_rows,
    load_records,
    persistence_cycles,
    persistence_log10,
)


class PersistenceTests(unittest.TestCase):
    def test_compton_frequency_uses_rest_energy_over_h(self) -> None:
        self.assertTrue(math.isclose(compton_frequency(1.0), 2.417989242084918e20))

    def test_persistence_index_is_dimensionless_cycle_count(self) -> None:
        self.assertEqual(persistence_cycles(12.5, 4.0), 50.0)
        self.assertEqual(persistence_log10(1000.0), 3.0)

    def test_invalid_values_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            compton_frequency(0.0)
        with self.assertRaises(ValueError):
            persistence_cycles(1.0, -1.0)

    def test_frequency_framework_preserves_frequency_type(self) -> None:
        frequency = characteristic_frequency(7.83, FrequencyType.SCHUMANN, "measured resonance")
        self.assertEqual(frequency.frequency_type, FrequencyType.SCHUMANN)

    def test_focused_dataset_has_requested_particles(self) -> None:
        records = load_records()
        self.assertEqual(
            {record.particle for record in records},
            {"electron", "muon", "tau", "neutron", "proton", "W", "Z", "Higgs"},
        )

    def test_gap_ranking_and_clusters_are_generated(self) -> None:
        records = load_records()
        gaps = gap_rows(records)
        clusters, threshold = cluster_rows(records)
        self.assertEqual(len(gaps), len(records) - 1)
        self.assertEqual(len(clusters), len(records))
        self.assertEqual(gaps[0]["Largest_gap_rank"], 1)
        self.assertGreater(threshold, 0.0)

    def test_landscape_metrics_compare_mass_lifetime_and_persistence(self) -> None:
        metrics = landscape_metric_rows(load_records())
        self.assertEqual(
            {row["Landscape"] for row in metrics},
            {"mass_mev", "lifetime_s", "persistence_N"},
        )


if __name__ == "__main__":
    unittest.main()
