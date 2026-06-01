"""Compare endpoint scales for the five strongest base-2 gap-center pairs.

The comparison is descriptive. Similar endpoint ratios would motivate further
testing, but would not establish a physical mechanism or a preferred base.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from src.complex_phase import DATA
from src.gap_geometry import load_inputs, select_largest_gaps

BASE = "2"
TOP_K = 5
PAIR_OUTPUT = DATA / "gap_pair_structure.csv"
SUMMARY_OUTPUT = DATA / "gap_pair_structure_summary.csv"
ENTITY_METRICS = [
    "mass_mev",
    "compton_frequency_hz",
    "lifetime_s",
    "compton_cycles",
    "log2_N",
    "log10_N",
]
RATIO_METRICS = [
    "mass_mev",
    "compton_frequency_hz",
    "lifetime_s",
    "compton_cycles",
]
SUMMARY_COLUMNS = [
    "metric",
    "minimum_ratio",
    "maximum_ratio",
    "geometric_mean_ratio",
    "multiplicative_spread",
    "log10_ratio_std",
    "shared_scale",
]


def selected_gap_pairs(gaps: pd.DataFrame) -> pd.DataFrame:
    """Return the five gap pairs that define the strongest base-2 polygon."""

    selected = select_largest_gaps(gaps, BASE, TOP_K)
    return selected.sort_values("gap_deg", ascending=False).reset_index(drop=True)


def _entity_lookup(phases: pd.DataFrame) -> pd.DataFrame:
    """Return one base-2 row per entity with the requested physical metrics."""

    frame = phases[phases["base"].astype(str) == BASE].copy()
    frame["log2_N"] = np.log2(frame["compton_cycles"])
    indexed = frame.set_index("name")
    if not indexed.index.is_unique:
        raise ValueError("base-2 phase table contains duplicate entity names")
    return indexed


def _geometric_mean(left: float, right: float) -> float:
    return float(math.sqrt(left * right))


def pair_structure_table(gaps: pd.DataFrame, phases: pd.DataFrame) -> pd.DataFrame:
    """Build endpoint, ratio, difference, and geometric-mean metrics per pair."""

    pairs = selected_gap_pairs(gaps)
    entities = _entity_lookup(phases)
    rows: list[dict[str, object]] = []
    for _, gap in pairs.iterrows():
        from_entity = str(gap["from_entity"])
        to_entity = str(gap["to_entity"])
        left = entities.loc[from_entity]
        right = entities.loc[to_entity]
        row: dict[str, object] = {
            "base": BASE,
            "from_entity": from_entity,
            "to_entity": to_entity,
            "gap_deg": float(gap["gap_deg"]),
            "gap_fraction": float(gap["gap_fraction"]),
            "midpoint_angle_deg": float(gap["mid_theta_deg"]),
        }
        for metric in ENTITY_METRICS:
            from_value = float(left[metric])
            to_value = float(right[metric])
            row[f"from_{metric}"] = from_value
            row[f"to_{metric}"] = to_value
            row[f"difference_{metric}"] = to_value - from_value
            if metric in RATIO_METRICS:
                row[f"ratio_{metric}"] = to_value / from_value
                row[f"geometric_mean_{metric}"] = _geometric_mean(from_value, to_value)
        row["absolute_difference_log2_N"] = abs(float(row["difference_log2_N"]))
        row["absolute_difference_log10_N"] = abs(float(row["difference_log10_N"]))
        rows.append(row)
    return pd.DataFrame(rows)


def ratio_summary_table(pairs: pd.DataFrame) -> pd.DataFrame:
    """Measure whether pair ratios cluster around a common multiplicative scale.

    Ratios are orientation-independent here: each ratio is folded to be at
    least one. ``multiplicative_spread`` is max/min. A shared-scale candidate
    requires that spread to stay within a factor of two.
    """

    rows: list[dict[str, object]] = []
    for metric in RATIO_METRICS:
        ratios = pairs[f"ratio_{metric}"].to_numpy(dtype=float)
        folded = np.maximum(ratios, 1.0 / ratios)
        minimum = float(folded.min())
        maximum = float(folded.max())
        spread = maximum / minimum
        rows.append(
            {
                "metric": metric,
                "minimum_ratio": minimum,
                "maximum_ratio": maximum,
                "geometric_mean_ratio": float(np.exp(np.mean(np.log(folded)))),
                "multiplicative_spread": spread,
                "log10_ratio_std": float(np.std(np.log10(folded))),
                "shared_scale": bool(spread <= 2.0),
            }
        )
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def generate() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Write the detailed five-pair comparison and ratio summary CSV files."""

    gaps, phases = load_inputs()
    pairs = pair_structure_table(gaps, phases)
    summary = ratio_summary_table(pairs)
    DATA.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(PAIR_OUTPUT, index=False)
    summary.to_csv(SUMMARY_OUTPUT, index=False)
    return pairs, summary


if __name__ == "__main__":
    generated_pairs, generated_summary = generate()
    print("Selected base-2 top-5 gap pairs:")
    print(generated_pairs[["from_entity", "to_entity", "gap_deg"]].to_string(index=False))
    print("\nRatio-scale comparison:")
    print(generated_summary.to_string(index=False))
