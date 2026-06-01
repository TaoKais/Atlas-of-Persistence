"""Exploratory polygon geometry for the largest complex-phase gaps.

The connected gap centers are descriptive constructions. They do not establish
a preferred logarithmic base, a physical mechanism, or a new physical law.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.complex_phase import BASES, DATA, TWO_PI

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "gap_geometry"
TOP_K_VALUES = (3, 4, 5, 6)
REQUIRED_GAP_COLUMNS = {
    "base",
    "from_entity",
    "to_entity",
    "from_theta_deg",
    "to_theta_deg",
    "gap_deg",
    "gap_fraction",
}
REQUIRED_PHASE_COLUMNS = {"base", "name", "theta_deg", "complex_x", "complex_y"}
SUMMARY_COLUMNS = [
    "base",
    "top_k",
    "selected_gaps",
    "midpoint_angles_deg",
    "polygon_area",
    "centroid_radius",
    "mean_angle_gap",
    "angle_gap_std",
    "symmetry_score",
    "qualitative_shape",
]


def load_inputs(
    gaps_path: Path = DATA / "complex_phase_gaps.csv",
    phases_path: Path = DATA / "complex_phase_table.csv",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load and validate the generated complex-phase input tables."""

    gaps = pd.read_csv(gaps_path, dtype={"base": str})
    phases = pd.read_csv(phases_path, dtype={"base": str})
    missing_gaps = REQUIRED_GAP_COLUMNS - set(gaps.columns)
    missing_phases = REQUIRED_PHASE_COLUMNS - set(phases.columns)
    if missing_gaps:
        raise ValueError(f"gap table is missing required columns: {sorted(missing_gaps)}")
    if missing_phases:
        raise ValueError(f"phase table is missing required columns: {sorted(missing_phases)}")
    return gaps, phases


def midpoint_angle_deg(from_theta_deg: float, gap_deg: float) -> float:
    """Return a circular gap midpoint in the half-open interval [0, 360)."""

    return float((from_theta_deg + gap_deg / 2.0) % 360.0)


def select_largest_gaps(gaps: pd.DataFrame, base: str, top_k: int) -> pd.DataFrame:
    """Select and annotate the largest available gaps for one base."""

    selected = (
        gaps[gaps["base"].astype(str) == str(base)]
        .sort_values("gap_deg", ascending=False)
        .head(top_k)
        .copy()
    )
    selected["mid_theta_deg"] = [
        midpoint_angle_deg(from_theta, gap)
        for from_theta, gap in zip(selected["from_theta_deg"], selected["gap_deg"])
    ]
    radians = np.radians(selected["mid_theta_deg"].to_numpy(dtype=float))
    selected["mid_x"] = np.cos(radians)
    selected["mid_y"] = np.sin(radians)
    return selected.sort_values("mid_theta_deg").reset_index(drop=True)


def polygon_metrics(midpoint_angles_deg: np.ndarray | list[float]) -> dict[str, object]:
    """Compute polygon and circular-spacing metrics for gap-center angles."""

    angles = np.sort(np.mod(np.asarray(midpoint_angles_deg, dtype=float), 360.0))
    number_of_vertices = len(angles)
    if number_of_vertices == 0:
        raise ValueError("at least one midpoint angle is required")
    radians = np.radians(angles)
    x = np.cos(radians)
    y = np.sin(radians)
    angle_gaps = np.diff(np.r_[angles, angles[0] + 360.0])
    next_x = np.roll(x, -1)
    next_y = np.roll(y, -1)
    side_lengths = np.hypot(next_x - x, next_y - y)
    polygon_area = 0.5 * abs(float(np.sum(x * next_y - y * next_x)))
    centroid_x = float(np.mean(x))
    centroid_y = float(np.mean(y))
    expected_gap = 360.0 / number_of_vertices
    symmetry_error = float(np.mean(np.abs(angle_gaps - expected_gap)))
    return {
        "number_of_vertices": number_of_vertices,
        "side_lengths": side_lengths.tolist(),
        "mean_side_length": float(np.mean(side_lengths)),
        "side_length_std": float(np.std(side_lengths)),
        "angle_gaps_between_vertices": angle_gaps.tolist(),
        "mean_angle_gap": float(np.mean(angle_gaps)),
        "angle_gap_std": float(np.std(angle_gaps)),
        "polygon_area": polygon_area,
        "centroid_x": centroid_x,
        "centroid_y": centroid_y,
        "centroid_radius": float(math.hypot(centroid_x, centroid_y)),
        "circular_symmetry_score": 1.0 / (1.0 + symmetry_error),
    }


def qualitative_shape(top_k: int, metrics: dict[str, object]) -> str:
    """Apply deliberately simple labels to a gap-center polygon."""

    symmetry_score = float(metrics["circular_symmetry_score"])
    centroid_radius = float(metrics["centroid_radius"])
    angle_gap_std = float(metrics["angle_gap_std"])
    angle_gaps = np.asarray(metrics["angle_gaps_between_vertices"], dtype=float)
    expected_gap = 360.0 / top_k
    if symmetry_score >= 0.08:
        if top_k == 3:
            return "triangle-like"
        if top_k == 4:
            return "quadrilateral / cross-like"
        if top_k == 5:
            return "pentagon-like"
    if centroid_radius < 0.15 and angle_gap_std < 20.0:
        return "balanced radial structure"
    if float(angle_gaps.max()) > 1.75 * expected_gap:
        return "asymmetric sectoring"
    return "irregular"


def geometry_records(gaps: pd.DataFrame) -> tuple[pd.DataFrame, dict[tuple[str, int], pd.DataFrame]]:
    """Compute summary records and selected gap details for every base and K."""

    rows: list[dict[str, object]] = []
    selections: dict[tuple[str, int], pd.DataFrame] = {}
    for base in BASES:
        for top_k in TOP_K_VALUES:
            selected = select_largest_gaps(gaps, base, top_k)
            metrics = polygon_metrics(selected["mid_theta_deg"].to_numpy())
            selections[(base, top_k)] = selected
            labels = (selected["from_entity"] + " -> " + selected["to_entity"]).tolist()
            rows.append(
                {
                    "base": base,
                    "top_k": top_k,
                    "selected_gaps": json.dumps(labels),
                    "midpoint_angles_deg": json.dumps(
                        [round(value, 12) for value in selected["mid_theta_deg"]]
                    ),
                    "polygon_area": metrics["polygon_area"],
                    "centroid_radius": metrics["centroid_radius"],
                    "mean_angle_gap": metrics["mean_angle_gap"],
                    "angle_gap_std": metrics["angle_gap_std"],
                    "symmetry_score": metrics["circular_symmetry_score"],
                    "qualitative_shape": qualitative_shape(top_k, metrics),
                }
            )
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS), selections


def _save_polygon_plot(phases: pd.DataFrame, selected: pd.DataFrame, base: str, top_k: int) -> None:
    frame = phases[phases["base"].astype(str) == str(base)]
    fig, ax = plt.subplots(figsize=(8, 8))
    circle = np.linspace(0.0, TWO_PI, 400)
    ax.plot(np.cos(circle), np.sin(circle), color="0.65", linewidth=1.0)
    ax.scatter(frame["complex_x"], frame["complex_y"], color="0.55", alpha=0.28, s=22)
    if not selected.empty:
        closed = pd.concat([selected, selected.iloc[[0]]], ignore_index=True)
        ax.plot(closed["mid_x"], closed["mid_y"], color="#c83e4d", linewidth=1.8)
        ax.scatter(selected["mid_x"], selected["mid_y"], color="#c83e4d", s=52, zorder=3)
        for _, row in selected.iterrows():
            label = f"{row['from_entity']} -> {row['to_entity']}"
            ax.annotate(label, (row["mid_x"], row["mid_y"]), xytext=(5, 5),
                        textcoords="offset points", fontsize=7)
    ax.set(
        aspect="equal",
        xlim=(-1.2, 1.2),
        ylim=(-1.2, 1.2),
        xlabel="cos(theta)",
        ylabel="sin(theta)",
        title=f"Gap-center polygon: base {base}, top {top_k}",
    )
    fig.tight_layout()
    fig.savefig(OUTPUT / f"gap_center_polygon_base_{base}_top_{top_k}.png", dpi=180)
    plt.close(fig)


def _save_grouped_bar(summary: pd.DataFrame, column: str, ylabel: str, filename: str) -> None:
    pivot = summary.pivot(index="base", columns="top_k", values=column).reindex(BASES)
    ax = pivot.plot(kind="bar", figsize=(9, 5), width=0.8)
    ax.set(xlabel="logarithmic base", ylabel=ylabel, title=f"{ylabel} by base and top K")
    ax.legend(title="top K")
    ax.figure.tight_layout()
    ax.figure.savefig(OUTPUT / filename, dpi=180)
    plt.close(ax.figure)


def _save_best_symmetry_radar(summary: pd.DataFrame) -> None:
    best = summary.groupby("base")["symmetry_score"].max().reindex(BASES)
    labels = best.index.tolist()
    values = best.to_numpy(dtype=float)
    angles = np.linspace(0.0, TWO_PI, len(labels), endpoint=False)
    closed_angles = np.r_[angles, angles[0]]
    closed_values = np.r_[values, values[0]]
    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={"projection": "polar"})
    ax.plot(closed_angles, closed_values, color="#2457a6", linewidth=2.0)
    ax.fill(closed_angles, closed_values, color="#2457a6", alpha=0.18)
    ax.set_xticks(angles, labels)
    ax.set_title("Best gap-center symmetry score by base")
    fig.tight_layout()
    fig.savefig(OUTPUT / "best_symmetry_score_radar.png", dpi=180)
    plt.close(fig)


def generate_plots(
    phases: pd.DataFrame,
    summary: pd.DataFrame,
    selections: dict[tuple[str, int], pd.DataFrame],
) -> None:
    """Write all requested gap-center polygon and comparison plots."""

    OUTPUT.mkdir(parents=True, exist_ok=True)
    for (base, top_k), selected in selections.items():
        _save_polygon_plot(phases, selected, base, top_k)
    _save_grouped_bar(summary, "symmetry_score", "symmetry score", "symmetry_score_by_base_and_k.png")
    _save_grouped_bar(summary, "polygon_area", "polygon area", "polygon_area_by_base_and_k.png")
    _save_best_symmetry_radar(summary)


def generate() -> pd.DataFrame:
    """Generate the gap-center geometry summary and visualizations."""

    gaps, phases = load_inputs()
    summary, selections = geometry_records(gaps)
    DATA.mkdir(parents=True, exist_ok=True)
    summary.to_csv(DATA / "gap_geometry_summary.csv", index=False)
    generate_plots(phases, summary, selections)
    return summary


if __name__ == "__main__":
    generated_summary = generate()
    print(f"Generated {len(generated_summary)} gap-center geometry rows in {OUTPUT}")
    print(generated_summary.to_string(index=False))
