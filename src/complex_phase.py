"""Exploratory complex-phase landscape for particle persistence cycles.

The phase transform is a numerical representation. It does not establish a
new physical law or imply that a selected logarithmic base is physically
preferred.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.persistence import compton_frequency

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUT = ROOT / "outputs" / "complex_phase"
TWO_PI = 2.0 * math.pi
BASES = {
    "2": 2.0,
    "e": math.e,
    "pi": math.pi,
    "phi": (1.0 + math.sqrt(5.0)) / 2.0,
    "10": 10.0,
}

PHASE_COLUMNS = [
    "name",
    "family",
    "dominant_interaction",
    "mass_mev",
    "compton_frequency_hz",
    "lifetime_s",
    "compton_cycles",
    "log10_N",
    "base",
    "log_base_N",
    "phase_fraction",
    "theta_rad",
    "theta_deg",
    "complex_x",
    "complex_y",
]

GAP_COLUMNS = [
    "base",
    "from_entity",
    "to_entity",
    "from_theta_deg",
    "to_theta_deg",
    "gap_deg",
    "gap_fraction",
    "from_log10_N",
    "to_log10_N",
    "from_frequency_hz",
    "to_frequency_hz",
    "from_cycles",
    "to_cycles",
]


def load_entities(path: Path = DATA / "particles.csv") -> pd.DataFrame:
    """Load finite-lifetime entities and derive their Compton cycle counts."""

    entities = pd.read_csv(path)
    required = {"name", "family", "dominant_interaction", "mass_mev", "lifetime_s"}
    missing = required - set(entities.columns)
    if missing:
        raise ValueError(f"particle table is missing required columns: {sorted(missing)}")
    entities = entities.copy()
    entities["mass_mev"] = pd.to_numeric(entities["mass_mev"], errors="coerce")
    entities["lifetime_s"] = pd.to_numeric(entities["lifetime_s"], errors="coerce")
    entities = entities[
        entities["mass_mev"].gt(0.0) & entities["lifetime_s"].gt(0.0)
    ].copy()
    entities["compton_frequency_hz"] = entities["mass_mev"].map(compton_frequency)
    entities["compton_cycles"] = (
        entities["compton_frequency_hz"] * entities["lifetime_s"]
    )
    entities["log10_N"] = np.log10(entities["compton_cycles"])
    return entities.reset_index(drop=True)


def phase_table(entities: pd.DataFrame) -> pd.DataFrame:
    """Map positive persistence counts to unit-circle phases for each base."""

    rows: list[pd.DataFrame] = []
    for label, base in BASES.items():
        frame = entities.copy()
        frame["base"] = label
        frame["log_base_N"] = np.log(frame["compton_cycles"]) / math.log(base)
        frame["phase_fraction"] = np.mod(frame["log_base_N"], 1.0)
        frame["theta_rad"] = TWO_PI * frame["phase_fraction"]
        frame["theta_deg"] = np.degrees(frame["theta_rad"])
        frame["complex_x"] = np.cos(frame["theta_rad"])
        frame["complex_y"] = np.sin(frame["theta_rad"])
        rows.append(frame[PHASE_COLUMNS])
    return pd.concat(rows, ignore_index=True)


def gap_table(phases: pd.DataFrame) -> pd.DataFrame:
    """Compute ranked adjacent circular gaps independently for each base."""

    rows: list[dict[str, object]] = []
    for base in BASES:
        ordered = phases[phases["base"] == base].sort_values("theta_rad").reset_index(drop=True)
        for index, current in ordered.iterrows():
            following = ordered.iloc[(index + 1) % len(ordered)]
            gap_rad = float((following["theta_rad"] - current["theta_rad"]) % TWO_PI)
            rows.append(
                {
                    "base": base,
                    "from_entity": current["name"],
                    "to_entity": following["name"],
                    "from_theta_deg": current["theta_deg"],
                    "to_theta_deg": following["theta_deg"],
                    "gap_deg": math.degrees(gap_rad),
                    "gap_fraction": gap_rad / TWO_PI,
                    "from_log10_N": current["log10_N"],
                    "to_log10_N": following["log10_N"],
                    "from_frequency_hz": current["compton_frequency_hz"],
                    "to_frequency_hz": following["compton_frequency_hz"],
                    "from_cycles": current["compton_cycles"],
                    "to_cycles": following["compton_cycles"],
                }
            )
    return (
        pd.DataFrame(rows, columns=GAP_COLUMNS)
        .sort_values(["base", "gap_deg"], ascending=[True, False])
        .reset_index(drop=True)
    )


def summary_table(phases: pd.DataFrame, gaps: pd.DataFrame) -> pd.DataFrame:
    """Return descriptive circular statistics for each selected base."""

    rows = []
    for base in BASES:
        base_phases = phases[phases["base"] == base]
        base_gaps = gaps[gaps["base"] == base]["gap_deg"]
        mean_z = np.mean(np.exp(1j * base_phases["theta_rad"].to_numpy()))
        r = float(abs(mean_z))
        rows.append(
            {
                "base": base,
                "R": r,
                "circular_mean_angle_rad": float(np.mod(np.angle(mean_z), TWO_PI)),
                "circular_mean_angle_deg": float(np.degrees(np.mod(np.angle(mean_z), TWO_PI))),
                "circular_variance": 1.0 - r,
                "largest_gap_deg": float(base_gaps.max()),
                "median_gap_deg": float(base_gaps.median()),
                "mean_gap_deg": float(base_gaps.mean()),
                "number_of_entities": len(base_phases),
            }
        )
    return pd.DataFrame(rows)


def width_ratio_from_cycles(compton_cycles: pd.Series | np.ndarray) -> pd.Series | np.ndarray:
    """Return E/Gamma = 2*pi*N when tau = hbar/Gamma is applicable."""

    return TWO_PI * compton_cycles


def _category_colors(values: pd.Series) -> dict[str, object]:
    labels = sorted(values.unique())
    palette = plt.get_cmap("tab10")
    return {label: palette(index % 10) for index, label in enumerate(labels)}


def _save_unit_circle(phases: pd.DataFrame, base: str) -> None:
    frame = phases[phases["base"] == base]
    colors = _category_colors(frame["family"])
    sizes = 35.0 + 5.0 * (frame["log10_N"] - frame["log10_N"].min())
    fig, ax = plt.subplots(figsize=(9, 9))
    angle = np.linspace(0.0, TWO_PI, 400)
    ax.plot(np.cos(angle), np.sin(angle), color="0.65", linewidth=1.0)
    for family, group in frame.groupby("family"):
        ax.scatter(group["complex_x"], group["complex_y"], s=sizes.loc[group.index], label=family,
                   color=colors[family], alpha=0.8)
    for _, row in frame.iterrows():
        ax.annotate(row["name"], (row["complex_x"], row["complex_y"]), xytext=(4, 4),
                    textcoords="offset points", fontsize=7)
    ax.set(aspect="equal", xlabel="cos(theta)", ylabel="sin(theta)",
           title=f"Complex persistence phase on the unit circle: base {base}")
    ax.legend(title="family", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT / f"unit_circle_base_{base}.png", dpi=180)
    plt.close(fig)


def _save_polar(phases: pd.DataFrame, base: str) -> None:
    frame = phases[phases["base"] == base].copy()
    colors = _category_colors(frame["family"])
    span = frame["log10_N"].max() - frame["log10_N"].min()
    frame["radius"] = (frame["log10_N"] - frame["log10_N"].min()) / (span or 1.0)
    fig, ax = plt.subplots(figsize=(9, 8), subplot_kw={"projection": "polar"})
    for family, group in frame.groupby("family"):
        ax.scatter(group["theta_rad"], group["radius"], s=50, label=family,
                   color=colors[family], alpha=0.8)
    for _, row in frame.iterrows():
        ax.annotate(row["name"], (row["theta_rad"], row["radius"]), fontsize=7)
    ax.set_title(f"Normalized persistence radius vs phase: base {base}")
    ax.legend(title="family", loc="upper right", bbox_to_anchor=(1.25, 1.1), fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT / f"polar_base_{base}.png", dpi=180)
    plt.close(fig)


def _save_scatter(
    phases: pd.DataFrame, x: str, xlabel: str, filename: str, *, log10_x: bool = False
) -> None:
    colors = _category_colors(phases["dominant_interaction"])
    fig, axes = plt.subplots(3, 2, figsize=(13, 14), sharey=True)
    for ax, base in zip(axes.flat, BASES):
        frame = phases[phases["base"] == base]
        for interaction, group in frame.groupby("dominant_interaction"):
            x_values = np.log10(group[x]) if log10_x else group[x]
            ax.scatter(x_values, group["theta_deg"], s=30, label=interaction,
                       color=colors[interaction], alpha=0.8)
        ax.set(title=f"base {base}", xlabel=xlabel, ylabel="theta (degrees)", ylim=(0, 360))
    axes.flat[-1].axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, title="dominant interaction", loc="lower right")
    fig.tight_layout()
    fig.savefig(OUTPUT / filename, dpi=180)
    plt.close(fig)


def _save_gap_rankings(gaps: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(15, 16))
    for ax, base in zip(axes.flat, BASES):
        frame = gaps[gaps["base"] == base].nlargest(10, "gap_deg").sort_values("gap_deg")
        labels = frame["from_entity"] + " -> " + frame["to_entity"]
        ax.barh(labels, frame["gap_deg"], color="#2457a6")
        ax.set(title=f"Largest circular gaps: base {base}", xlabel="gap (degrees)")
    axes.flat[-1].axis("off")
    fig.tight_layout()
    fig.savefig(OUTPUT / "gap_ranking_top10.png", dpi=180)
    plt.close(fig)


def _save_heatmap(phases: pd.DataFrame) -> None:
    matrix = phases.pivot(index="name", columns="base", values="phase_fraction").reindex(columns=BASES)
    fig, ax = plt.subplots(figsize=(8, 11))
    image = ax.imshow(matrix.to_numpy(), cmap="viridis", aspect="auto", vmin=0.0, vmax=1.0)
    ax.set_xticks(range(len(matrix.columns)), matrix.columns)
    ax.set_yticks(range(len(matrix.index)), matrix.index)
    ax.set(xlabel="logarithmic base", ylabel="entity", title="Persistence-cycle phase fractions")
    fig.colorbar(image, ax=ax, label="phase fraction")
    fig.tight_layout()
    fig.savefig(OUTPUT / "phase_fraction_heatmap.png", dpi=180)
    plt.close(fig)


def generate_plots(phases: pd.DataFrame, gaps: pd.DataFrame) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for base in BASES:
        _save_unit_circle(phases, base)
        _save_polar(phases, base)
    _save_scatter(
        phases,
        "compton_frequency_hz",
        "log10(Compton frequency / Hz)",
        "frequency_vs_phase.png",
        log10_x=True,
    )
    _save_scatter(phases, "log10_N", "log10(compton cycles)", "persistence_vs_phase.png")
    _save_gap_rankings(gaps)
    _save_heatmap(phases)


def generate() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Generate the complex phase tables and plots."""

    entities = load_entities()
    phases = phase_table(entities)
    gaps = gap_table(phases)
    summary = summary_table(phases, gaps)
    DATA.mkdir(parents=True, exist_ok=True)
    phases.to_csv(DATA / "complex_phase_table.csv", index=False)
    gaps.to_csv(DATA / "complex_phase_gaps.csv", index=False)
    summary.to_csv(DATA / "complex_phase_summary.csv", index=False)
    generate_plots(phases, gaps)
    return phases, gaps, summary


if __name__ == "__main__":
    generated_phases, _, generated_summary = generate()
    print(f"Generated {len(generated_phases)} phase rows in {OUTPUT}")
    print(generated_summary.to_string(index=False))
