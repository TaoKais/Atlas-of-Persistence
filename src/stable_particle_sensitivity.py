"""Sensitivity analysis for stable-particle lifetime handling.

Stable-particle truncations are controlled numerical assumptions, not measured
lifetimes. The resulting comparisons are exploratory and do not claim physics.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.dynamic_phase_structure import DATA, PHI, ROOT, TWO_PI, circular_gaps
from src.gap_geometry import polygon_metrics
from src.persistence import compton_frequency

OUTPUT = ROOT / "outputs" / "stable_sensitivity"
BASES = {"2": 2.0, "e": math.e, "pi": math.pi, "phi": PHI, "10": 10.0}
STABLE_TAU_VALUES = (1e25, 1e30, 1e35, 1e40)
DEFAULT_STABLE_TAU = 1e35
CONTINUOUS_BASES = np.linspace(1.2, 50.0, 1000)


def _variant_name(stable_tau: float | None) -> str:
    return "finite_lifetime_only" if stable_tau is None else f"stable_truncated_{stable_tau:.0e}"


def load_variant(stable_tau: float | None = None, path: Path = DATA / "particles.csv") -> pd.DataFrame:
    """Load finite particles or include stable particles with a declared truncation."""

    entities = pd.read_csv(path).copy()
    entities["mass_mev"] = pd.to_numeric(entities["mass_mev"], errors="coerce")
    entities["lifetime_s"] = pd.to_numeric(entities["lifetime_s"], errors="coerce")
    stable = entities["stability"].eq("stable")
    if stable_tau is None:
        entities = entities[~stable & entities["lifetime_s"].gt(0.0)].copy()
        entities["lifetime_handling"] = "finite_value"
    else:
        entities.loc[stable, "lifetime_s"] = stable_tau
        entities = entities[entities["lifetime_s"].gt(0.0)].copy()
        entities["lifetime_handling"] = np.where(
            entities["stability"].eq("stable"), "stable_truncation", "finite_value"
        )
    entities["variant"] = _variant_name(stable_tau)
    entities["tau_stable_s"] = stable_tau
    entities["compton_frequency_hz"] = entities["mass_mev"].map(compton_frequency)
    entities["compton_cycles"] = entities["compton_frequency_hz"] * entities["lifetime_s"]
    entities["log10_N"] = np.log10(entities["compton_cycles"])
    entities["ln_N"] = np.log(entities["compton_cycles"])
    return entities.reset_index(drop=True)


def phases_for_base(entities: pd.DataFrame, base: float) -> np.ndarray:
    """Return circular persistence phases for one logarithmic base."""

    return TWO_PI * np.mod(entities["ln_N"].to_numpy(dtype=float) / math.log(base), 1.0)


def _circular_centroid(theta: np.ndarray) -> tuple[float, float]:
    mean = np.mean(np.exp(1j * theta))
    return float(np.degrees(np.mod(np.angle(mean), TWO_PI))), float(1.0 - abs(mean))


def _gap_metrics(theta: np.ndarray) -> dict[str, object]:
    ordered, _, gaps = circular_gaps(theta)
    gaps_deg = np.degrees(gaps)
    row: dict[str, object] = {
        "largest_gap_deg": float(gaps_deg.max()),
        "gap_distribution_deg": json.dumps(np.sort(gaps_deg)[::-1].tolist()),
    }
    for top_k in (3, 4, 5, 6):
        selected = np.argsort(gaps)[::-1][:top_k]
        midpoint = np.mod(ordered[selected] + gaps[selected] / 2.0, TWO_PI)
        row[f"top_{top_k}_symmetry_score"] = polygon_metrics(
            np.degrees(midpoint)
        )["circular_symmetry_score"]
    return row


def sensitivity_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute summary, centroid, and particle phase tables for all variants."""

    summary_rows = []
    centroid_rows = []
    phase_rows = []
    for stable_tau in (None, *STABLE_TAU_VALUES):
        entities = load_variant(stable_tau)
        for base_label, base in BASES.items():
            theta = phases_for_base(entities, base)
            summary_rows.append(
                {
                    "variant": _variant_name(stable_tau),
                    "tau_stable_s": stable_tau,
                    "base": base_label,
                    "base_value": base,
                    "particle_count": len(entities),
                    **_gap_metrics(theta),
                }
            )
            for category_type in ("family", "dominant_interaction"):
                for category, group in entities.assign(theta_rad=theta).groupby(category_type):
                    centroid, variance = _circular_centroid(group["theta_rad"].to_numpy())
                    centroid_rows.append(
                        {
                            "variant": _variant_name(stable_tau),
                            "tau_stable_s": stable_tau,
                            "base": base_label,
                            "base_value": base,
                            "category_type": category_type,
                            "category": category,
                            "particle_count": len(group),
                            "angular_centroid_deg": centroid,
                            "angular_variance": variance,
                        }
                    )
            for (_, particle), angle in zip(entities.iterrows(), theta):
                phase_rows.append(
                    {
                        "variant": _variant_name(stable_tau),
                        "tau_stable_s": stable_tau,
                        "base": base_label,
                        "base_value": base,
                        "name": particle["name"],
                        "family": particle["family"],
                        "dominant_interaction": particle["dominant_interaction"],
                        "lifetime_handling": particle["lifetime_handling"],
                        "lifetime_s": particle["lifetime_s"],
                        "compton_frequency_hz": particle["compton_frequency_hz"],
                        "compton_cycles": particle["compton_cycles"],
                        "log10_N": particle["log10_N"],
                        "theta_rad": angle,
                        "theta_deg": math.degrees(angle),
                        "x": math.cos(angle),
                        "y": math.sin(angle),
                    }
                )
    return pd.DataFrame(summary_rows), pd.DataFrame(centroid_rows), pd.DataFrame(phase_rows)


def continuous_comparison() -> pd.DataFrame:
    """Return finite-only and default-truncation top-4 symmetry curves."""

    rows = []
    for stable_tau in (None, DEFAULT_STABLE_TAU):
        entities = load_variant(stable_tau)
        for base in CONTINUOUS_BASES:
            theta = phases_for_base(entities, float(base))
            rows.append(
                {
                    "variant": _variant_name(stable_tau),
                    "base": base,
                    "top_4_symmetry_score": _gap_metrics(theta)["top_4_symmetry_score"],
                }
            )
    return pd.DataFrame(rows)


def _family_frame(centroids: pd.DataFrame, variant: str) -> pd.DataFrame:
    return centroids[
        (centroids["variant"] == variant) & (centroids["category_type"] == "family")
    ].copy()


def _save_family_centroids(centroids: pd.DataFrame, variant: str, filename: str) -> None:
    frame = _family_frame(centroids, variant)
    fig, ax = plt.subplots(figsize=(9, 5))
    for category, group in frame.groupby("category"):
        ax.plot(group["base"], group["angular_centroid_deg"], marker="o", label=category)
    ax.set(xlabel="logarithmic base", ylabel="angular centroid (degrees)",
           title=f"Family centroid evolution: {variant}", ylim=(0, 360))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT / filename, dpi=170)
    plt.close(fig)


def _circular_difference(left: pd.Series, right: pd.Series) -> np.ndarray:
    return (left.to_numpy() - right.to_numpy() + 180.0) % 360.0 - 180.0


def _save_centroid_difference(centroids: pd.DataFrame) -> None:
    finite = _family_frame(centroids, "finite_lifetime_only")
    stable = _family_frame(centroids, _variant_name(DEFAULT_STABLE_TAU))
    merged = stable.merge(finite, on=["base", "category"], suffixes=("_b", "_a"))
    merged["difference_deg"] = _circular_difference(
        merged["angular_centroid_deg_b"], merged["angular_centroid_deg_a"]
    )
    fig, ax = plt.subplots(figsize=(9, 5))
    for category, group in merged.groupby("category"):
        ax.plot(group["base"], group["difference_deg"], marker="o", label=category)
    ax.axhline(0.0, color="0.5", linewidth=0.8)
    ax.set(xlabel="logarithmic base", ylabel="centroid B - centroid A (degrees)",
           title="Family centroid sensitivity: stable truncation minus finite-only")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT / "family_centroid_difference.png", dpi=170)
    plt.close(fig)


def _save_unit_circles(phases: pd.DataFrame) -> None:
    variants = ("finite_lifetime_only", _variant_name(DEFAULT_STABLE_TAU))
    fig, axes = plt.subplots(2, len(BASES), figsize=(18, 7))
    circle = np.linspace(0.0, TWO_PI, 400)
    for row_index, variant in enumerate(variants):
        for ax, base in zip(axes[row_index], BASES):
            frame = phases[(phases["variant"] == variant) & (phases["base"] == base)]
            ax.plot(np.cos(circle), np.sin(circle), color="0.75")
            for family, group in frame.groupby("family"):
                ax.scatter(group["x"], group["y"], s=18, label=family)
            ax.set(aspect="equal", title=f"{variant}\nbase {base}")
            ax.axis("off")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, fontsize=8)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(OUTPUT / "unit_circles_variants.png", dpi=170)
    plt.close(fig)


def _save_base_scan(comparison: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    for variant, group in comparison.groupby("variant"):
        ax.plot(group["base"], group["top_4_symmetry_score"], label=variant, linewidth=0.9)
    ax.set(xlabel="logarithmic base", ylabel="top-4 symmetry score",
           title="Base scan sensitivity to stable-particle handling")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT / "base_scan_comparison.png", dpi=170)
    plt.close(fig)


def _save_lepton_trajectories() -> None:
    entities = load_variant(DEFAULT_STABLE_TAU)
    leptons = entities[entities["family"] == "lepton"].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(10, 5))
    for _, particle in leptons.iterrows():
        theta = TWO_PI * np.mod(particle["ln_N"] / np.log(CONTINUOUS_BASES), 1.0)
        ax.plot(CONTINUOUS_BASES, np.degrees(theta), label=particle["name"], linewidth=0.8)
    ax.set(xlabel="logarithmic base", ylabel="theta (degrees)",
           title="Lepton-only phase trajectories: stable truncation 1e35 s", ylim=(0, 360))
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT / "lepton_phase_trajectories.png", dpi=170)
    plt.close(fig)


def generate_plots(
    centroids: pd.DataFrame, phases: pd.DataFrame, comparison: pd.DataFrame
) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    _save_family_centroids(centroids, "finite_lifetime_only", "family_centroids_finite_only.png")
    _save_family_centroids(
        centroids, _variant_name(DEFAULT_STABLE_TAU), "family_centroids_stable_truncated.png"
    )
    _save_centroid_difference(centroids)
    _save_unit_circles(phases)
    _save_base_scan(comparison)
    _save_lepton_trajectories()


def _category_shift(centroids: pd.DataFrame, category: str) -> float:
    finite = _family_frame(centroids, "finite_lifetime_only")
    stable = _family_frame(centroids, _variant_name(DEFAULT_STABLE_TAU))
    merged = stable.merge(finite, on=["base", "category"], suffixes=("_b", "_a"))
    selected = merged[merged["category"] == category]
    return float(np.mean(np.abs(_circular_difference(
        selected["angular_centroid_deg_b"], selected["angular_centroid_deg_a"]
    ))))


def _report(summary: pd.DataFrame, centroids: pd.DataFrame, comparison: pd.DataFrame) -> str:
    finite_curve = comparison[comparison["variant"] == "finite_lifetime_only"][
        "top_4_symmetry_score"
    ].to_numpy()
    stable_curve = comparison[comparison["variant"] == _variant_name(DEFAULT_STABLE_TAU)][
        "top_4_symmetry_score"
    ].to_numpy()
    correlation = float(np.corrcoef(finite_curve, stable_curve)[0, 1])
    mean_difference = float(np.mean(np.abs(finite_curve - stable_curve)))
    lepton_shift = _category_shift(centroids, "lepton")
    baryon_shift = _category_shift(centroids, "baryon")
    lines = [
        "# Stable Particle Sensitivity",
        "",
        "## Scope",
        "",
        "This sensitivity test compares finite-lifetime particles with controlled stable-particle "
        "truncations. A truncation is a numerical assumption, not a measured lifetime. No physical "
        "significance is claimed.",
        "",
        "## Variants",
        "",
        "- `finite_lifetime_only`: excludes rows marked stable and excludes missing lifetimes.",
        "- `stable_truncated`: includes stable particles with `tau_stable` values of `1e25`, `1e30`, "
        "`1e35`, and `1e40` seconds. Plots use `1e35 s` unless stated otherwise.",
        "",
        "## Answers",
        "",
        "- **Was electron included in previous centroid plots?** Yes. The preceding dynamic base scan "
        "included electron and proton using documented lifetime lower bounds. This test replaces that "
        "handling with explicit finite-only and controlled-truncation variants.",
        f"- **Does including electron dominate the lepton centroid?** It materially changes the lepton "
        f"centroid. The mean absolute circular shift across the five bases is `{lepton_shift:.3f}` degrees. "
        "With only muon and tau in the finite-only variant, adding electron can substantially redirect "
        "the three-particle centroid.",
        f"- **Does proton dominate the baryon centroid?** It affects the baryon centroid, with a mean "
        f"absolute circular shift of `{baryon_shift:.3f}` degrees across the five bases. The term "
        "`dominate` is not consistently justified because the direction and magnitude vary by base.",
        f"- **Are previous structures robust without stable particles?** They are only partially robust. "
        f"The continuous top-4 symmetry curves have Pearson correlation `{correlation:.3f}` and mean "
        f"absolute score difference `{mean_difference:.6f}`. Stable handling changes some peak details.",
        "- **Which conclusions survive both variants?** The phase layout remains strongly base-dependent; "
        "selected-base geometry remains exploratory; and stable-particle treatment must be stated "
        "explicitly whenever family centroids or gap symmetry are interpreted.",
        "",
        "## Interpretation Limit",
        "",
        "This is a numerical sensitivity analysis. It does not establish new physics, a preferred base, "
        "or a stable-particle lifetime model.",
        "",
        "## Reproduction",
        "",
        "```powershell",
        "python -m src.stable_particle_sensitivity",
        "```",
    ]
    return "\n".join(lines) + "\n"


def generate() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Generate sensitivity CSVs, plots, and report."""

    summary, centroids, phases = sensitivity_tables()
    comparison = continuous_comparison()
    DATA.mkdir(parents=True, exist_ok=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    summary.to_csv(DATA / "stable_sensitivity_summary.csv", index=False)
    centroids.to_csv(DATA / "stable_sensitivity_centroids.csv", index=False)
    phases.to_csv(DATA / "stable_sensitivity_phases.csv", index=False)
    comparison.to_csv(DATA / "stable_sensitivity_base_scan.csv", index=False)
    generate_plots(centroids, phases, comparison)
    (ROOT / "docs" / "stable_particle_sensitivity.md").write_text(
        _report(summary, centroids, comparison), encoding="ascii"
    )
    return summary, centroids


if __name__ == "__main__":
    generated = generate()
    print(f"Generated {len(generated[0])} summary rows and {len(generated[1])} centroid rows in {OUTPUT}")
