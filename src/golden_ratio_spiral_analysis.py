"""Falsification-oriented golden-ratio spiral robustness analysis.

The plots in this module are coordinate constructions. A visual spiral, a
high score, or a small phi error does not establish physical significance.
"""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import pandas as pd

from src.cylindrical_helicoid_bases import dbscan, distance_matrix
from src.dynamic_phase_structure import ROOT, TWO_PI, circular_gaps
from src.exploratory_3d_geometry import normalized
from src.stable_particle_sensitivity import load_variant, phases_for_base

DATA = ROOT / "data"
OUTPUT = ROOT / "outputs" / "golden_ratio_spiral"
PHI = (1.0 + math.sqrt(5.0)) / 2.0
GOLDEN_K = 2.0 * math.log(PHI) / math.pi
BASE_SAMPLES = 5000
SCAN_BASES = np.linspace(1.2, 50.0, BASE_SAMPLES)
SPECIAL_BASES = {
    "phi": PHI, "pi": math.pi, "e": math.e, "2": 2.0, "10": 10.0,
    "2*pi": TWO_PI, "sqrt(2*pi)": math.sqrt(TWO_PI), "pi^2": math.pi**2,
    "e^pi": math.e**math.pi,
}
PHI_NEIGHBORS = np.array([PHI + delta for delta in (-.20, -.10, -.05, -.01, 0, .01, .05, .10, .20)])
STABLE_MODES = {
    "finite_lifetime_only": None,
    "stable_truncated_1e25": 1e25,
    "stable_truncated_1e30": 1e30,
    "stable_truncated_1e35": 1e35,
}
FIT_ORDERINGS = ("sorted_theta", "sorted_log10_N", "persistent_neighbor_path", "nearest_neighbor_3d")
DEFAULT_FIT_ORDERING = "sorted_log10_N"
RANDOM_SEED = 20260602
CONTROL_REPLICATES = 100


def load_entities(mode: str = "finite_lifetime_only") -> pd.DataFrame:
    """Load a declared stable-particle handling mode.

    Infinity markers cannot enter logarithmic fits, so that mode keeps finite
    rows for metrics and exposes marker names in ``attrs`` for callers.
    """

    if mode == "stable_as_infinity_marker":
        finite = load_variant()
        source = pd.read_csv(DATA / "particles.csv")
        finite.attrs["infinity_markers"] = source.loc[source["stability"].eq("stable"), "name"].tolist()
        finite.attrs["stable_mode"] = mode
        return finite
    if mode not in STABLE_MODES:
        raise ValueError(f"Unknown stable-particle mode: {mode}")
    entities = load_variant(STABLE_MODES[mode])
    entities.attrs["stable_mode"] = mode
    entities.attrs["infinity_markers"] = []
    return entities


def radial_coordinates(entities: pd.DataFrame, base: float) -> pd.DataFrame:
    """Return cylindrical, radial, and conical helicoid coordinates."""

    theta = (entities["phase_override"].to_numpy(dtype=float)
             if "phase_override" in entities else phases_for_base(entities, base))
    log_n = entities["log10_N"].to_numpy(dtype=float)
    radius = normalized(log_n)
    frame = entities[["name", "family", "dominant_interaction", "mass_mev",
                      "compton_frequency_hz", "lifetime_s", "stability",
                      "compton_cycles", "log10_N"]].copy()
    frame["base"] = float(base)
    frame["theta"] = theta
    if "phase_override" in entities:
        frame["phase_override"] = theta
    frame["r"] = radius
    frame["cylindrical_x"], frame["cylindrical_y"] = np.cos(theta), np.sin(theta)
    frame["radial_x"], frame["radial_y"] = radius * np.cos(theta), radius * np.sin(theta)
    frame["conical_x"], frame["conical_y"] = log_n * np.cos(theta), log_n * np.sin(theta)
    frame["z"] = log_n
    return frame


def _greedy_path(points: np.ndarray, start: int) -> np.ndarray:
    remaining, path = set(range(len(points))), [start]
    remaining.remove(start)
    while remaining:
        current = path[-1]
        nxt = min(remaining, key=lambda index: float(np.linalg.norm(points[current] - points[index])))
        path.append(nxt)
        remaining.remove(nxt)
    return np.asarray(path, dtype=int)


def ordering_indices(frame: pd.DataFrame, ordering: str) -> np.ndarray:
    """Return a deterministic path for sparse, unordered spiral-fit points."""

    if ordering == "sorted_theta":
        return np.argsort(frame["theta"].to_numpy())
    if ordering == "sorted_log10_N":
        return np.argsort(frame["log10_N"].to_numpy())
    points = frame[["radial_x", "radial_y", "z"]].to_numpy()
    if ordering == "nearest_neighbor_3d":
        return _greedy_path(points, int(np.argmin(frame["log10_N"].to_numpy())))
    if ordering == "persistent_neighbor_path":
        # Start with the most persistent close pair over a local, fixed base window.
        entities = frame.copy()
        bases = np.linspace(max(1.200001, float(frame["base"].iloc[0]) - .05), float(frame["base"].iloc[0]) + .05, 21)
        distance_sum = np.zeros((len(frame), len(frame)))
        for base in bases:
            theta = TWO_PI * np.mod(np.log(entities["compton_cycles"].to_numpy()) / math.log(base), 1.0)
            points_at_base = np.column_stack((normalized(entities["log10_N"]) * np.cos(theta),
                                              normalized(entities["log10_N"]) * np.sin(theta),
                                              entities["z"]))
            distance_sum += distance_matrix(points_at_base)
        np.fill_diagonal(distance_sum, np.inf)
        start, second = np.unravel_index(np.argmin(distance_sum), distance_sum.shape)
        path, remaining = [int(start), int(second)], set(range(len(frame))) - {int(start), int(second)}
        while remaining:
            nxt = min(remaining, key=lambda index: distance_sum[path[-1], index])
            path.append(int(nxt))
            remaining.remove(nxt)
        return np.asarray(path)
    raise ValueError(f"Unknown ordering: {ordering}")


def spiral_fit(frame: pd.DataFrame, ordering: str = DEFAULT_FIT_ORDERING) -> dict[str, object]:
    """Fit log(r) = intercept + k * unwrapped(theta) on positive radii."""

    ordered = frame.iloc[ordering_indices(frame, ordering)].copy()
    ordered = ordered[ordered["r"].gt(0.0)].copy()
    if len(ordered) < 3:
        raise ValueError("Logarithmic spiral fit requires at least three positive radii")
    theta = np.unwrap(ordered["theta"].to_numpy(dtype=float))
    log_r = np.log(ordered["r"].to_numpy(dtype=float))
    design = np.column_stack((np.ones(len(theta)), theta))
    intercept, fitted_k = np.linalg.lstsq(design, log_r, rcond=None)[0]
    predicted = design @ np.array([intercept, fitted_k])
    residuals = log_r - predicted
    rss = float(np.sum(residuals**2))
    tss = float(np.sum((log_r - np.mean(log_r))**2))
    r_squared = 1.0 - rss / tss if tss else 0.0
    n, parameters = len(theta), 2
    variance = max(rss / n, np.finfo(float).tiny)
    return {
        "ordering": ordering, "fitted_k": float(fitted_k), "golden_k": GOLDEN_K,
        "golden_spiral_error": abs(float(fitted_k) - GOLDEN_K),
        "golden_spiral_relative_error": abs(float(fitted_k) - GOLDEN_K) / GOLDEN_K,
        "spiral_fit_R2": float(r_squared), "a": float(math.exp(intercept)),
        "aic": float(n * math.log(variance) + 2 * parameters),
        "bic": float(n * math.log(variance) + parameters * math.log(n)),
        "theta_unwrapped": theta, "log_r": log_r, "predicted_log_r": predicted,
        "residuals": residuals, "names": ordered["name"].tolist(),
    }


def _entropy(theta: np.ndarray) -> float:
    counts, _ = np.histogram(theta, bins=12, range=(0.0, TWO_PI))
    probabilities = counts[counts > 0] / counts.sum()
    return float(-np.sum(probabilities * np.log2(probabilities)))


def _category_compactness(frame: pd.DataFrame, category: str) -> float:
    values = []
    for _, group in frame.groupby(category):
        points = group[["radial_x", "radial_y", "z"]].to_numpy()
        if len(points) > 1:
            values.extend(distance_matrix(points)[np.triu_indices(len(points), 1)])
    mean = float(np.mean(values)) if values else 0.0
    return 1.0 / (1.0 + mean)


def _persistent_neighbor_stability(frame: pd.DataFrame) -> float:
    """Measure nearest-neighbor identity stability under a fixed local base perturbation."""

    log_n = frame["log10_N"].to_numpy(dtype=float)
    radius = normalized(log_n)
    nearest = []
    base = float(frame["base"].iloc[0])
    for candidate in (max(1.200001, base - .01), base, base + .01):
        theta = (frame["phase_override"].to_numpy()
                 if "phase_override" in frame else
                 TWO_PI * np.mod(np.log(frame["compton_cycles"].to_numpy()) / math.log(candidate), 1.0))
        points = np.column_stack((radius * np.cos(theta), radius * np.sin(theta), log_n))
        distances = distance_matrix(points)
        np.fill_diagonal(distances, np.inf)
        nearest.append(np.argmin(distances, axis=1))
    return float(np.mean((nearest[0] == nearest[1]) & (nearest[1] == nearest[2])))


def base_metrics(entities: pd.DataFrame, base: float, ordering: str = DEFAULT_FIT_ORDERING) -> dict[str, float]:
    """Compute geometry metrics for one base."""

    frame = radial_coordinates(entities, base)
    theta = frame["theta"].to_numpy()
    _, _, gaps = circular_gaps(theta)
    gaps = np.sort(gaps)[::-1]
    points = frame[["radial_x", "radial_y", "z"]].to_numpy()
    distances = distance_matrix(points)
    np.fill_diagonal(distances, np.inf)
    labels = dbscan(points)
    fit = spiral_fit(frame, ordering)
    return {
        "base": float(base), "angular_entropy": _entropy(theta),
        "circular_resultant_R": float(abs(np.mean(np.exp(1j * theta)))),
        "largest_angular_gap": float(gaps[0]), "top_3_gap_sum": float(gaps[:3].sum()),
        "top_5_gap_sum": float(gaps[:5].sum()), "dbscan_cluster_count": int(len(set(labels)) - (-1 in labels)),
        "average_nearest_neighbor_distance_3d": float(np.min(distances, axis=1).mean()),
        "persistent_neighbor_stability": _persistent_neighbor_stability(frame),
        "radial_spiral_fit_R2": float(fit["spiral_fit_R2"]),
        "fitted_k": float(fit["fitted_k"]), "golden_k": GOLDEN_K,
        "golden_spiral_error": float(fit["golden_spiral_error"]),
        "golden_spiral_relative_error": float(fit["golden_spiral_relative_error"]),
        "family_compactness": _category_compactness(frame, "family"),
        "interaction_compactness": _category_compactness(frame, "dominant_interaction"),
    }


def _unit(values: pd.Series) -> pd.Series:
    span = values.max() - values.min()
    return (values - values.min()) / (span if span else 1.0)


def add_coherence_score(metrics: pd.DataFrame) -> pd.DataFrame:
    """Add the declared heuristic score after scanning all bases."""

    result = metrics.copy()
    product = (_unit(result["radial_spiral_fit_R2"]) * _unit(result["family_compactness"])
               * _unit(result["interaction_compactness"]) * _unit(result["persistent_neighbor_stability"]))
    result["geometry_coherence_score"] = product - _unit(result["golden_spiral_error"]) - _unit(result["angular_entropy"])
    return result


def scan_base_metrics(entities: pd.DataFrame, bases: np.ndarray = SCAN_BASES) -> pd.DataFrame:
    """Compute metrics on the continuous scan plus exact named and neighborhood bases."""

    all_bases = np.unique(np.r_[bases, list(SPECIAL_BASES.values()), PHI_NEIGHBORS])
    table = add_coherence_score(pd.DataFrame(base_metrics(entities, float(base)) for base in all_bases))
    table["base_name"] = ""
    for name, value in SPECIAL_BASES.items():
        table.loc[np.isclose(table["base"], value, rtol=0.0, atol=1e-14), "base_name"] = name
    table["is_phi_neighbor"] = table["base"].isin(PHI_NEIGHBORS)
    return table.sort_values("base").reset_index(drop=True)


def randomized_entities(entities: pd.DataFrame, mode: str, random: np.random.Generator) -> pd.DataFrame:
    """Create a randomized control while retaining the metadata columns."""

    frame = entities.copy()
    if mode == "shuffle_lifetimes":
        frame["lifetime_s"] = random.permutation(frame["lifetime_s"].to_numpy())
        frame["compton_cycles"] = frame["compton_frequency_hz"] * frame["lifetime_s"]
    elif mode == "shuffle_masses":
        order = random.permutation(len(frame))
        frame["mass_mev"] = frame["mass_mev"].to_numpy()[order]
        frame["compton_frequency_hz"] = frame["compton_frequency_hz"].to_numpy()[order]
        frame["compton_cycles"] = frame["compton_frequency_hz"] * frame["lifetime_s"]
    elif mode == "shuffle_N_values":
        frame["compton_cycles"] = random.permutation(frame["compton_cycles"].to_numpy())
    elif mode == "random_logN_distribution":
        low, high = frame["log10_N"].min(), frame["log10_N"].max()
        frame["compton_cycles"] = 10.0 ** random.uniform(low, high, len(frame))
    elif mode == "random_angular_phase":
        frame["phase_override"] = random.uniform(0.0, TWO_PI, len(frame))
    else:
        raise ValueError(f"Unknown randomized control: {mode}")
    frame["log10_N"] = np.log10(frame["compton_cycles"])
    frame["ln_N"] = np.log(frame["compton_cycles"])
    return frame.reset_index(drop=True)


def control_metrics(entities: pd.DataFrame, base: float = PHI, replicates: int = CONTROL_REPLICATES) -> pd.DataFrame:
    """Compare observed phi metrics with randomized control distributions."""

    observed = base_metrics(entities, base)
    rows = [{"control": "observed", "replicate": -1, **observed}]
    random = np.random.default_rng(RANDOM_SEED)
    for mode in ("shuffle_lifetimes", "shuffle_masses", "shuffle_N_values",
                 "random_logN_distribution", "random_angular_phase"):
        control_rows = []
        for replicate in range(replicates):
            randomized = randomized_entities(entities, mode, random)
            control_rows.append({"control": mode, "replicate": replicate, **base_metrics(randomized, base)})
        rows.extend(control_rows)
    table = pd.DataFrame(rows)
    observed_score = float(table.loc[table["control"].eq("observed"), "geometry_coherence_score"].iloc[0]) if "geometry_coherence_score" in table else math.nan
    # Score normalizations are meaningful only across the assembled control table.
    table = add_coherence_score(table.drop(columns=["geometry_coherence_score"], errors="ignore"))
    observed_score = float(table.loc[table["control"].eq("observed"), "geometry_coherence_score"].iloc[0])
    table["observed_phi_score"] = observed_score
    table["empirical_p_value_score_ge_observed"] = np.nan
    for mode in table["control"].unique():
        if mode == "observed":
            continue
        mask = table["control"].eq(mode)
        table.loc[mask, "empirical_p_value_score_ge_observed"] = (
            1.0 + np.sum(table.loc[mask, "geometry_coherence_score"] >= observed_score)
        ) / (1.0 + mask.sum())
    return table


def stable_mode_metrics() -> pd.DataFrame:
    rows = []
    for mode in (*STABLE_MODES, "stable_as_infinity_marker"):
        metrics = base_metrics(load_entities(mode), PHI)
        rows.append({"stable_mode": mode, **metrics})
    return pd.DataFrame(rows)


def phi_ordering_fits(entities: pd.DataFrame) -> pd.DataFrame:
    rows = []
    frame = radial_coordinates(entities, PHI)
    for ordering in FIT_ORDERINGS:
        fit = spiral_fit(frame, ordering)
        rows.append({key: value for key, value in fit.items() if not isinstance(value, (np.ndarray, list))})
    return pd.DataFrame(rows)


def _colors(frame: pd.DataFrame) -> list[object]:
    labels = sorted(frame["family"].unique())
    palette = {label: plt.get_cmap("tab10")(index) for index, label in enumerate(labels)}
    return [palette[label] for label in frame["family"]]


def _save_3d(frame: pd.DataFrame, prefix: str, stem: str, title: str) -> None:
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(projection="3d")
    order = np.argsort(frame["z"])
    ax.plot(frame[f"{prefix}_x"].iloc[order], frame[f"{prefix}_y"].iloc[order], frame["z"].iloc[order], color="#2457a6", alpha=.65)
    ax.scatter(frame[f"{prefix}_x"], frame[f"{prefix}_y"], frame["z"], c=_colors(frame), s=28)
    ax.set(xlabel="x", ylabel="y", zlabel="log10(N)", title=title)
    fig.tight_layout()
    fig.savefig(OUTPUT / f"{stem}.png", dpi=170)
    plt.close(fig)


def save_visualizations(entities: pd.DataFrame, metrics: pd.DataFrame, controls: pd.DataFrame) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    frames = {name: radial_coordinates(entities, base) for name, base in SPECIAL_BASES.items() if name in ("phi", "pi", "e", "10")}
    for name, frame in frames.items():
        _save_3d(frame, "radial", f"radial_helicoid_{name}", f"Radial helicoid: base {name}")
    phi_frame = frames["phi"]
    _save_3d(phi_frame, "cylindrical", "cylindrical_helicoid_phi", "Cylindrical helicoid: base phi")
    _save_3d(phi_frame, "conical", "conical_helicoid_phi", "Conical helicoid: base phi")

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(phi_frame["radial_x"], phi_frame["radial_y"], c=_colors(phi_frame))
    ax.set(aspect="equal", xlabel="radial x", ylabel="radial y", title="Top-down radial projection: base phi")
    fig.tight_layout(); fig.savefig(OUTPUT / "top_down_polar_projection_phi.png", dpi=170); plt.close(fig)

    fit = spiral_fit(phi_frame)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(fit["theta_unwrapped"], fit["log_r"], label="particles")
    ax.plot(fit["theta_unwrapped"], fit["predicted_log_r"], label="fitted logarithmic spiral")
    ax.set(xlabel="unwrapped theta", ylabel="log(r)", title="Exploratory spiral fit: base phi")
    ax.legend(); fig.tight_layout(); fig.savefig(OUTPUT / "spiral_fit_overlay_phi.png", dpi=170); plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axhline(0, color=".5", linewidth=.8); ax.scatter(fit["theta_unwrapped"], fit["residuals"])
    ax.set(xlabel="unwrapped theta", ylabel="log(r) residual", title="Phi spiral-fit residuals")
    fig.tight_layout(); fig.savefig(OUTPUT / "spiral_fit_residuals_phi.png", dpi=170); plt.close(fig)

    for column, stem, ylabel in (
        ("golden_spiral_error", "base_vs_golden_spiral_error", "|fitted k - golden k|"),
        ("radial_spiral_fit_R2", "base_vs_spiral_fit_R2", "spiral fit R^2"),
        ("geometry_coherence_score", "base_vs_geometry_coherence_score", "heuristic coherence score"),
    ):
        fig, ax = plt.subplots(figsize=(10, 5)); ax.plot(metrics["base"], metrics[column], linewidth=.8)
        ax.axvline(PHI, color="#c83e4d", label="phi"); ax.set(xlabel="base", ylabel=ylabel, title=ylabel)
        ax.legend(); fig.tight_layout(); fig.savefig(OUTPUT / f"{stem}.png", dpi=170); plt.close(fig)
    zoom = metrics[metrics["base"].between(PHI - .21, PHI + .21)]
    fig, ax = plt.subplots(figsize=(9, 5)); ax.plot(zoom["base"], zoom["geometry_coherence_score"])
    ax.scatter(PHI_NEIGHBORS, [metrics.iloc[np.abs(metrics["base"] - value).argmin()]["geometry_coherence_score"] for value in PHI_NEIGHBORS])
    ax.axvline(PHI, color="#c83e4d"); ax.set(xlabel="base", ylabel="heuristic score", title="Phi-neighborhood scan")
    fig.tight_layout(); fig.savefig(OUTPUT / "phi_neighborhood_zoom.png", dpi=170); plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 5))
    controls.boxplot(column="geometry_coherence_score", by="control", ax=ax, rot=20)
    ax.set(title="Real phi score versus randomized controls", xlabel="control", ylabel="heuristic score")
    fig.suptitle(""); fig.tight_layout(); fig.savefig(OUTPUT / "real_vs_randomized_controls.png", dpi=170); plt.close(fig)

    sweep_bases = np.linspace(PHI - .20, PHI + .20, 45)
    fig = plt.figure(figsize=(6, 6)); ax = fig.add_subplot(projection="3d")
    def sweep(index: int):
        ax.clear(); frame = radial_coordinates(entities, float(sweep_bases[index]))
        ax.scatter(frame["radial_x"], frame["radial_y"], frame["z"], c=_colors(frame), s=20)
        ax.set(xlabel="x", ylabel="y", zlabel="log10(N)", title=f"Radial helicoid base {sweep_bases[index]:.5f}")
    animation = FuncAnimation(fig, sweep, frames=len(sweep_bases), interval=90)
    animation.save(OUTPUT / "animated_base_sweep_around_phi.gif", writer=PillowWriter(fps=10), dpi=90); plt.close(fig)

    fig = plt.figure(figsize=(6, 6)); ax = fig.add_subplot(projection="3d")
    ax.scatter(phi_frame["radial_x"], phi_frame["radial_y"], phi_frame["z"], c=_colors(phi_frame), s=24)
    ax.set(xlabel="x", ylabel="y", zlabel="log10(N)", title="Rotating radial helicoid: base phi")
    def rotate(index: int):
        ax.view_init(elev=25, azim=index * 6)
    animation = FuncAnimation(fig, rotate, frames=60, interval=70)
    animation.save(OUTPUT / "rotating_3d_phi.gif", writer=PillowWriter(fps=12), dpi=90); plt.close(fig)


def report(metrics: pd.DataFrame, controls: pd.DataFrame, stable: pd.DataFrame, fits: pd.DataFrame) -> str:
    phi_row = metrics.loc[metrics["base_name"].eq("phi")].iloc[0]
    rank_score = int(metrics["geometry_coherence_score"].rank(ascending=False, method="min")[phi_row.name])
    rank_error = int(metrics["golden_spiral_error"].rank(method="min")[phi_row.name])
    better = metrics.nsmallest(8, "golden_spiral_error")[["base", "golden_spiral_error"]]
    control_summary = controls[controls["control"].ne("observed")].groupby("control")["empirical_p_value_score_ge_observed"].first()
    def markdown_table(frame: pd.DataFrame) -> str:
        headers = frame.columns.tolist()
        lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
        for _, row in frame.iterrows():
            lines.append("| " + " | ".join(f"{value:.6g}" if isinstance(value, float) else str(value) for value in row) + " |")
        return "\n".join(lines)

    fit_rows = markdown_table(fits[["ordering", "fitted_k", "golden_spiral_relative_error", "spiral_fit_R2"]])
    return f"""# Golden Ratio Spiral Robustness Analysis

## Scope

The golden ratio is not assumed physically relevant. This analysis is falsification-oriented:
the goal is to test whether phi is special, not to prove it. The score is heuristic only and
must not be presented as physical evidence.

## Questions

1. **Does phi produce unusually coherent geometry?** Phi ranks {rank_score} of {len(metrics)}
   scanned bases by the declared heuristic score. This does not establish unusual coherence.
2. **Is the apparent spiral actually logarithmic?** Under the fixed `{DEFAULT_FIT_ORDERING}`
   ordering, the sparse radial projection has R^2 = {phi_row["radial_spiral_fit_R2"]:.4f}.
   The fit depends on ordering and should not be overinterpreted.
3. **Is it close to a golden spiral?** Phi's fitted k is {phi_row["fitted_k"]:.6f}; golden k is
   {GOLDEN_K:.6f}; relative error is {phi_row["golden_spiral_relative_error"]:.4f}. Phi ranks
   {rank_error} of {len(metrics)} by golden-spiral error, so the scan directly tests whether
   other bases fit at least as well.
4. **Does the pattern survive randomized controls?** The control table records empirical
   score p-values: {", ".join(f"{name}={value:.3f}" for name, value in control_summary.items())}.
   Similar randomized scores support an artifact interpretation.
5. **Does the result depend on stable particle handling?** Yes. See
   `data/golden_ratio_stable_mode_metrics.csv`; truncations are numerical assumptions, while
   infinity markers are excluded from finite logarithmic fits.
6. **Which bases outperform phi?** The smallest golden-spiral errors include:

{markdown_table(better)}

7. **Geometry, artifact, or potentially interesting signal?** The conservative interpretation
   is visualization geometry whose appearance is sensitive to base, point ordering, and stable
   handling. This exploratory analysis does not claim discovery or physical significance.

## Ordering Sensitivity At Phi

{fit_rows}

## Reproduce

```bash
python -m src.golden_ratio_spiral_analysis
```
"""


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    entities = load_entities()
    metrics = scan_base_metrics(entities)
    controls = control_metrics(entities)
    stable = stable_mode_metrics()
    fits = phi_ordering_fits(entities)
    metrics.to_csv(DATA / "golden_ratio_base_metrics.csv", index=False)
    controls.to_csv(DATA / "golden_ratio_control_comparison.csv", index=False)
    stable.to_csv(DATA / "golden_ratio_stable_mode_metrics.csv", index=False)
    fits.to_csv(DATA / "golden_ratio_phi_ordering_fits.csv", index=False)
    save_visualizations(entities, metrics, controls)
    (ROOT / "docs" / "golden_ratio_spiral_analysis.md").write_text(
        report(metrics, controls, stable, fits), encoding="ascii"
    )


if __name__ == "__main__":
    main()
