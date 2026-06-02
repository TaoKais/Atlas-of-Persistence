"""Continuous logarithmic-base scan for exploratory persistence phases.

The metrics in this module describe a finite dataset under a chosen circular
projection. They do not establish a preferred base, a physical resonance, or
new physics.
"""

from __future__ import annotations

import itertools
import json
import math
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import pandas as pd

from src.gap_geometry import polygon_metrics
from src.persistence import compton_frequency

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUT = ROOT / "outputs" / "base_scan"
TWO_PI = 2.0 * math.pi
ALPHA = 7.297_352_5643e-3
PHI = (1.0 + math.sqrt(5.0)) / 2.0
BASE_SAMPLES = 5000
BASE_MIN = 1.2
BASE_MAX = 50.0
NEIGHBOR_THRESHOLD_DEG = 20.0
CLUSTER_THRESHOLD_DEG = 28.0
ENTROPY_BINS = 12
SPECIAL_BASES = {
    "2": 2.0,
    "e": math.e,
    "pi": math.pi,
    "phi": PHI,
    "2*pi": TWO_PI,
    "sqrt(2*pi)": math.sqrt(TWO_PI),
    "pi^2": math.pi**2,
    "e^pi": math.e**math.pi,
    "10": 10.0,
    "1/alpha": 1.0 / ALPHA,
}
FAMILIES = ("lepton", "meson", "baryon", "boson", "quark")
INTERACTIONS = ("weak", "strong", "electromagnetic", "mixed", "stable")


def _json(values: list[object] | np.ndarray) -> str:
    return json.dumps([round(float(value), 12) for value in values])


def _circular_distance_deg(left: np.ndarray | float, right: np.ndarray | float) -> np.ndarray:
    return np.abs((np.asarray(left) - np.asarray(right) + 180.0) % 360.0 - 180.0)


def load_entities(
    particle_path: Path = DATA / "particles.csv",
    persistence_path: Path = DATA / "persistence_particles.csv",
) -> pd.DataFrame:
    """Load finite lifetimes and substitute documented stable-particle bounds."""

    particles = pd.read_csv(particle_path)
    bounds = (
        pd.read_csv(persistence_path)
        .rename(columns={"particle": "name", "lifetime_s": "bound_lifetime_s"})
        [["name", "bound_lifetime_s", "lifetime_qualifier", "lifetime_basis"]]
    )
    entities = particles.merge(bounds, on="name", how="left")
    entities["lifetime_s"] = pd.to_numeric(entities["lifetime_s"], errors="coerce")
    entities["bound_lifetime_s"] = pd.to_numeric(entities["bound_lifetime_s"], errors="coerce")
    stable_bound = entities["lifetime_s"].isna() & entities["bound_lifetime_s"].gt(0.0)
    entities.loc[stable_bound, "lifetime_s"] = entities.loc[stable_bound, "bound_lifetime_s"]
    entities["lifetime_handling"] = np.where(stable_bound, "documented_lower_bound", "finite_value")
    entities = entities[entities["lifetime_s"].gt(0.0) & entities["mass_mev"].gt(0.0)].copy()
    entities["compton_frequency_hz"] = entities["mass_mev"].map(compton_frequency)
    entities["compton_cycles"] = entities["compton_frequency_hz"] * entities["lifetime_s"]
    entities["ln_N"] = np.log(entities["compton_cycles"])
    return entities.reset_index(drop=True)


def phase_matrix(entities: pd.DataFrame, bases: np.ndarray) -> np.ndarray:
    """Return particle phases in radians with shape ``(number_of_bases, entities)``."""

    return TWO_PI * np.mod(entities["ln_N"].to_numpy()[None, :] / np.log(bases)[:, None], 1.0)


def circular_gaps(theta: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return sorted positions, source indices, and adjacent circular gaps."""

    order = np.argsort(theta)
    ordered = theta[order]
    gaps = np.diff(np.r_[ordered, ordered[0] + TWO_PI])
    return ordered, order, gaps


def circular_clusters(theta: np.ndarray, threshold_deg: float = CLUSTER_THRESHOLD_DEG) -> list[list[int]]:
    """Split circularly ordered points at gaps larger than ``threshold_deg``."""

    _, order, gaps = circular_gaps(theta)
    breaks = np.flatnonzero(np.degrees(gaps) > threshold_deg)
    if not len(breaks):
        return [order.tolist()]
    start = (int(breaks[0]) + 1) % len(order)
    rotated = np.roll(order, -start)
    rotated_gaps = np.roll(gaps, -start)
    split_after = np.flatnonzero(np.degrees(rotated_gaps[:-1]) > threshold_deg) + 1
    return [part.tolist() for part in np.split(rotated, split_after) if len(part)]


def _entropy(theta: np.ndarray) -> float:
    counts, _ = np.histogram(theta, bins=ENTROPY_BINS, range=(0.0, TWO_PI))
    probabilities = counts[counts > 0] / counts.sum()
    return float(-np.sum(probabilities * np.log2(probabilities)))


def _gap_geometry(
    gaps: np.ndarray, ordered: np.ndarray, order: np.ndarray, names: list[str], top_k: int
) -> dict[str, object]:
    selected = np.argsort(gaps)[::-1][:top_k]
    midpoints = np.mod(ordered[selected] + gaps[selected] / 2.0, TWO_PI)
    metrics = polygon_metrics(np.degrees(midpoints))
    labels = [
        f"{names[order[index]]} -> {names[order[(index + 1) % len(order)]]}"
        for index in selected
    ]
    return {
        "top_k": top_k,
        "selected_gaps": json.dumps(labels),
        "midpoint_angles_deg": _json(np.degrees(np.sort(midpoints))),
        "midpoint_x": _json(np.cos(np.sort(midpoints))),
        "midpoint_y": _json(np.sin(np.sort(midpoints))),
        "area": metrics["polygon_area"],
        "centroid_x": metrics["centroid_x"],
        "centroid_y": metrics["centroid_y"],
        "centroid_radius": metrics["centroid_radius"],
        "side_lengths": _json(metrics["side_lengths"]),
        "angular_spacing_deg": _json(metrics["angle_gaps_between_vertices"]),
        "symmetry_score": metrics["circular_symmetry_score"],
    }


def scan_metrics(
    entities: pd.DataFrame, bases: np.ndarray, phases: np.ndarray
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute base-level circular, gap, cluster, and polygon metrics."""

    metric_rows: list[dict[str, object]] = []
    geometry_rows: list[dict[str, object]] = []
    names = entities["name"].tolist()
    for index, (base, theta) in enumerate(zip(bases, phases)):
        ordered, order, gaps = circular_gaps(theta)
        gaps_deg = np.degrees(gaps)
        clusters = circular_clusters(theta)
        z_mean = np.mean(np.exp(1j * theta))
        geometries = [
            _gap_geometry(gaps, ordered, order, names, top_k) for top_k in (3, 4, 5, 6)
        ]
        top4 = next(row for row in geometries if row["top_k"] == 4)
        largest_gap_index = int(np.argmax(gaps))
        metric_rows.append(
            {
                "base_index": index,
                "base": base,
                "circular_mean_rad": float(np.mod(np.angle(z_mean), TWO_PI)),
                "circular_mean_deg": float(np.degrees(np.mod(np.angle(z_mean), TWO_PI))),
                "circular_variance": 1.0 - abs(z_mean),
                "resultant_length_R": abs(z_mean),
                "entropy_bits": _entropy(theta),
                "largest_gap_deg": float(np.max(gaps_deg)),
                "largest_gap_from": names[order[largest_gap_index]],
                "largest_gap_to": names[order[(largest_gap_index + 1) % len(order)]],
                "top_3_gap_sum_deg": float(np.sum(np.sort(gaps_deg)[-3:])),
                "top_5_gap_sum_deg": float(np.sum(np.sort(gaps_deg)[-5:])),
                "gap_distribution_deg": _json(np.sort(gaps_deg)[::-1]),
                "cluster_count": len(clusters),
                "mean_cluster_density": float(np.mean([len(cluster) for cluster in clusters])),
                "isolated_points": sum(len(cluster) == 1 for cluster in clusters),
                "top4_symmetry_score": top4["symmetry_score"],
            }
        )
        for geometry in geometries:
            geometry_rows.append({"base_index": index, "base": base, **geometry})
    return pd.DataFrame(metric_rows), pd.DataFrame(geometry_rows)


def _category_alignment(
    entities: pd.DataFrame,
    bases: np.ndarray,
    phases: np.ndarray,
    category_column: str,
    categories: tuple[str, ...],
) -> pd.DataFrame:
    rows = []
    values = entities[category_column].to_numpy()
    for base_index, (base, theta) in enumerate(zip(bases, phases)):
        theta_deg = np.degrees(theta)
        for category in categories:
            mask = values == category
            if not np.any(mask):
                continue
            selected = theta[mask]
            z_mean = np.mean(np.exp(1j * selected))
            outside = theta_deg[~mask]
            overlap = (
                float(np.mean(np.min(_circular_distance_deg(theta_deg[mask, None], outside), axis=1)
                              < NEIGHBOR_THRESHOLD_DEG))
                if len(outside)
                else 0.0
            )
            rows.append(
                {
                    "base_index": base_index,
                    "base": base,
                    "category_type": category_column,
                    "category": category,
                    "particle_count": int(np.sum(mask)),
                    "angular_centroid_deg": float(np.degrees(np.mod(np.angle(z_mean), TWO_PI))),
                    "angular_variance": 1.0 - abs(z_mean),
                    "overlap_score": overlap,
                }
            )
    return pd.DataFrame(rows)


def alignment_metrics(entities: pd.DataFrame, bases: np.ndarray, phases: np.ndarray) -> pd.DataFrame:
    """Compute family and dominant-interaction alignment tables."""

    return pd.concat(
        [
            _category_alignment(entities, bases, phases, "family", FAMILIES),
            _category_alignment(
                entities, bases, phases, "dominant_interaction", INTERACTIONS
            ),
        ],
        ignore_index=True,
    )


def _neighbor_masks(phases: np.ndarray) -> dict[tuple[int, int], np.ndarray]:
    masks = {}
    for left, right in itertools.combinations(range(phases.shape[1]), 2):
        masks[(left, right)] = _circular_distance_deg(
            np.degrees(phases[:, left]), np.degrees(phases[:, right])
        ) < NEIGHBOR_THRESHOLD_DEG
    return masks


def persistent_relationships(
    entities: pd.DataFrame, bases: np.ndarray, phases: np.ndarray
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return persistent pairs, triplets, recurring clusters, and particle stability."""

    names = entities["name"].tolist()
    masks = _neighbor_masks(phases)
    rows: list[dict[str, object]] = []
    for (left, right), mask in masks.items():
        rows.append(
            {
                "record_type": "pair",
                "members": json.dumps([names[left], names[right]]),
                "member_count": 2,
                "persistence_neighbor_score": float(np.mean(mask)),
            }
        )
    for members in itertools.combinations(range(len(names)), 3):
        mask = masks[(members[0], members[1])] & masks[(members[0], members[2])] & masks[
            (members[1], members[2])
        ]
        score = float(np.mean(mask))
        if score > 0.0:
            rows.append(
                {
                    "record_type": "triplet",
                    "members": json.dumps([names[index] for index in members]),
                    "member_count": 3,
                    "persistence_neighbor_score": score,
                }
            )
    recurring: Counter[tuple[str, ...]] = Counter()
    for theta in phases:
        for cluster in circular_clusters(theta, NEIGHBOR_THRESHOLD_DEG):
            if len(cluster) >= 3:
                recurring[tuple(sorted(names[index] for index in cluster))] += 1
    for members, count in recurring.items():
        rows.append(
            {
                "record_type": "cluster",
                "members": json.dumps(members),
                "member_count": len(members),
                "persistence_neighbor_score": count / len(phases),
            }
        )
    relationships = pd.DataFrame(rows).sort_values(
        ["record_type", "persistence_neighbor_score"], ascending=[True, False]
    )

    unwrapped = np.unwrap(phases, axis=0)
    velocity = np.gradient(unwrapped, bases, axis=0)
    rank_history = np.argsort(np.argsort(phases, axis=1), axis=1)
    stability_rows = []
    for particle_index, name in enumerate(names):
        pair_scores = [
            float(np.mean(mask))
            for members, mask in masks.items()
            if particle_index in members
        ]
        stability_rows.append(
            {
                "name": name,
                "angular_drift_deg": float(np.degrees(np.ptp(unwrapped[:, particle_index]))),
                "mean_abs_angular_velocity_deg_per_base": float(
                    np.degrees(np.mean(np.abs(velocity[:, particle_index])))
                ),
                "max_abs_angular_velocity_deg_per_base": float(
                    np.degrees(np.max(np.abs(velocity[:, particle_index])))
                ),
                "crossing_events": int(np.count_nonzero(np.diff(rank_history[:, particle_index]))),
                "neighborhood_persistence": float(np.mean(pair_scores)),
            }
        )
    return relationships.reset_index(drop=True), pd.DataFrame(stability_rows)


def special_base_metrics(
    entities: pd.DataFrame, scan: pd.DataFrame
) -> pd.DataFrame:
    """Evaluate named special bases, including values outside the scan interval."""

    bases = np.array(list(SPECIAL_BASES.values()))
    special_scan, _ = scan_metrics(entities, bases, phase_matrix(entities, bases))
    special_scan.insert(0, "special_base", list(SPECIAL_BASES))
    scan_bases = scan["base"].to_numpy()
    special_scan["within_continuous_scan"] = special_scan["base"].between(BASE_MIN, BASE_MAX)
    special_scan["nearest_scan_base"] = [
        scan_bases[np.argmin(np.abs(scan_bases - base))] if base <= BASE_MAX else np.nan
        for base in bases
    ]
    special_scan["neighbor_largest_gap_percentile"] = [
        _local_percentile(scan, base, "largest_gap_deg") for base in bases
    ]
    special_scan["neighbor_symmetry_percentile"] = [
        _local_percentile(scan, base, "top4_symmetry_score") for base in bases
    ]
    return special_scan


def _local_percentile(scan: pd.DataFrame, base: float, column: str, radius: float = 0.5) -> float:
    local = scan[np.abs(scan["base"] - base) <= radius]
    if local.empty:
        return float("nan")
    nearest = scan.iloc[(scan["base"] - base).abs().argmin()][column]
    return float(np.mean(local[column] <= nearest))


def _save_line(scan: pd.DataFrame, column: str, ylabel: str, filename: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(scan["base"], scan[column], linewidth=1.0)
    ax.set(xlabel="logarithmic base", ylabel=ylabel, title=f"Base vs {ylabel}")
    fig.tight_layout()
    fig.savefig(OUTPUT / filename, dpi=160)
    plt.close(fig)


def _save_relationship_network(entities: pd.DataFrame, relationships: pd.DataFrame) -> None:
    pairs = relationships[relationships["record_type"] == "pair"].nlargest(
        18, "persistence_neighbor_score"
    )
    names = sorted(set(itertools.chain.from_iterable(pairs["members"].map(json.loads))))
    angles = np.linspace(0.0, TWO_PI, len(names), endpoint=False)
    positions = dict(zip(names, zip(np.cos(angles), np.sin(angles))))
    fig, ax = plt.subplots(figsize=(8, 8))
    for _, row in pairs.iterrows():
        left, right = json.loads(row["members"])
        x = [positions[left][0], positions[right][0]]
        y = [positions[left][1], positions[right][1]]
        ax.plot(x, y, color="#2457a6", alpha=0.15 + 0.8 * row["persistence_neighbor_score"],
                linewidth=0.7 + 4.0 * row["persistence_neighbor_score"])
    for name, (x, y) in positions.items():
        ax.scatter(x, y, color="#c83e4d", s=45, zorder=3)
        ax.annotate(name, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.set(aspect="equal", title="Persistent neighbor network: top 18 pairs")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(OUTPUT / "persistent_neighbor_network.png", dpi=180)
    plt.close(fig)


def _save_centroid_evolution(alignment: pd.DataFrame, category_type: str, filename: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    frame = alignment[alignment["category_type"] == category_type]
    for category, group in frame.groupby("category"):
        ax.plot(group["base"], group["angular_centroid_deg"], label=category, linewidth=0.9)
    ax.set(xlabel="logarithmic base", ylabel="angular centroid (degrees)",
           title=f"{category_type.replace('_', ' ').title()} centroid evolution", ylim=(0, 360))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT / filename, dpi=160)
    plt.close(fig)


def _save_structured_bases(scan: pd.DataFrame, largest: bool, filename: str) -> None:
    frame = scan.nlargest(10, "top4_symmetry_score") if largest else scan.nsmallest(
        10, "top4_symmetry_score"
    )
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.barh(frame["base"].map(lambda value: f"{value:.5f}"), frame["top4_symmetry_score"])
    ax.set(xlabel="top-4 symmetry score", ylabel="base",
           title=f"Top 10 {'most' if largest else 'least'} structured bases")
    fig.tight_layout()
    fig.savefig(OUTPUT / filename, dpi=160)
    plt.close(fig)


def _animation(
    entities: pd.DataFrame, bases: np.ndarray, phases: np.ndarray, geometry: pd.DataFrame
) -> None:
    frame_indices = np.linspace(0, len(bases) - 1, 120, dtype=int)
    family_colors = {
        family: plt.get_cmap("tab10")(index)
        for index, family in enumerate(sorted(entities["family"].unique()))
    }
    colors = [family_colors[family] for family in entities["family"]]
    fig, ax = plt.subplots(figsize=(7, 7))
    circle = np.linspace(0.0, TWO_PI, 400)
    ax.plot(np.cos(circle), np.sin(circle), color="0.7")
    initial_theta = phases[frame_indices[0]]
    points = ax.scatter(np.cos(initial_theta), np.sin(initial_theta), c=colors, s=34)
    polygon, = ax.plot([], [], color="#c83e4d", linewidth=1.8)
    gap_lines = [ax.plot([], [], color="#c83e4d", alpha=0.5)[0] for _ in range(4)]
    title = ax.set_title("")
    ax.set(aspect="equal", xlim=(-1.15, 1.15), ylim=(-1.15, 1.15),
           xlabel="cos(theta)", ylabel="sin(theta)")

    def update(frame_number: int):
        base_index = frame_indices[frame_number]
        theta = phases[base_index]
        points.set_offsets(np.c_[np.cos(theta), np.sin(theta)])
        row = geometry[(geometry["base_index"] == base_index) & (geometry["top_k"] == 4)].iloc[0]
        polygon_theta = np.radians(json.loads(row["midpoint_angles_deg"]))
        x = np.cos(polygon_theta)
        y = np.sin(polygon_theta)
        polygon.set_data(np.r_[x, x[0]], np.r_[y, y[0]])
        ordered, _, gaps = circular_gaps(theta)
        for line, gap_index in zip(gap_lines, np.argsort(gaps)[::-1][:4]):
            start = ordered[gap_index]
            end = start + gaps[gap_index]
            arc = np.linspace(start, end, 30)
            line.set_data(np.cos(arc), np.sin(arc))
        title.set_text(f"Dynamic persistence phase: base {bases[base_index]:.4f}")
        return [points, polygon, title, *gap_lines]

    animation = FuncAnimation(fig, update, frames=len(frame_indices), interval=90, blit=False)
    animation.save(OUTPUT / "dynamic_phase_scan.gif", writer=PillowWriter(fps=12), dpi=100)
    try:
        import imageio_ffmpeg

        plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
        animation.save(OUTPUT / "dynamic_phase_scan.mp4", writer="ffmpeg", fps=12, dpi=100)
    except (ImportError, RuntimeError, FileNotFoundError):
        (OUTPUT / "dynamic_phase_scan.mp4.unavailable.txt").write_text(
            "MP4 export requires imageio-ffmpeg. Install requirements.txt and rerun.\n",
            encoding="ascii",
        )
    plt.close(fig)


def generate_plots(
    entities: pd.DataFrame,
    bases: np.ndarray,
    phases: np.ndarray,
    scan: pd.DataFrame,
    geometry: pd.DataFrame,
    alignment: pd.DataFrame,
    relationships: pd.DataFrame,
) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    _save_line(scan, "largest_gap_deg", "largest gap (degrees)", "base_vs_largest_gap.png")
    _save_line(scan, "entropy_bits", "entropy (bits)", "base_vs_entropy.png")
    _save_line(scan, "top4_symmetry_score", "top-4 symmetry score", "base_vs_symmetry_score.png")
    _save_line(scan, "resultant_length_R", "resultant length R", "base_vs_resultant_length.png")
    _save_relationship_network(entities, relationships)
    _save_centroid_evolution(alignment, "family", "family_centroid_evolution.png")
    _save_centroid_evolution(
        alignment, "dominant_interaction", "interaction_centroid_evolution.png"
    )
    _save_structured_bases(scan, True, "top_10_most_structured_bases.png")
    _save_structured_bases(scan, False, "top_10_least_structured_bases.png")
    _animation(entities, bases, phases, geometry)


def _report(
    entities: pd.DataFrame,
    scan: pd.DataFrame,
    special: pd.DataFrame,
    relationships: pd.DataFrame,
) -> str:
    pair_rows = relationships[relationships["record_type"] == "pair"].nlargest(
        8, "persistence_neighbor_score"
    )
    strongest = scan.nlargest(1, "top4_symmetry_score").iloc[0]
    weakest = scan.nsmallest(1, "top4_symmetry_score").iloc[0]
    selected = special.set_index("special_base")
    lines = [
        "# Dynamic Phase Structure Scan",
        "",
        "## Scope",
        "",
        f"This exploratory scan evaluates {len(entities)} particles across {len(scan)} evenly spaced "
        f"bases from `{BASE_MIN}` to `{BASE_MAX}`. Electron and proton use documented lifetime lower "
        "bounds. Their phases are therefore bound-based coordinates, not measured total-lifetime phases.",
        "",
        "The projection is descriptive. It does not establish a preferred logarithmic base, resonance "
        "mechanism, physical significance, or new physics.",
        "",
        "## Answers",
        "",
        f"- **Is pi genuinely special?** No statistically justified special status is established. Its "
        f"local top-4 symmetry percentile is `{selected.loc['pi', 'neighbor_symmetry_percentile']:.3f}`.",
        f"- **Is phi genuinely special?** No statistically justified special status is established. Its "
        f"local top-4 symmetry percentile is `{selected.loc['phi', 'neighbor_symmetry_percentile']:.3f}`.",
        f"- **Is base 2 genuinely special?** No statistically justified special status is established. Its "
        f"local top-4 symmetry percentile is `{selected.loc['2', 'neighbor_symmetry_percentile']:.3f}`.",
        "- **Are observed structures robust?** The metrics vary repeatedly under continuous base changes. "
        "Any selected-base geometry must be interpreted against the full scan.",
        "- **Are the largest gaps stable?** No. Their magnitudes and bounding particles change across the scan.",
        "- **Do persistent clusters exist?** Recurring neighbor groups exist descriptively at the declared "
        f"`{NEIGHBOR_THRESHOLD_DEG}` degree threshold; they are dataset relationships, not physical resonances.",
        f"- **Does any base maximize organization beyond neighboring bases?** The strongest sampled top-4 "
        f"symmetry occurs at `{strongest['base']:.8f}` with score `{strongest['top4_symmetry_score']:.6f}`. "
        "This is a scan maximum, not evidence of significance; nearby-base and null-model testing would be "
        "required for a stronger claim.",
        "",
        "## Scan Extremes",
        "",
        f"- Highest sampled top-4 symmetry: base `{strongest['base']:.8f}`, score "
        f"`{strongest['top4_symmetry_score']:.6f}`.",
        f"- Lowest sampled top-4 symmetry: base `{weakest['base']:.8f}`, score "
        f"`{weakest['top4_symmetry_score']:.6f}`.",
        "",
        "## Most Persistent Neighbor Pairs",
        "",
        "| Pair | Fraction of bases |",
        "| --- | ---: |",
    ]
    for _, row in pair_rows.iterrows():
        lines.append(f"| {', '.join(json.loads(row['members']))} | {row['persistence_neighbor_score']:.6f} |")
    lines += [
        "",
        "## Bonus Projection Comparison",
        "",
        "`theta = 2*pi*frac(ln(N))` and `theta = arg(exp(i*ln(N)))` are mathematically equivalent "
        "up to angle wrapping. They are also exactly the base-`e` member of "
        "`theta = 2*pi*frac(log_b(N))`. Circular geometry is therefore caused by the circular projection; "
        "its detailed layout changes with logarithmic base.",
        "",
        "## Reproduction",
        "",
        "```powershell",
        "python -m src.dynamic_phase_structure",
        "```",
    ]
    return "\n".join(lines) + "\n"


def generate() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Generate CSV tables, figures, animations, and the report."""

    entities = load_entities()
    bases = np.linspace(BASE_MIN, BASE_MAX, BASE_SAMPLES)
    phases = phase_matrix(entities, bases)
    scan, geometry = scan_metrics(entities, bases, phases)
    alignment = alignment_metrics(entities, bases, phases)
    relationships, stability = persistent_relationships(entities, bases, phases)
    special = special_base_metrics(entities, scan)
    DATA.mkdir(parents=True, exist_ok=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    scan.to_csv(DATA / "base_scan_metrics.csv", index=False)
    geometry.to_csv(DATA / "gap_geometry_metrics.csv", index=False)
    alignment.to_csv(DATA / "family_alignment_metrics.csv", index=False)
    relationships.to_csv(DATA / "persistent_neighbors.csv", index=False)
    stability.to_csv(DATA / "particle_phase_stability.csv", index=False)
    special.to_csv(DATA / "special_base_metrics.csv", index=False)
    generate_plots(entities, bases, phases, scan, geometry, alignment, relationships)
    (ROOT / "docs" / "dynamic_phase_structure_report.md").write_text(
        _report(entities, scan, special, relationships), encoding="ascii"
    )
    return scan, geometry, alignment, relationships


if __name__ == "__main__":
    generated = generate()
    print(
        f"Generated {len(generated[0])} base-scan rows and {len(generated[1])} geometry rows "
        f"in {OUTPUT}"
    )
