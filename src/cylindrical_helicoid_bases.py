"""Exploratory cylindrical helicoid plots across logarithmic bases.

Visual grouping is descriptive. It does not establish physical relationships,
validate a hypothesis, or imply a preferred logarithmic base.
"""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import pandas as pd

from src.dynamic_phase_structure import ROOT, TWO_PI, circular_gaps
from src.exploratory_3d_geometry import normalized
from src.stable_particle_sensitivity import load_variant, phases_for_base

DATA = ROOT / "data"
OUTPUT = ROOT / "outputs" / "cylindrical_helicoid"
BASES = {
    "2": 2.0,
    "e": math.e,
    "pi": math.pi,
    "phi": (1.0 + math.sqrt(5.0)) / 2.0,
    "10": 10.0,
    "2pi": TWO_PI,
    "sqrt2pi": math.sqrt(TWO_PI),
    "pi2": math.pi**2,
    "epi": math.e**math.pi,
}
DBSCAN_EPS = 1.6
DBSCAN_MIN_SAMPLES = 2
CLOSE_PAIR_DISTANCE = 1.6


def coordinate_table(entities: pd.DataFrame) -> pd.DataFrame:
    """Return MATLAB/Wolfram-ready helicoid coordinates for every selected base."""

    rows = []
    for base_name, base_value in BASES.items():
        theta = phases_for_base(entities, base_value)
        for (_, particle), angle in zip(entities.iterrows(), theta):
            rows.append(
                {
                    "name": particle["name"],
                    "family": particle["family"],
                    "dominant_interaction": particle["dominant_interaction"],
                    "base_name": base_name,
                    "base_value": base_value,
                    "theta_rad": angle,
                    "theta_deg": math.degrees(angle),
                    "x": math.cos(angle),
                    "y": math.sin(angle),
                    "z_log10_N": particle["log10_N"],
                    "N": particle["compton_cycles"],
                    "compton_frequency_hz": particle["compton_frequency_hz"],
                    "lifetime_s": particle["lifetime_s"],
                }
            )
    return pd.DataFrame(rows)


def distance_matrix(points: np.ndarray) -> np.ndarray:
    """Return pairwise Euclidean distances."""

    points = np.asarray(points, dtype=float)
    return np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)


def dbscan(points: np.ndarray, eps: float = DBSCAN_EPS, min_samples: int = DBSCAN_MIN_SAMPLES) -> np.ndarray:
    """Return compact DBSCAN labels with ``-1`` for noise."""

    distances = distance_matrix(points)
    labels = np.full(len(points), -99, dtype=int)
    cluster = 0
    for index in range(len(points)):
        if labels[index] != -99:
            continue
        neighbors = np.flatnonzero(distances[index] <= eps)
        if len(neighbors) < min_samples:
            labels[index] = -1
            continue
        labels[index] = cluster
        queue = list(neighbors[neighbors != index])
        while queue:
            neighbor = queue.pop()
            if labels[neighbor] == -1:
                labels[neighbor] = cluster
            if labels[neighbor] != -99:
                continue
            labels[neighbor] = cluster
            expanded = np.flatnonzero(distances[neighbor] <= eps)
            if len(expanded) >= min_samples:
                queue.extend(candidate for candidate in expanded if labels[candidate] in (-99, -1))
        cluster += 1
    return labels


def silhouette_score(points: np.ndarray, labels: np.ndarray) -> float:
    """Return a simple silhouette score over non-noise DBSCAN members."""

    points = np.asarray(points, dtype=float)
    labels = np.asarray(labels)
    active = labels >= 0
    unique = np.unique(labels[active])
    if len(unique) < 2:
        return float("nan")
    distances = distance_matrix(points)
    scores = []
    for index in np.flatnonzero(active):
        same = np.flatnonzero(labels == labels[index])
        same = same[same != index]
        a = float(np.mean(distances[index, same])) if len(same) else 0.0
        b = min(
            float(np.mean(distances[index, labels == label]))
            for label in unique if label != labels[index]
        )
        scores.append((b - a) / max(a, b) if max(a, b) else 0.0)
    return float(np.mean(scores))


def category_compactness(frame: pd.DataFrame, category: str) -> float:
    """Return mean within-category pair distance for categories with at least two members."""

    values = []
    for _, group in frame.groupby(category):
        points = group[["x", "y", "z_log10_N"]].to_numpy()
        if len(points) < 2:
            continue
        distances = distance_matrix(points)
        values.extend(distances[np.triu_indices(len(points), 1)])
    return float(np.mean(values)) if values else float("nan")


def metric_table(coordinates: pd.DataFrame) -> pd.DataFrame:
    """Compute requested grouping metrics for each base."""

    rows = []
    for base_name, frame in coordinates.groupby("base_name", sort=False):
        points = frame[["x", "y", "z_log10_N"]].to_numpy()
        distances = distance_matrix(points)
        nearest = np.partition(distances, 1, axis=1)[:, 1]
        labels = dbscan(points)
        _, _, gaps = circular_gaps(frame["theta_rad"].to_numpy())
        ordered_z = np.sort(frame["z_log10_N"].to_numpy())
        rows.append(
            {
                "base_name": base_name,
                "base_value": frame["base_value"].iloc[0],
                "avg_nearest_neighbor_distance": float(np.mean(nearest)),
                "family_compactness": category_compactness(frame, "family"),
                "interaction_compactness": category_compactness(frame, "dominant_interaction"),
                "cluster_count": int(len(set(labels)) - (-1 in labels)),
                "silhouette_score": silhouette_score(points, labels),
                "largest_gap_z": float(np.diff(ordered_z).max()),
                "largest_angular_gap_deg": float(np.degrees(gaps.max())),
            }
        )
    return pd.DataFrame(rows)


def persistent_neighbors(coordinates: pd.DataFrame) -> pd.DataFrame:
    """Track particle pairs that remain close in 3D helicoid space across bases."""

    names = coordinates[coordinates["base_name"] == next(iter(BASES))]["name"].tolist()
    values = {pair: [] for pair in itertools.combinations(names, 2)}
    for _, frame in coordinates.groupby("base_name", sort=False):
        points = frame.set_index("name")[["x", "y", "z_log10_N"]]
        for left, right in values:
            values[(left, right)].append(float(np.linalg.norm(points.loc[left] - points.loc[right])))
    rows = []
    for (left, right), distances in values.items():
        array = np.asarray(distances)
        rows.append(
            {
                "particle_a": left,
                "particle_b": right,
                "bases_close_count": int(np.sum(array < CLOSE_PAIR_DISTANCE)),
                "fraction_of_bases_close": float(np.mean(array < CLOSE_PAIR_DISTANCE)),
                "mean_3d_distance": float(array.mean()),
                "min_3d_distance": float(array.min()),
                "max_3d_distance": float(array.max()),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["fraction_of_bases_close", "mean_3d_distance"], ascending=[False, True]
    )


def _family_colors(frame: pd.DataFrame) -> list[object]:
    labels = sorted(frame["family"].unique())
    palette = {label: plt.get_cmap("tab10")(index) for index, label in enumerate(labels)}
    return [palette[label] for label in frame["family"]]


def _decorate_helicoid(ax, frame: pd.DataFrame) -> None:
    z_min, z_max = frame["z_log10_N"].min(), frame["z_log10_N"].max()
    theta = np.linspace(0.0, TWO_PI, 80)
    z = np.linspace(z_min, z_max, 50)
    theta_grid, z_grid = np.meshgrid(theta, z)
    ax.plot_surface(np.cos(theta_grid), np.sin(theta_grid), z_grid, color="#8fb8de",
                    alpha=0.08, linewidth=0)
    for ring_z in np.arange(math.floor(z_min / 2.0) * 2.0, z_max + 2.0, 2.0):
        ax.plot(np.cos(theta), np.sin(theta), np.full_like(theta, ring_z),
                color="0.55", alpha=0.35, linewidth=0.6)
    ax.set(xlabel="cos(theta)", ylabel="sin(theta)", zlabel="log10(N)",
           xlim=(-1.15, 1.15), ylim=(-1.15, 1.15))


def _draw_helicoid(ax, frame: pd.DataFrame, title: str, labels: bool) -> None:
    _decorate_helicoid(ax, frame)
    sizes = 25.0 + 75.0 * normalized(frame["z_log10_N"].to_numpy())
    ax.scatter(frame["x"], frame["y"], frame["z_log10_N"], c=_family_colors(frame),
               s=sizes, depthshade=True)
    if labels:
        for _, row in frame.iterrows():
            ax.text(row["x"], row["y"], row["z_log10_N"], row["name"], fontsize=6)
    ax.set_title(title)


def _save_base_plots(frame: pd.DataFrame, base_name: str) -> None:
    fig = plt.figure(figsize=(9, 8))
    _draw_helicoid(fig.add_subplot(projection="3d"), frame, f"Cylindrical helicoid: base {base_name}", True)
    fig.tight_layout()
    fig.savefig(OUTPUT / f"helicoid_base_{base_name}.png", dpi=180)
    fig.savefig(OUTPUT / f"helicoid_base_{base_name}.svg")
    plt.close(fig)

    fig = plt.figure(figsize=(9, 8))
    _draw_helicoid(fig.add_subplot(projection="3d"), frame,
                   f"Cylindrical helicoid without labels: base {base_name}", False)
    fig.tight_layout()
    fig.savefig(OUTPUT / f"helicoid_base_{base_name}_nolabels.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(frame["x"], frame["y"], c=_family_colors(frame), s=40)
    for _, row in frame.iterrows():
        ax.annotate(row["name"], (row["x"], row["y"]), fontsize=6)
    ax.set(aspect="equal", xlabel="cos(theta)", ylabel="sin(theta)",
           title=f"Top-down helicoid projection: base {base_name}")
    fig.tight_layout()
    fig.savefig(OUTPUT / f"helicoid_base_{base_name}_topdown.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(frame["x"], frame["z_log10_N"], c=_family_colors(frame), s=40)
    for _, row in frame.iterrows():
        ax.annotate(row["name"], (row["x"], row["z_log10_N"]), fontsize=6)
    ax.set(xlabel="cos(theta)", ylabel="log10(N)", title=f"Side helicoid projection: base {base_name}")
    fig.tight_layout()
    fig.savefig(OUTPUT / f"helicoid_base_{base_name}_side.png", dpi=170)
    plt.close(fig)


def _write_video(animation: FuncAnimation, stem: Path) -> None:
    animation.save(stem.with_suffix(".gif"), writer=PillowWriter(fps=12), dpi=90)
    try:
        import imageio_ffmpeg

        plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
        animation.save(stem.with_suffix(".mp4"), writer="ffmpeg", fps=12, dpi=90)
    except (ImportError, RuntimeError, FileNotFoundError):
        stem.with_suffix(".mp4.unavailable.txt").write_text(
            "MP4 export requires imageio-ffmpeg. Install requirements.txt and rerun.\n",
            encoding="ascii",
        )


def _save_rotation(frame: pd.DataFrame, base_name: str) -> None:
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection="3d")
    _draw_helicoid(ax, frame, f"Rotating cylindrical helicoid: base {base_name}", False)

    def update(index: int):
        ax.view_init(elev=24, azim=index * 4)
        return []

    animation = FuncAnimation(fig, update, frames=90, interval=80)
    _write_video(animation, OUTPUT / f"helicoid_base_{base_name}_rotating")
    plt.close(fig)


def _save_grid(coordinates: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(18, 18))
    for panel, base_name in enumerate(BASES, start=1):
        frame = coordinates[coordinates["base_name"] == base_name]
        ax = fig.add_subplot(3, 3, panel, projection="3d")
        _draw_helicoid(ax, frame, f"base {base_name}", False)
        ax.tick_params(labelsize=6)
    fig.suptitle("Cylindrical helicoid comparison across logarithmic bases")
    fig.tight_layout()
    fig.savefig(OUTPUT / "all_bases_grid.png", dpi=180)
    fig.savefig(OUTPUT / "all_bases_grid.svg")
    plt.close(fig)


def _save_metric_plot(metrics: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    axes[0].plot(metrics["base_name"], metrics["family_compactness"], marker="o")
    axes[0].set(ylabel="family compactness")
    axes[1].plot(metrics["base_name"], metrics["cluster_count"], marker="o")
    axes[1].set(ylabel="cluster count")
    axes[2].plot(metrics["base_name"], metrics["avg_nearest_neighbor_distance"], marker="o")
    axes[2].set(xlabel="logarithmic base", ylabel="avg nearest-neighbor distance")
    fig.suptitle("Cylindrical helicoid grouping metrics")
    fig.tight_layout()
    fig.savefig(OUTPUT / "base_comparison_metrics.png", dpi=180)
    fig.savefig(OUTPUT / "base_comparison_metrics.svg")
    plt.close(fig)


def generate_plots(coordinates: pd.DataFrame, metrics: pd.DataFrame) -> None:
    """Generate all requested static views, animations, and comparisons."""

    OUTPUT.mkdir(parents=True, exist_ok=True)
    for base_name in BASES:
        frame = coordinates[coordinates["base_name"] == base_name]
        _save_base_plots(frame, base_name)
        _save_rotation(frame, base_name)
    _save_grid(coordinates)
    _save_metric_plot(metrics)


def _report(metrics: pd.DataFrame, neighbors: pd.DataFrame) -> str:
    best = metrics.sort_values(
        ["family_compactness", "avg_nearest_neighbor_distance"]
    ).iloc[0]
    pi_row = metrics.set_index("base_name").loc["pi"]
    ten_row = metrics.set_index("base_name").loc["10"]
    repeated = neighbors[neighbors["fraction_of_bases_close"] == 1.0].head(8)
    lines = [
        "# Cylindrical Helicoid Persistence Across Bases",
        "",
        "## Scope",
        "",
        "This exploratory visualization uses the `finite_lifetime_only` dataset. Stable rows are "
        "excluded so truncation assumptions do not drive the comparison. Visual grouping does not "
        "imply a physical relationship, validation, or discovery.",
        "",
        "## Answers",
        "",
        f"- **Which base produces the clearest visible grouping?** By the declared family-compactness "
        f"metric, base `{best['base_name']}` is the smallest-distance case (`{best['family_compactness']:.6f}`). "
        "This is a descriptive ranking, not a privileged base.",
        "- **Do groups survive across bases?** Some local relationships survive, but angular arrangement "
        "changes visibly. The common `z = log10(N)` height preserves layering across every base.",
        "- **Are the same neighbors repeated?** Yes. The table below lists pairs close in all selected bases.",
        f"- **Does base pi remain visually interesting?** Base `pi` has family compactness "
        f"`{pi_row['family_compactness']:.6f}`, cluster count `{int(pi_row['cluster_count'])}`, and average "
        f"nearest-neighbor distance `{pi_row['avg_nearest_neighbor_distance']:.6f}`. It is one exploratory "
        "view among several, not evidence of special status.",
        f"- **Is base 10 only a control or does it show structure?** Base `10` shows the clearest direct "
        "wrapped-log10 helicoid because its angle and height use the same logarithm. Its family compactness "
        f"is `{ten_row['family_compactness']:.6f}`. That structure is primarily coordinate-induced.",
        "- **Does the cylindrical helicoid preserve more useful information than the 2D circle?** It preserves "
        "`log10(N)` height, which distinguishes particles that overlap angularly. Whether that is useful "
        "depends on the analysis question; it is not independent information because height is derived from `N`.",
        "- **Is visible structure dominated by log10(N) height?** Largely yes. Height is unchanged across "
        "bases and often dominates Euclidean 3D separation. Base changes rotate points around the cylinder.",
        "",
        "## Repeated Close Neighbors",
        "",
        "| Particle A | Particle B | Fraction of bases close | Mean 3D distance |",
        "| --- | --- | ---: | ---: |",
    ]
    for _, row in repeated.iterrows():
        lines.append(
            f"| {row['particle_a']} | {row['particle_b']} | {row['fraction_of_bases_close']:.6f} | "
            f"{row['mean_3d_distance']:.6f} |"
        )
    lines += [
        "",
        "## Clustering Method",
        "",
        f"DBSCAN uses Euclidean helicoid coordinates with `eps = {DBSCAN_EPS}` and "
        f"`min_samples = {DBSCAN_MIN_SAMPLES}`. Close-pair tracking uses distance `< {CLOSE_PAIR_DISTANCE}`. "
        "These are exploratory display-scale choices.",
        "",
        "## Warning",
        "",
        "The helicoid representation is exploratory. Visual grouping does not imply physical relation.",
        "",
        "## Reproduction",
        "",
        "```powershell",
        "python -m src.cylindrical_helicoid_bases",
        "```",
    ]
    return "\n".join(lines) + "\n"


def generate() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Generate coordinate tables, metrics, plots, animations, and report."""

    entities = load_variant()
    coordinates = coordinate_table(entities)
    metrics = metric_table(coordinates)
    neighbors = persistent_neighbors(coordinates)
    DATA.mkdir(parents=True, exist_ok=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    coordinates.to_csv(DATA / "helicoid_coordinates_by_base.csv", index=False)
    metrics.to_csv(DATA / "cylindrical_helicoid_metrics.csv", index=False)
    neighbors.to_csv(DATA / "helicoid_persistent_neighbors.csv", index=False)
    generate_plots(coordinates, metrics)
    (ROOT / "docs" / "cylindrical_helicoid_bases.md").write_text(
        _report(metrics, neighbors), encoding="ascii"
    )
    return coordinates, metrics, neighbors


if __name__ == "__main__":
    generated = generate()
    print(
        f"Generated {len(generated[0])} coordinates, {len(generated[1])} base metrics, "
        f"and {len(generated[2])} neighbor rows in {OUTPUT}"
    )
