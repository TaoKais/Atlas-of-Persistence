"""Exploratory harmonic matching and 3D geometry for persistence-phase peaks.

Close numerical matches are descriptive coincidences unless independently
validated against an explicit null model. This module does not claim physical
significance or a preferred logarithmic base.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import pandas as pd

from src.dynamic_phase_structure import (
    BASE_MAX,
    BASE_MIN,
    BASE_SAMPLES,
    DATA,
    NEIGHBOR_THRESHOLD_DEG,
    PHI,
    ROOT,
    TWO_PI,
    circular_gaps,
    load_entities,
    persistent_relationships,
    phase_matrix,
)
from src.gap_geometry import polygon_metrics

OUTPUT = ROOT / "outputs" / "harmonic_base_geometry"
SELECTED_OUTPUT = OUTPUT / "selected_bases"
TOP_BASES = 20


def detect_peaks(scan: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Detect local extrema and estimate peak prominence and half-height width."""

    bases = scan["base"].to_numpy(dtype=float)
    scores = scan["top4_symmetry_score"].to_numpy(dtype=float)
    maxima = np.flatnonzero((scores[1:-1] > scores[:-2]) & (scores[1:-1] >= scores[2:])) + 1
    minima = np.flatnonzero((scores[1:-1] < scores[:-2]) & (scores[1:-1] <= scores[2:])) + 1
    peak_rows = []
    for index in maxima:
        left_minima = minima[minima < index]
        right_minima = minima[minima > index]
        left = int(left_minima[-1]) if len(left_minima) else 0
        right = int(right_minima[0]) if len(right_minima) else len(scores) - 1
        baseline = max(scores[left], scores[right])
        prominence = max(0.0, scores[index] - baseline)
        half_height = scores[index] - prominence / 2.0
        width_left = index
        width_right = index
        while width_left > left and scores[width_left] > half_height:
            width_left -= 1
        while width_right < right and scores[width_right] > half_height:
            width_right += 1
        peak_rows.append(
            {
                "base": bases[index],
                "score": scores[index],
                "prominence": prominence,
                "width": bases[width_right] - bases[width_left],
                "base_index": index,
            }
        )
    peaks = pd.DataFrame(peak_rows).sort_values(
        ["prominence", "score"], ascending=False
    ).reset_index(drop=True)
    peaks.insert(0, "rank", np.arange(1, len(peaks) + 1))
    minima_frame = pd.DataFrame(
        {
            "base": bases[minima],
            "score": scores[minima],
            "base_index": minima,
        }
    ).sort_values("score").reset_index(drop=True)
    return peaks[["base", "score", "prominence", "width", "rank"]], minima_frame


def harmonic_constants() -> pd.DataFrame:
    """Return the requested simple harmonic constant candidates."""

    rows = []
    for n in range(1, 21):
        candidates = {
            "n*pi": n * math.pi,
            "n*e": n * math.e,
            "n*phi": n * PHI,
            "pi^n": math.pi**n,
            "e^n": math.e**n,
            "phi^n": PHI**n,
            "sqrt(n*pi)": math.sqrt(n * math.pi),
        }
        rows.extend(
            {"expression": expression, "n": n, "constant_value": value}
            for expression, value in candidates.items()
        )
    return pd.DataFrame(rows)


def constant_matches(peaks: pd.DataFrame) -> pd.DataFrame:
    """Compare every detected peak with every requested harmonic candidate."""

    constants = harmonic_constants()
    rows = []
    for _, peak in peaks.iterrows():
        for _, candidate in constants.iterrows():
            absolute_error = abs(peak["base"] - candidate["constant_value"])
            rows.append(
                {
                    "peak_rank": int(peak["rank"]),
                    "peak_base": peak["base"],
                    "peak_score": peak["score"],
                    "expression": candidate["expression"],
                    "n": int(candidate["n"]),
                    "constant_value": candidate["constant_value"],
                    "absolute_error": absolute_error,
                    "relative_error": absolute_error / peak["base"],
                }
            )
    matches = pd.DataFrame(rows).sort_values(
        ["peak_rank", "relative_error", "absolute_error"]
    ).reset_index(drop=True)
    matches["match_rank_for_peak"] = matches.groupby("peak_rank").cumcount() + 1
    return matches


def fft_spectrum(scan: pd.DataFrame) -> pd.DataFrame:
    """Compute a one-sided FFT power spectrum for the evenly sampled score signal."""

    bases = scan["base"].to_numpy(dtype=float)
    score = scan["top4_symmetry_score"].to_numpy(dtype=float)
    centered = score - np.mean(score)
    window = np.hanning(len(centered))
    fft = np.fft.rfft(centered * window)
    frequency = np.fft.rfftfreq(len(centered), d=float(bases[1] - bases[0]))
    power = np.abs(fft) ** 2
    return pd.DataFrame(
        {
            "frequency_cycles_per_base": frequency,
            "period_base_units": np.divide(
                1.0, frequency, out=np.full_like(frequency, np.nan), where=frequency > 0.0
            ),
            "power": power,
            "normalized_power": power / (power.max() or 1.0),
        }
    )


def gap_geometry_records(
    entities: pd.DataFrame, selected: pd.DataFrame
) -> tuple[pd.DataFrame, dict[tuple[int, int], np.ndarray]]:
    """Compute top-3 through top-6 gap-center polygon metrics for selected bases."""

    names = entities["name"].tolist()
    phases = phase_matrix(entities, selected["base"].to_numpy(dtype=float))
    rows = []
    midpoints_by_key = {}
    for (_, peak), theta in zip(selected.iterrows(), phases):
        ordered, order, gaps = circular_gaps(theta)
        for top_k in (3, 4, 5, 6):
            indices = np.argsort(gaps)[::-1][:top_k]
            midpoints = np.mod(ordered[indices] + gaps[indices] / 2.0, TWO_PI)
            metrics = polygon_metrics(np.degrees(midpoints))
            midpoints_by_key[(int(peak["rank"]), top_k)] = np.sort(midpoints)
            labels = [
                f"{names[order[index]]} -> {names[order[(index + 1) % len(order)]]}"
                for index in indices
            ]
            perimeter = float(np.sum(metrics["side_lengths"]))
            area = float(metrics["polygon_area"])
            rows.append(
                {
                    "peak_rank": int(peak["rank"]),
                    "base": peak["base"],
                    "top_k": top_k,
                    "selected_gaps": json.dumps(labels),
                    "midpoint_angles_deg": json.dumps(np.degrees(np.sort(midpoints)).tolist()),
                    "area": area,
                    "perimeter": perimeter,
                    "symmetry_score": metrics["circular_symmetry_score"],
                    "compactness": 4.0 * math.pi * area / perimeter**2 if perimeter else 0.0,
                }
            )
    return pd.DataFrame(rows), midpoints_by_key


def _colors(entities: pd.DataFrame) -> list[object]:
    families = sorted(entities["family"].unique())
    palette = {family: plt.get_cmap("tab10")(index) for index, family in enumerate(families)}
    return [palette[family] for family in entities["family"]]


def _save(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=170)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)


def _save_unit_circle(entities: pd.DataFrame, theta: np.ndarray, path: Path, base: float) -> None:
    fig, ax = plt.subplots(figsize=(7, 7))
    circle = np.linspace(0.0, TWO_PI, 400)
    ax.plot(np.cos(circle), np.sin(circle), color="0.55", linewidth=1.0)
    ax.scatter(np.cos(theta), np.sin(theta), c=_colors(entities), s=36)
    for name, x, y in zip(entities["name"], np.cos(theta), np.sin(theta)):
        ax.annotate(name, (x, y), xytext=(3, 3), textcoords="offset points", fontsize=6)
    ax.set(aspect="equal", xlabel="cos(theta)", ylabel="sin(theta)",
           title=f"Persistence phase unit circle: base {base:.6f}")
    _save(fig, path)


def _save_3d(
    entities: pd.DataFrame, theta: np.ndarray, path: Path, base: float, radial: bool
) -> None:
    log_n = np.log10(entities["compton_cycles"].to_numpy(dtype=float))
    normalized = (log_n - log_n.min()) / (np.ptp(log_n) or 1.0)
    radius = normalized if radial else np.ones_like(normalized)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(projection="3d")
    order = np.argsort(log_n)
    ax.plot(x[order], y[order], log_n[order], color="#2457a6", alpha=0.65, linewidth=1.0)
    ax.scatter(x, y, log_n, c=_colors(entities), s=30, depthshade=True)
    ax.set(xlabel="x", ylabel="y", zlabel="log10(N)",
           title=f"{'Radial spiral' if radial else 'Persistence helix'}: base {base:.6f}")
    _save(fig, path)


def _save_gap_polygons(
    theta: np.ndarray,
    path: Path,
    base: float,
    rank: int,
    midpoints_by_key: dict[tuple[int, int], np.ndarray],
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    circle = np.linspace(0.0, TWO_PI, 400)
    for ax, top_k in zip(axes.flat, (3, 4, 5, 6)):
        midpoint = midpoints_by_key[(rank, top_k)]
        x = np.cos(midpoint)
        y = np.sin(midpoint)
        ax.plot(np.cos(circle), np.sin(circle), color="0.75")
        ax.scatter(np.cos(theta), np.sin(theta), color="0.55", alpha=0.25, s=16)
        ax.plot(np.r_[x, x[0]], np.r_[y, y[0]], color="#c83e4d")
        ax.scatter(x, y, color="#c83e4d", s=32)
        ax.set(aspect="equal", title=f"top {top_k} gap centers")
    fig.suptitle(f"Gap-center polygons: base {base:.6f}")
    _save(fig, path)


def _write_video(animation: FuncAnimation, stem: Path, fps: int = 12) -> None:
    animation.save(stem.with_suffix(".gif"), writer=PillowWriter(fps=fps), dpi=95)
    try:
        import imageio_ffmpeg

        plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
        animation.save(stem.with_suffix(".mp4"), writer="ffmpeg", fps=fps, dpi=95)
    except (ImportError, RuntimeError, FileNotFoundError):
        stem.with_suffix(".mp4.unavailable.txt").write_text(
            "MP4 export requires imageio-ffmpeg. Install requirements.txt and rerun.\n",
            encoding="ascii",
        )


def _save_rotation_animation(entities: pd.DataFrame, base: float, theta: np.ndarray) -> None:
    log_n = np.log10(entities["compton_cycles"].to_numpy(dtype=float))
    radius = (log_n - log_n.min()) / (np.ptp(log_n) or 1.0)
    x, y = radius * np.cos(theta), radius * np.sin(theta)
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection="3d")
    order = np.argsort(log_n)
    ax.plot(x[order], y[order], log_n[order], color="#2457a6", alpha=0.7)
    ax.scatter(x, y, log_n, c=_colors(entities), s=30)
    ax.set(xlabel="x", ylabel="y", zlabel="log10(N)",
           title=f"Rotating radial persistence spiral: base {base:.6f}")

    def update(frame: int):
        ax.view_init(elev=24, azim=frame * 3)
        return []

    animation = FuncAnimation(fig, update, frames=120, interval=80)
    _write_video(animation, OUTPUT / "rotating_3d_view")
    plt.close(fig)


def _save_base_sweep_animation(entities: pd.DataFrame) -> None:
    bases = np.linspace(BASE_MIN, BASE_MAX, 120)
    phases = phase_matrix(entities, bases)
    log_n = np.log10(entities["compton_cycles"].to_numpy(dtype=float))
    radius = (log_n - log_n.min()) / (np.ptp(log_n) or 1.0)
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection="3d")
    theta = phases[0]
    scatter = ax.scatter(radius * np.cos(theta), radius * np.sin(theta), log_n,
                         c=_colors(entities), s=30)
    title = ax.set_title("")
    ax.set(xlabel="x", ylabel="y", zlabel="log10(N)")

    def update(frame: int):
        theta = phases[frame]
        scatter._offsets3d = (radius * np.cos(theta), radius * np.sin(theta), log_n)
        title.set_text(f"Radial persistence spiral sweep: base {bases[frame]:.4f}")
        return [scatter, title]

    animation = FuncAnimation(fig, update, frames=len(bases), interval=80)
    _write_video(animation, OUTPUT / "base_sweep_3d")
    plt.close(fig)


def generate_visualizations(
    entities: pd.DataFrame,
    selected: pd.DataFrame,
    midpoints_by_key: dict[tuple[int, int], np.ndarray],
) -> None:
    """Export Wolfram-style and MATLAB-like static and animated figures."""

    OUTPUT.mkdir(parents=True, exist_ok=True)
    SELECTED_OUTPUT.mkdir(parents=True, exist_ok=True)
    phases = phase_matrix(entities, selected["base"].to_numpy(dtype=float))
    for (_, peak), theta in zip(selected.iterrows(), phases):
        rank = int(peak["rank"])
        directory = SELECTED_OUTPUT / f"rank_{rank:02d}_base_{peak['base']:.6f}"
        directory.mkdir(parents=True, exist_ok=True)
        _save_unit_circle(entities, theta, directory / "unit_circle", peak["base"])
        _save_3d(entities, theta, directory / "persistence_helix", peak["base"], radial=False)
        _save_3d(entities, theta, directory / "radial_spiral", peak["base"], radial=True)
        _save_gap_polygons(theta, directory / "gap_center_polygons", peak["base"], rank,
                           midpoints_by_key)
    _save_rotation_animation(entities, selected.iloc[0]["base"], phases[0])
    _save_base_sweep_animation(entities)


def _save_fft_plot(spectrum: pd.DataFrame) -> None:
    frame = spectrum.iloc[1:].copy()
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(frame["frequency_cycles_per_base"], frame["normalized_power"], linewidth=0.9)
    ax.set(xlabel="frequency (cycles per base unit)", ylabel="normalized power",
           title="FFT power spectrum of top-4 symmetry score")
    ax.set_xlim(0.0, min(25.0, frame["frequency_cycles_per_base"].max()))
    fig.tight_layout()
    fig.savefig(OUTPUT / "fft_spectrum.png", dpi=180)
    fig.savefig(OUTPUT / "fft_spectrum.svg")
    plt.close(fig)


def _report(
    peaks: pd.DataFrame,
    matches: pd.DataFrame,
    spectrum: pd.DataFrame,
    geometry: pd.DataFrame,
    relationships: pd.DataFrame,
) -> str:
    closest = matches[matches["match_rank_for_peak"] == 1].head(10)
    dominant = spectrum.iloc[1:].nlargest(5, "power")
    pair_rows = relationships[relationships["record_type"] == "pair"].nlargest(
        8, "persistence_neighbor_score"
    )
    top_geometry = geometry[geometry["top_k"] == 4].nlargest(5, "symmetry_score")
    lines = [
        "# Harmonic Base and 3D Persistence Geometry",
        "",
        "## Scope",
        "",
        "This report compares local maxima in the continuous top-4 gap-center symmetry score with "
        "simple mathematical expressions. Close matches are expected when many peaks and candidate "
        "expressions are tested. They are exploratory numerical coincidences unless supported by a "
        "predeclared null model and independent validation. No physical significance is claimed.",
        "",
        "## Answers",
        "",
        "- **Are peaks random?** This scan alone cannot establish randomness. It detects many local "
        "maxima and a structured FFT spectrum, but both can arise from deterministic phase wrapping.",
        "- **Are peaks related to pi?** Some leading peaks are numerically close to multiples of `pi`, "
        "including the peak near `10*pi`. The candidate search tests many alternatives, so proximity "
        "alone is not evidence of a privileged constant.",
        "- **Are peaks related to e?** Some peaks are numerically close to expressions involving `e`, "
        "including the peak near `e^3`. This remains an exploratory coincidence.",
        "- **Are peaks related to phi?** Some peaks are numerically close to multiples of `phi`. This "
        "also remains exploratory and is not evidence of a privileged constant.",
        "- **Do peaks correspond to geometric organization?** Yes by construction: the ranked signal "
        "is the top-4 gap-center polygon symmetry score. The 3D views show the same phase alignments "
        "with log10(N) added as height.",
        f"- **Which bases generate the strongest 3D structures?** The leading sampled peak is "
        f"`{peaks.iloc[0]['base']:.8f}` with score `{peaks.iloc[0]['score']:.6f}` and prominence "
        f"`{peaks.iloc[0]['prominence']:.6f}`. The table below lists additional candidates.",
        "- **Which relationships persist across base changes?** Stable neighbor pairs are listed below. "
        "They are threshold-dependent descriptive relationships.",
        "",
        "## Top Peaks",
        "",
        "| Rank | Base | Score | Prominence | Width |",
        "| ---: | ---: | ---: | ---: | ---: |",
    ]
    for _, row in peaks.head(20).iterrows():
        lines.append(
            f"| {int(row['rank'])} | {row['base']:.8f} | {row['score']:.6f} | "
            f"{row['prominence']:.6f} | {row['width']:.6f} |"
        )
    lines += ["", "## Closest Harmonic Matches", "",
              "| Peak rank | Peak base | Expression | n | Candidate | Relative error |",
              "| ---: | ---: | --- | ---: | ---: | ---: |"]
    for _, row in closest.iterrows():
        lines.append(
            f"| {int(row['peak_rank'])} | {row['peak_base']:.8f} | `{row['expression']}` | "
            f"{int(row['n'])} | {row['constant_value']:.8f} | {row['relative_error']:.6g} |"
        )
    lines += ["", "## Dominant FFT Components", "",
              "| Frequency (cycles/base) | Period (base units) | Normalized power |",
              "| ---: | ---: | ---: |"]
    for _, row in dominant.iterrows():
        lines.append(
            f"| {row['frequency_cycles_per_base']:.6f} | {row['period_base_units']:.6f} | "
            f"{row['normalized_power']:.6f} |"
        )
    lines += ["", "## Strongest Top-4 Polygon Geometry", "",
              "| Peak rank | Base | Area | Perimeter | Symmetry | Compactness |",
              "| ---: | ---: | ---: | ---: | ---: | ---: |"]
    for _, row in top_geometry.iterrows():
        lines.append(
            f"| {int(row['peak_rank'])} | {row['base']:.8f} | {row['area']:.6f} | "
            f"{row['perimeter']:.6f} | {row['symmetry_score']:.6f} | {row['compactness']:.6f} |"
        )
    lines += ["", "## Persistent Neighbor Pairs", "", "| Pair | Fraction of bases |",
              "| --- | ---: |"]
    for _, row in pair_rows.iterrows():
        lines.append(
            f"| {', '.join(json.loads(row['members']))} | {row['persistence_neighbor_score']:.6f} |"
        )
    lines += [
        "",
        "## Interpretation Limit",
        "",
        "The score signal is generated by logarithmic scaling followed by circular wrapping. Harmonic "
        "matches and FFT components may therefore reflect the transform and finite sample. They do not "
        "demonstrate a resonance mechanism or new physics.",
        "",
        "## Reproduction",
        "",
        "```powershell",
        "python -m src.harmonic_base_geometry",
        "```",
    ]
    return "\n".join(lines) + "\n"


def generate() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Generate peak tables, harmonic matches, geometry, visualizations, and report."""

    scan = pd.read_csv(DATA / "base_scan_metrics.csv")
    entities = load_entities()
    peaks, minima = detect_peaks(scan)
    selected = peaks.head(TOP_BASES).copy()
    matches = constant_matches(peaks)
    spectrum = fft_spectrum(scan)
    geometry, midpoints = gap_geometry_records(entities, selected)
    bases = np.linspace(BASE_MIN, BASE_MAX, BASE_SAMPLES)
    relationships, stability = persistent_relationships(entities, bases, phase_matrix(entities, bases))
    DATA.mkdir(parents=True, exist_ok=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    peaks.to_csv(DATA / "base_peaks.csv", index=False)
    minima.to_csv(DATA / "base_minima.csv", index=False)
    matches.to_csv(DATA / "peak_constant_matches.csv", index=False)
    spectrum.to_csv(DATA / "fft_power_spectrum.csv", index=False)
    geometry.to_csv(DATA / "harmonic_gap_geometry.csv", index=False)
    relationships.to_csv(DATA / "harmonic_persistent_neighbors.csv", index=False)
    stability.to_csv(DATA / "harmonic_particle_stability.csv", index=False)
    _save_fft_plot(spectrum)
    generate_visualizations(entities, selected, midpoints)
    (ROOT / "docs" / "harmonic_base_geometry.md").write_text(
        _report(peaks, matches, spectrum, geometry, relationships), encoding="ascii"
    )
    return peaks, matches, spectrum, geometry


if __name__ == "__main__":
    generated = generate()
    print(
        f"Generated {len(generated[0])} peaks, {len(generated[1])} constant matches, "
        f"and {len(generated[3])} geometry rows in {OUTPUT}"
    )
