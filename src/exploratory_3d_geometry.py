"""Exploratory 3D persistence geometry and randomized controls.

The helix-like views in this module are coordinate constructions. They are
useful for inspecting occupancy and relationships, but they do not imply new
physics or validate a hypothesis.
"""

from __future__ import annotations

import html
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import pandas as pd

from src.dynamic_phase_structure import (
    NEIGHBOR_THRESHOLD_DEG,
    ROOT,
    TWO_PI,
    circular_gaps,
    persistent_relationships,
    phase_matrix,
)
from src.gap_geometry import polygon_metrics
from src.stable_particle_sensitivity import load_variant, phases_for_base

OUTPUT = ROOT / "outputs" / "exploratory_3d"
DATA = ROOT / "data"
REFERENCE_BASE = 10.0
SCAN_BASES = np.linspace(1.2, 50.0, 240)
RANDOM_SEED = 20260602


def normalized(values: np.ndarray) -> np.ndarray:
    """Normalize values to the closed unit interval."""

    values = np.asarray(values, dtype=float)
    span = np.ptp(values)
    return (values - values.min()) / (span or 1.0)


def geometry_coordinates(entities: pd.DataFrame, base: float = REFERENCE_BASE) -> pd.DataFrame:
    """Return all requested coordinate systems for one base."""

    theta = phases_for_base(entities, base)
    log_n = entities["log10_N"].to_numpy(dtype=float)
    radius = normalized(log_n)
    spherical_phi = radius * math.pi
    major_radius, minor_radius = 1.0, 0.38
    frame = entities[["name", "family", "dominant_interaction", "compton_frequency_hz",
                      "lifetime_s", "compton_cycles", "log10_N"]].copy()
    frame["base"] = base
    frame["theta_rad"] = theta
    frame["cylindrical_x"] = np.cos(theta)
    frame["cylindrical_y"] = np.sin(theta)
    frame["cylindrical_z"] = log_n
    frame["radial_x"] = radius * np.cos(theta)
    frame["radial_y"] = radius * np.sin(theta)
    frame["radial_z"] = log_n
    frame["conical_x"] = log_n * np.cos(theta)
    frame["conical_y"] = log_n * np.sin(theta)
    frame["conical_z"] = log_n
    frame["toroidal_x"] = (major_radius + minor_radius * np.cos(TWO_PI * radius)) * np.cos(theta)
    frame["toroidal_y"] = (major_radius + minor_radius * np.cos(TWO_PI * radius)) * np.sin(theta)
    frame["toroidal_z"] = minor_radius * np.sin(TWO_PI * radius)
    frame["spherical_x"] = np.sin(spherical_phi) * np.cos(theta)
    frame["spherical_y"] = np.sin(spherical_phi) * np.sin(theta)
    frame["spherical_z"] = np.cos(spherical_phi)
    frame["landscape_x"] = np.log10(frame["compton_frequency_hz"])
    frame["landscape_y"] = np.log10(frame["lifetime_s"])
    frame["landscape_z"] = theta
    return frame


def randomized_entities(entities: pd.DataFrame, mode: str) -> pd.DataFrame:
    """Return a deterministic shuffled-mass or shuffled-lifetime control."""

    random = np.random.default_rng(RANDOM_SEED + (0 if mode == "lifetime" else 1))
    frame = entities.copy()
    if mode == "lifetime":
        frame["lifetime_s"] = random.permutation(frame["lifetime_s"].to_numpy())
    elif mode == "mass":
        frame["compton_frequency_hz"] = random.permutation(
            frame["compton_frequency_hz"].to_numpy()
        )
    else:
        raise ValueError("mode must be 'lifetime' or 'mass'")
    frame["compton_cycles"] = frame["compton_frequency_hz"] * frame["lifetime_s"]
    frame["log10_N"] = np.log10(frame["compton_cycles"])
    frame["ln_N"] = np.log(frame["compton_cycles"])
    return frame


def gap_constellations(entities: pd.DataFrame, bases: np.ndarray = SCAN_BASES) -> pd.DataFrame:
    """Compute gap-center polygon stability over base space."""

    rows = []
    for base in bases:
        theta = phases_for_base(entities, float(base))
        ordered, _, gaps = circular_gaps(theta)
        for top_k in (3, 4, 5, 6):
            selected = np.argsort(gaps)[::-1][:top_k]
            midpoints = np.mod(ordered[selected] + gaps[selected] / 2.0, TWO_PI)
            metrics = polygon_metrics(np.degrees(midpoints))
            perimeter = float(np.sum(metrics["side_lengths"]))
            rows.append(
                {
                    "base": base,
                    "top_k": top_k,
                    "area": metrics["polygon_area"],
                    "perimeter": perimeter,
                    "symmetry_score": metrics["circular_symmetry_score"],
                    "compactness": (
                        4.0 * math.pi * metrics["polygon_area"] / perimeter**2
                        if perimeter else 0.0
                    ),
                    "midpoint_angles_deg": json.dumps(np.degrees(np.sort(midpoints)).tolist()),
                }
            )
    return pd.DataFrame(rows)


def density_cloud(entities: pd.DataFrame, bases: np.ndarray = SCAN_BASES) -> pd.DataFrame:
    """Accumulate radial persistence positions visited during the base scan."""

    log_n = entities["log10_N"].to_numpy(dtype=float)
    radius = normalized(log_n)
    phases = phase_matrix(entities, bases)
    rows = []
    for base, theta in zip(bases, phases):
        for index, particle in entities.iterrows():
            rows.append(
                {
                    "base": base,
                    "name": particle["name"],
                    "family": particle["family"],
                    "x": radius[index] * math.cos(theta[index]),
                    "y": radius[index] * math.sin(theta[index]),
                    "z": log_n[index],
                }
            )
    return pd.DataFrame(rows)


def family_tubes(entities: pd.DataFrame, bases: np.ndarray = SCAN_BASES) -> pd.DataFrame:
    """Track circular family centroids through base space."""

    phases = phase_matrix(entities, bases)
    rows = []
    for base, theta in zip(bases, phases):
        for family, indices in entities.groupby("family").groups.items():
            selected = theta[list(indices)]
            mean = np.mean(np.exp(1j * selected))
            rows.append(
                {
                    "base": base,
                    "family": family,
                    "theta_rad": float(np.mod(np.angle(mean), TWO_PI)),
                    "x": float(np.real(mean)),
                    "y": float(np.imag(mean)),
                    "resultant_length": abs(mean),
                }
            )
    return pd.DataFrame(rows)


def _colors(entities: pd.DataFrame) -> list[object]:
    labels = sorted(entities["family"].unique())
    palette = {label: plt.get_cmap("tab10")(index) for index, label in enumerate(labels)}
    return [palette[label] for label in entities["family"]]


def _save(fig: plt.Figure, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(OUTPUT / f"{stem}.png", dpi=175)
    fig.savefig(OUTPUT / f"{stem}.svg")
    plt.close(fig)


def _scatter_3d(ax, frame: pd.DataFrame, prefix: str, title: str) -> None:
    order = np.argsort(frame[f"{prefix}_z"])
    ax.plot(frame[f"{prefix}_x"].to_numpy()[order], frame[f"{prefix}_y"].to_numpy()[order],
            frame[f"{prefix}_z"].to_numpy()[order], color="#2457a6", alpha=0.65, linewidth=0.9)
    ax.scatter(frame[f"{prefix}_x"], frame[f"{prefix}_y"], frame[f"{prefix}_z"],
               c=_colors(frame), s=28)
    ax.set(xlabel="x", ylabel="y", zlabel="z", title=title)


def _save_core_views(frame: pd.DataFrame) -> None:
    views = [
        ("cylindrical", "Cylindrical helix"),
        ("radial", "Radial helix"),
        ("conical", "Conical helix"),
        ("toroidal", "Toroidal projection"),
        ("spherical", "Spherical projection"),
        ("landscape", "Persistence landscape"),
    ]
    for prefix, title in views:
        fig = plt.figure(figsize=(8, 7))
        ax = fig.add_subplot(projection="3d")
        _scatter_3d(ax, frame, prefix, f"{title}: base {REFERENCE_BASE:g}")
        _save(fig, prefix)


def _save_density(cloud: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(projection="3d")
    sample = cloud.iloc[::3]
    ax.scatter(sample["x"], sample["y"], sample["z"], c=sample["z"], cmap="viridis",
               s=4, alpha=0.16)
    ax.set(xlabel="x", ylabel="y", zlabel="log10(N)", title="Radial density cloud over base scan")
    _save(fig, "density_cloud")


def _save_worldlines(entities: pd.DataFrame, cloud: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(projection="3d")
    for name, group in cloud.groupby("name"):
        ax.plot(group["base"], group["x"], group["y"], linewidth=0.6, alpha=0.65)
    ax.set(xlabel="base", ylabel="radial x", zlabel="radial y",
           title="Particle worldlines through base space")
    _save(fig, "particle_worldlines")


def _save_constellations(constellations: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    for top_k, group in constellations.groupby("top_k"):
        ax.plot(group["base"], group["symmetry_score"], label=f"top {top_k}", linewidth=0.9)
    ax.set(xlabel="base", ylabel="symmetry score", title="Gap-center constellation stability")
    ax.legend()
    _save(fig, "gap_center_constellations")


def _save_gap_polygon_panel(entities: pd.DataFrame) -> None:
    theta = phases_for_base(entities, REFERENCE_BASE)
    ordered, _, gaps = circular_gaps(theta)
    circle = np.linspace(0.0, TWO_PI, 400)
    fig, axes = plt.subplots(2, 2, figsize=(9, 9))
    for ax, top_k in zip(axes.flat, (3, 4, 5, 6)):
        selected = np.argsort(gaps)[::-1][:top_k]
        midpoint = np.mod(ordered[selected] + gaps[selected] / 2.0, TWO_PI)
        x, y = np.cos(midpoint), np.sin(midpoint)
        ax.plot(np.cos(circle), np.sin(circle), color="0.75")
        ax.scatter(np.cos(theta), np.sin(theta), color="0.55", alpha=0.28, s=16)
        ax.plot(np.r_[x, x[0]], np.r_[y, y[0]], color="#c83e4d")
        ax.scatter(x, y, color="#c83e4d", s=30)
        ax.set(aspect="equal", title=f"top {top_k} gap centers")
        ax.axis("off")
    fig.suptitle(f"Connected gap-center constellations: base {REFERENCE_BASE:g}")
    _save(fig, "gap_center_polygon_panel")


def _save_family_tubes(tubes: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(projection="3d")
    for family, group in tubes.groupby("family"):
        ax.plot(group["base"], group["x"], group["y"], label=family, linewidth=1.3)
    ax.set(xlabel="base", ylabel="centroid x", zlabel="centroid y",
           title="Family centroid tubes through base space")
    ax.legend(fontsize=8)
    _save(fig, "family_tubes")


def _save_neighbor_worldlines(entities: pd.DataFrame, relationships: pd.DataFrame) -> None:
    phases = phase_matrix(entities, SCAN_BASES)
    names = entities["name"].tolist()
    pairs = relationships[relationships["record_type"] == "pair"].nlargest(
        8, "persistence_neighbor_score"
    )
    fig, ax = plt.subplots(figsize=(10, 5))
    for _, row in pairs.iterrows():
        left, right = json.loads(row["members"])
        separation = np.abs(
            (np.degrees(phases[:, names.index(left)] - phases[:, names.index(right)]) + 180.0)
            % 360.0 - 180.0
        )
        ax.plot(SCAN_BASES, separation, label=f"{left} - {right}", linewidth=0.9)
    ax.axhline(NEIGHBOR_THRESHOLD_DEG, color="0.45", linestyle="--", linewidth=0.8)
    ax.set(xlabel="base", ylabel="angular separation (degrees)",
           title="Persistent neighbor worldlines")
    ax.legend(fontsize=7, ncol=2)
    _save(fig, "neighbor_worldlines")


def _save_randomized_controls(entities: pd.DataFrame) -> pd.DataFrame:
    variants = {
        "observed": entities,
        "shuffled_lifetimes": randomized_entities(entities, "lifetime"),
        "shuffled_masses": randomized_entities(entities, "mass"),
    }
    rows = []
    fig = plt.figure(figsize=(15, 15))
    for row, (label, variant) in enumerate(variants.items()):
        frame = geometry_coordinates(variant)
        for column, prefix in enumerate(("cylindrical", "radial", "conical")):
            ax = fig.add_subplot(3, 3, row * 3 + column + 1, projection="3d")
            _scatter_3d(ax, frame, prefix, f"{label.replace('_', ' ')}: {prefix}")
        ordered = frame.sort_values("log10_N")
        turns = float((ordered["theta_rad"].iloc[-1] - ordered["theta_rad"].iloc[0]) / TWO_PI)
        rows.append({"variant": label, "particle_count": len(frame), "endpoint_turn_difference": turns})
    fig.suptitle("Randomized controls: wrapped-logarithm geometry remains coordinate-induced")
    _save(fig, "randomized_controls")
    return pd.DataFrame(rows)


def _write_video(animation: FuncAnimation, stem: str) -> None:
    animation.save(OUTPUT / f"{stem}.gif", writer=PillowWriter(fps=12), dpi=95)
    try:
        import imageio_ffmpeg

        plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
        animation.save(OUTPUT / f"{stem}.mp4", writer="ffmpeg", fps=12, dpi=95)
    except (ImportError, RuntimeError, FileNotFoundError):
        (OUTPUT / f"{stem}.mp4.unavailable.txt").write_text(
            "MP4 export requires imageio-ffmpeg. Install requirements.txt and rerun.\n",
            encoding="ascii",
        )


def _save_rotation(frame: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection="3d")
    _scatter_3d(ax, frame, "radial", "Rotating radial persistence helix")

    def update(index: int):
        ax.view_init(elev=24, azim=index * 3)
        return []

    animation = FuncAnimation(fig, update, frames=120, interval=80)
    _write_video(animation, "rotating_radial_helix")
    plt.close(fig)


def _save_base_sweep(entities: pd.DataFrame) -> None:
    bases = np.linspace(1.2, 50.0, 120)
    phases = phase_matrix(entities, bases)
    log_n = entities["log10_N"].to_numpy(dtype=float)
    radius = normalized(log_n)
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection="3d")
    points = ax.scatter(radius * np.cos(phases[0]), radius * np.sin(phases[0]), log_n,
                        c=_colors(entities), s=28)
    title = ax.set_title("")
    ax.set(xlabel="x", ylabel="y", zlabel="log10(N)")

    def update(index: int):
        points._offsets3d = (
            radius * np.cos(phases[index]), radius * np.sin(phases[index]), log_n
        )
        title.set_text(f"Base sweep radial helix: base {bases[index]:.4f}")
        return [points, title]

    animation = FuncAnimation(fig, update, frames=len(bases), interval=80)
    _write_video(animation, "base_sweep_radial_helix")
    plt.close(fig)


def _save_interactive_html(entities: pd.DataFrame) -> None:
    bases = np.linspace(1.2, 50.0, 120)
    phases = phase_matrix(entities, bases)
    log_n = entities["log10_N"].to_numpy(dtype=float)
    radius = normalized(log_n)
    payload = {
        "bases": bases.tolist(),
        "names": entities["name"].tolist(),
        "z": log_n.tolist(),
        "r": radius.tolist(),
        "theta": phases.tolist(),
    }
    document = f"""<!doctype html>
<meta charset="utf-8"><title>Interactive radial persistence geometry</title>
<style>body{{font-family:Arial;margin:20px}}svg{{border:1px solid #ccc}}label{{display:block}}</style>
<h1>Interactive radial persistence geometry</h1>
<p>Exploratory finite-lifetime-only view. Drag the base slider and rotation slider.</p>
<label>Base <input id="base" type="range" min="0" max="119" value="0"> <span id="baseValue"></span></label>
<label>Rotation <input id="rotation" type="range" min="0" max="360" value="25"></label>
<svg id="plot" width="760" height="620"></svg>
<script>
const data={json.dumps(payload)}; const svg=document.getElementById("plot");
const ns="http://www.w3.org/2000/svg"; const base=document.getElementById("base");
const rotation=document.getElementById("rotation"); const label=document.getElementById("baseValue");
function draw(){{svg.innerHTML=""; const i=+base.value, a=+rotation.value*Math.PI/180; label.textContent=data.bases[i].toFixed(4);
data.names.forEach((name,j)=>{{const t=data.theta[i][j]+a, x=data.r[j]*Math.cos(t), y=data.r[j]*Math.sin(t), z=(data.z[j]-Math.min(...data.z))/(Math.max(...data.z)-Math.min(...data.z));
const cx=380+210*x, cy=560-330*z-70*y; const c=document.createElementNS(ns,"circle"); c.setAttribute("cx",cx);c.setAttribute("cy",cy);c.setAttribute("r",4);c.setAttribute("fill","#2457a6");
const title=document.createElementNS(ns,"title");title.textContent=name;c.appendChild(title);svg.appendChild(c);}});}}
base.oninput=draw;rotation.oninput=draw;draw();
</script>"""
    (OUTPUT / "interactive_radial_geometry.html").write_text(document, encoding="ascii")


def _report(relationships: pd.DataFrame, controls: pd.DataFrame) -> str:
    pairs = relationships[relationships["record_type"] == "pair"].nlargest(
        8, "persistence_neighbor_score"
    )
    lines = [
        "# Exploratory 3D Persistence Geometry",
        "",
        "## Scope",
        "",
        "This study searches for unexpected geometry without validating a hypothesis. It uses the "
        "`finite_lifetime_only` dataset so stable-particle truncation does not drive the result. "
        "All structures are exploratory visual geometry only. No new physics is claimed.",
        "",
        "## Why a Helix Appears",
        "",
        "A cylindrical helix is induced by the coordinates. At base `10`, `z = log10(N)` and "
        "`theta = 2*pi*frac(log10(N))`, so increasing `z` necessarily winds around the unit circle. "
        "Other bases rescale the winding rate. Randomized controls retain a helix-like trace because "
        "the same wrapping rule is applied after shuffling. The helix itself is therefore primarily a "
        "visualization, logarithmic, and circular-projection artifact.",
        "",
        "## What Geometries Appear?",
        "",
        "- Cylindrical, radial, and conical views show winding layers because phase is a wrapped logarithm.",
        "- Toroidal and spherical projections reveal occupancy bands and voids, but their locations move "
        "with the chosen projection and base.",
        "- The density cloud shows frequently visited radial regions as base changes.",
        "- Particle worldlines and family tubes orbit and cross through base space.",
        "- Gap-center constellation symmetry changes with base rather than remaining fixed.",
        "",
        "## What Disappears Under Controls?",
        "",
        "Shuffling lifetimes or masses changes point ordering, local clusters, gaps, and neighbor identities. "
        "The broad helix-like winding survives because it is encoded by the transform. Dataset-specific "
        "occupancy patterns do not survive unchanged.",
        "",
        "## Which Structures Are Robust?",
        "",
        "The robust structure is the wrapped-logarithm winding itself. Local gaps, polygons, bands, and "
        "family trajectories are conditional on the dataset, base, and projection. Persistent neighbor "
        "pairs are the most reproducible dataset-specific relationships in this scan.",
        "",
        "## Persistent Neighbor Pairs",
        "",
        "| Pair | Fraction of bases |",
        "| --- | ---: |",
    ]
    for _, row in pairs.iterrows():
        lines.append(f"| {', '.join(json.loads(row['members']))} | {row['persistence_neighbor_score']:.6f} |")
    lines += [
        "",
        "## Randomized Controls",
        "",
        "| Variant | Particle count | Endpoint turn difference |",
        "| --- | ---: | ---: |",
    ]
    for _, row in controls.iterrows():
        lines.append(
            f"| {row['variant']} | {int(row['particle_count'])} | {row['endpoint_turn_difference']:.6f} |"
        )
    lines += [
        "",
        "## Interpretation Limit",
        "",
        "This is an exploratory geometry study. It does not establish a physical helix, preferred base, "
        "resonance, or new physics.",
        "",
        "## Reproduction",
        "",
        "```powershell",
        "python -m src.exploratory_3d_geometry",
        "```",
    ]
    return "\n".join(lines) + "\n"


def generate() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Generate exploratory geometry tables, figures, controls, animations, and report."""

    OUTPUT.mkdir(parents=True, exist_ok=True)
    entities = load_variant()
    frame = geometry_coordinates(entities)
    cloud = density_cloud(entities)
    constellations = gap_constellations(entities)
    tubes = family_tubes(entities)
    relationships, stability = persistent_relationships(
        entities, SCAN_BASES, phase_matrix(entities, SCAN_BASES)
    )
    controls = _save_randomized_controls(entities)
    frame.to_csv(DATA / "exploratory_geometry_coordinates.csv", index=False)
    cloud.to_csv(DATA / "exploratory_density_cloud.csv", index=False)
    constellations.to_csv(DATA / "exploratory_gap_constellations.csv", index=False)
    tubes.to_csv(DATA / "exploratory_family_tubes.csv", index=False)
    relationships.to_csv(DATA / "exploratory_neighbor_worldlines.csv", index=False)
    stability.to_csv(DATA / "exploratory_particle_stability.csv", index=False)
    controls.to_csv(DATA / "exploratory_randomized_controls.csv", index=False)
    _save_core_views(frame)
    _save_density(cloud)
    _save_worldlines(entities, cloud)
    _save_constellations(constellations)
    _save_gap_polygon_panel(entities)
    _save_family_tubes(tubes)
    _save_neighbor_worldlines(entities, relationships)
    _save_rotation(frame)
    _save_base_sweep(entities)
    _save_interactive_html(entities)
    (ROOT / "docs" / "exploratory_geometry.md").write_text(
        _report(relationships, controls), encoding="ascii"
    )
    return frame, cloud, constellations


if __name__ == "__main__":
    generated = generate()
    print(
        f"Generated {len(generated[0])} reference coordinates, {len(generated[1])} cloud rows, "
        f"and {len(generated[2])} constellation rows in {OUTPUT}"
    )
