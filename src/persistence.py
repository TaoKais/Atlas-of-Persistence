"""Measured-input persistence landscape calculations.

The generated results are descriptive. They do not establish a new physical law.
"""

from __future__ import annotations

import csv
import html
import math
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from statistics import fmean, pstdev
from typing import Callable, Iterable

H = 6.626_070_15e-34
EV_J = 1.602_176_634e-19
MEV_J = 1.0e6 * EV_J

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUT = ROOT / "output" / "persistence"


class FrequencyType(str, Enum):
    """Supported characteristic-frequency definitions."""

    COMPTON = "compton"
    SCHUMANN = "schumann"
    ORBITAL = "orbital"
    PLASMA = "plasma"
    CUSTOM = "custom"


@dataclass(frozen=True)
class CharacteristicFrequency:
    value_hz: float
    frequency_type: FrequencyType
    definition: str

    def __post_init__(self) -> None:
        _require_positive("value_hz", self.value_hz)


@dataclass(frozen=True)
class PersistenceRecord:
    particle: str
    mass_mev: float
    lifetime_s: float
    lifetime_qualifier: str
    lifetime_basis: str
    source: str

    @property
    def frequency_hz(self) -> float:
        return compton_frequency(self.mass_mev)

    @property
    def persistence_index(self) -> float:
        return persistence_cycles(self.frequency_hz, self.lifetime_s)

    @property
    def log10_persistence(self) -> float:
        return persistence_log10(self.persistence_index)


def _require_positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")


def compton_frequency(mass_mev: float) -> float:
    """Return f_C = E / h for a rest energy expressed in MeV."""

    _require_positive("mass_mev", mass_mev)
    return mass_mev * MEV_J / H


def persistence_cycles(frequency_hz: float, lifetime_s: float) -> float:
    """Return the dimensionless persistence index N = f * tau."""

    _require_positive("frequency_hz", frequency_hz)
    _require_positive("lifetime_s", lifetime_s)
    return frequency_hz * lifetime_s


def persistence_log10(n: float) -> float:
    """Return log10(N) for a positive persistence index."""

    _require_positive("N", n)
    return math.log10(n)


def characteristic_frequency(
    value_hz: float, frequency_type: FrequencyType, definition: str
) -> CharacteristicFrequency:
    """Create a traceable measured characteristic frequency for any system."""

    return CharacteristicFrequency(value_hz, frequency_type, definition)


def load_records() -> list[PersistenceRecord]:
    with (DATA / "persistence_particles.csv").open(newline="", encoding="utf-8") as handle:
        return [
            PersistenceRecord(
                particle=row["particle"],
                mass_mev=float(row["mass_mev"]),
                lifetime_s=float(row["lifetime_s"]),
                lifetime_qualifier=row["lifetime_qualifier"],
                lifetime_basis=row["lifetime_basis"],
                source=row["source"],
            )
            for row in csv.DictReader(handle)
        ]


def table_rows(records: Iterable[PersistenceRecord]) -> list[dict[str, object]]:
    return [
        {
            "Particle": record.particle,
            "Mass_MeV": record.mass_mev,
            "Lifetime_s": record.lifetime_s,
            "Lifetime_qualifier": record.lifetime_qualifier,
            "Frequency_Hz": record.frequency_hz,
            "N": record.persistence_index,
            "log10_N": record.log10_persistence,
            "Lifetime_basis": record.lifetime_basis,
            "Source": record.source,
        }
        for record in records
    ]


def gap_rows(records: Iterable[PersistenceRecord]) -> list[dict[str, object]]:
    ordered = sorted(records, key=lambda record: record.log10_persistence)
    gaps = [
        {
            "Lower_particle": lower.particle,
            "Upper_particle": upper.particle,
            "Lower_log10_N": lower.log10_persistence,
            "Upper_log10_N": upper.log10_persistence,
            "Delta_log10_N": upper.log10_persistence - lower.log10_persistence,
        }
        for lower, upper in zip(ordered, ordered[1:])
    ]
    for rank, gap in enumerate(
        sorted(gaps, key=lambda row: float(row["Delta_log10_N"]), reverse=True), start=1
    ):
        gap["Largest_gap_rank"] = rank
    return sorted(gaps, key=lambda row: int(row["Largest_gap_rank"]))


def cluster_rows(records: Iterable[PersistenceRecord]) -> tuple[list[dict[str, object]], float]:
    """Split the ordered sample at unusually large gaps for visual exploration."""

    ordered = sorted(records, key=lambda record: record.log10_persistence)
    adjacent = [
        upper.log10_persistence - lower.log10_persistence
        for lower, upper in zip(ordered, ordered[1:])
    ]
    threshold = fmean(adjacent) + pstdev(adjacent)
    cluster = 1
    rows = []
    for index, record in enumerate(ordered):
        if index and adjacent[index - 1] > threshold:
            cluster += 1
        rows.append(
            {
                "Particle": record.particle,
                "log10_N": record.log10_persistence,
                "Cluster": cluster,
                "Method": "adjacent gap > mean(gaps) + population_stddev(gaps)",
                "Threshold": threshold,
            }
        )
    return rows, threshold


def landscape_metric_rows(records: Iterable[PersistenceRecord]) -> list[dict[str, object]]:
    """Summarize each landscape without treating lower bounds as exact observations."""

    records = list(records)
    landscapes = {
        "mass_mev": [(record.particle, math.log10(record.mass_mev)) for record in records],
        "lifetime_s": [(record.particle, math.log10(record.lifetime_s)) for record in records],
        "persistence_N": [(record.particle, record.log10_persistence) for record in records],
    }
    rows = []
    for landscape, values in landscapes.items():
        ordered = sorted(values, key=lambda item: item[1])
        gaps = [
            (upper[1] - lower[1], lower[0], upper[0])
            for lower, upper in zip(ordered, ordered[1:])
        ]
        largest_gap, lower, upper = max(gaps)
        rows.append(
            {
                "Landscape": landscape,
                "Log10_range": ordered[-1][1] - ordered[0][1],
                "Largest_adjacent_log10_gap": largest_gap,
                "Largest_gap_lower_particle": lower,
                "Largest_gap_upper_particle": upper,
                "Interpretation": "descriptive_only; electron and proton lifetime-derived values are lower bounds",
            }
        )
    return rows


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def _svg(
    path: Path,
    title: str,
    items: list[tuple[str, float, float]],
    x_label: str,
    y_label: str,
    *,
    x_log: bool = False,
    y_log: bool = False,
) -> None:
    width, height = 840, 520
    left, right, top, bottom = 90, 30, 55, 75
    inner_width, inner_height = width - left - right, height - top - bottom
    transform_x: Callable[[float], float] = math.log10 if x_log else lambda value: value
    transform_y: Callable[[float], float] = math.log10 if y_log else lambda value: value
    points = [(name, transform_x(x), transform_y(y)) for name, x, y in items]
    x_values = [point[1] for point in points]
    y_values = [point[2] for point in points]
    x_min, x_max = min(x_values), max(x_values)
    y_min, y_max = min(y_values), max(y_values)
    x_pad = (x_max - x_min or 1.0) * 0.06
    y_pad = (y_max - y_min or 1.0) * 0.10
    x_min, x_max = x_min - x_pad, x_max + x_pad
    y_min, y_max = y_min - y_pad, y_max + y_pad

    def sx(value: float) -> float:
        return left + (value - x_min) * inner_width / (x_max - x_min)

    def sy(value: float) -> float:
        return top + inner_height - (value - y_min) * inner_height / (y_max - y_min)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="28" text-anchor="middle" font-size="18">{html.escape(title)}</text>',
        f'<line x1="{left}" y1="{top + inner_height}" x2="{left + inner_width}" y2="{top + inner_height}" stroke="black"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + inner_height}" stroke="black"/>',
        f'<text x="{width / 2}" y="{height - 18}" text-anchor="middle">{html.escape(x_label)}</text>',
        f'<text x="22" y="{height / 2}" transform="rotate(-90 22 {height / 2})" text-anchor="middle">{html.escape(y_label)}</text>',
    ]
    for index in range(6):
        x_value = x_min + index * (x_max - x_min) / 5
        y_value = y_min + index * (y_max - y_min) / 5
        parts.extend(
            [
                f'<text x="{sx(x_value):.1f}" y="{top + inner_height + 20}" text-anchor="middle" font-size="11">{x_value:.2f}</text>',
                f'<text x="{left - 9}" y="{sy(y_value) + 4:.1f}" text-anchor="end" font-size="11">{y_value:.2f}</text>',
            ]
        )
    for name, x_value, y_value in points:
        x, y = sx(x_value), sy(y_value)
        parts.extend(
            [
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="#2457a6"/>',
                f'<text x="{x + 7:.1f}" y="{y - 7:.1f}" font-size="11">{html.escape(name)}</text>',
            ]
        )
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def _bar_svg(path: Path, title: str, labels: list[str], values: list[float], y_label: str) -> None:
    width, height = 840, 520
    left, right, top, bottom = 90, 30, 55, 115
    inner_width, inner_height = width - left - right, height - top - bottom
    low, high = min(0.0, min(values)), max(values)
    span = high - low or 1.0
    bar_width = inner_width / len(values) * 0.68

    def sy(value: float) -> float:
        return top + inner_height - (value - low) * inner_height / span

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="28" text-anchor="middle" font-size="18">{html.escape(title)}</text>',
        f'<line x1="{left}" y1="{sy(low):.1f}" x2="{left + inner_width}" y2="{sy(low):.1f}" stroke="black"/>',
        f'<text x="22" y="{height / 2}" transform="rotate(-90 22 {height / 2})" text-anchor="middle">{html.escape(y_label)}</text>',
    ]
    for index, (label, value) in enumerate(zip(labels, values)):
        x = left + (index + 0.5) * inner_width / len(values)
        y = sy(value)
        parts.extend(
            [
                f'<rect x="{x - bar_width / 2:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{sy(low) - y:.1f}" fill="#2457a6"/>',
                f'<text x="{x:.1f}" y="{sy(low) + 18:.1f}" transform="rotate(45 {x:.1f} {sy(low) + 18:.1f})" font-size="11">{html.escape(label)}</text>',
            ]
        )
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def _gap_histogram_svg(path: Path, gaps: list[dict[str, object]]) -> None:
    values = [float(gap["Delta_log10_N"]) for gap in gaps]
    bins = 6
    high = max(values)
    width = high / bins
    counts = [0] * bins
    for value in values:
        counts[min(int(value / width), bins - 1)] += 1
    labels = [f"{index * width:.1f}-{(index + 1) * width:.1f}" for index in range(bins)]
    _bar_svg(path, "Adjacent persistence-gap histogram", labels, counts, "count")


def _summary(
    records: list[PersistenceRecord],
    gaps: list[dict[str, object]],
    metrics: list[dict[str, object]],
    threshold: float,
) -> str:
    ranked = "\n".join(
        f"| {gap['Largest_gap_rank']} | {gap['Lower_particle']} | {gap['Upper_particle']} | {float(gap['Delta_log10_N']):.6g} |"
        for gap in gaps
    )
    comparisons = "\n".join(
        f"| {row['Landscape']} | {float(row['Log10_range']):.6g} | {float(row['Largest_adjacent_log10_gap']):.6g} | {row['Largest_gap_lower_particle']} | {row['Largest_gap_upper_particle']} |"
        for row in metrics
    )
    return f"""# Persistence landscape numerical results

This report is descriptive and exploratory. It does not claim discovery of a
new law, natural classification, or stability island.

## Sample

- Focused measured-input particle records: {len(records)}
- Exact or central lifetime values: {sum(record.lifetime_qualifier == "central_value" for record in records)}
- Experimental lower bounds: {sum(record.lifetime_qualifier == "lower_bound" for record in records)}
- Clustering split threshold: `{threshold:.6g}` log10 units

Electron and proton persistence values are lower bounds. The proton input is a
mode-specific partial-lifetime bound and is not a measured total proton
lifetime.

## Largest adjacent gaps

| Rank | Lower particle | Upper particle | Delta log10(N) |
| --- | --- | --- | ---: |
{ranked}

## Landscape comparison

| Landscape | log10 range | Largest adjacent log10 gap | Lower particle | Upper particle |
| --- | ---: | ---: | --- | --- |
{comparisons}

## Interpretation limit

The sample is small and curated. The gaps and clusters can change when the
catalogue or experimental bounds change. The generated plots support inspection
of the research question; they do not by themselves answer it.
"""


def generate() -> None:
    records = load_records()
    rows = table_rows(records)
    gaps = gap_rows(records)
    clusters, threshold = cluster_rows(records)
    metrics = landscape_metric_rows(records)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    _write_csv(DATA / "persistence_table.csv", rows)
    _write_csv(OUTPUT / "persistence_gaps.csv", gaps)
    _write_csv(OUTPUT / "persistence_clusters.csv", clusters)
    _write_csv(OUTPUT / "landscape_metrics.csv", metrics)
    ordered = sorted(records, key=lambda record: record.log10_persistence)
    _bar_svg(
        OUTPUT / "a_log10_n_sorted.svg",
        "Persistence index sorted ascending",
        [record.particle for record in ordered],
        [record.log10_persistence for record in ordered],
        "log10(N)",
    )
    _svg(
        OUTPUT / "b_frequency_vs_lifetime.svg",
        "Frequency vs lifetime",
        [(record.particle, record.frequency_hz, record.lifetime_s) for record in records],
        "log10(frequency / Hz)",
        "log10(lifetime / s)",
        x_log=True,
        y_log=True,
    )
    _svg(
        OUTPUT / "c_frequency_vs_persistence.svg",
        "Frequency vs persistence index",
        [(record.particle, record.frequency_hz, record.persistence_index) for record in records],
        "log10(frequency / Hz)",
        "log10(N)",
        x_log=True,
        y_log=True,
    )
    _gap_histogram_svg(OUTPUT / "d_gap_histogram.svg", gaps)
    _svg(
        OUTPUT / "e_cluster_landscape.svg",
        "Exploratory adjacent-gap clusters",
        [(str(row["Particle"]), float(row["log10_N"]), float(row["Cluster"])) for row in clusters],
        "log10(N)",
        "cluster identifier",
    )
    (OUTPUT / "summary.md").write_text(
        _summary(records, gaps, metrics, threshold), encoding="utf-8"
    )


if __name__ == "__main__":
    generate()
    print(f"Generated persistence reports in {OUTPUT}")
