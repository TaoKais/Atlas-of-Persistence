from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

# Exact SI values for c and h. G is the CODATA 2022 recommended value.
C = 299_792_458.0
H = 6.626_070_15e-34
HBAR = H / (2.0 * math.pi)
G = 6.674_30e-11
EV_J = 1.602_176_634e-19
MEV_J = 1.0e6 * EV_J

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUT = ROOT / "output"


@dataclass(frozen=True)
class Particle:
    name: str
    family: str
    mass_mev: float
    lifetime_s: float | None
    stability: str
    decay_channel_count: int
    source: str

    @property
    def rest_energy_j(self) -> float:
        return self.mass_mev * MEV_J

    @property
    def compton_frequency_hz(self) -> float:
        return self.rest_energy_j / H

    @property
    def compton_cycles(self) -> float | None:
        return (
            self.compton_frequency_hz * self.lifetime_s
            if self.lifetime_s is not None
            else None
        )

    @property
    def quality_factor(self) -> float | None:
        return (
            2.0 * math.pi * self.compton_cycles
            if self.compton_cycles is not None
            else None
        )

    @property
    def relative_width(self) -> float | None:
        return (
            1.0 / self.quality_factor
            if self.quality_factor not in (None, 0.0)
            else None
        )


def read_csv(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_particles() -> list[Particle]:
    particles = []
    for row in read_csv("particles.csv"):
        particles.append(
            Particle(
                name=row["name"],
                family=row["family"],
                mass_mev=float(row["mass_mev"]),
                lifetime_s=float(row["lifetime_s"]) if row["lifetime_s"] else None,
                stability=row["stability"],
                decay_channel_count=int(row["decay_channel_count"]),
                source=row["source"],
            )
        )
    return particles


def compactness(mass_kg: float, radius_m: float) -> float:
    return G * mass_kg / (radius_m * C**2)


def rho_energy(total_energy_j: float, rest_energy_j: float) -> float:
    return total_energy_j / rest_energy_j


def destroyed_exergy(entropy_generated_j_per_k: float, ambient_temperature_k: float) -> float:
    return entropy_generated_j_per_k * ambient_temperature_k


def exergy_destruction_fraction(
    input_exergy_j: float, entropy_generated_j_per_k: float, ambient_temperature_k: float
) -> float:
    return destroyed_exergy(entropy_generated_j_per_k, ambient_temperature_k) / input_exergy_j


def logarithmic_gaps(
    particles: Iterable[Particle],
) -> list[dict[str, float | str]]:
    ordered = sorted(particles, key=lambda particle: particle.compton_frequency_hz)
    gaps = []
    for lower, upper in zip(ordered, ordered[1:]):
        frequency_ratio = upper.compton_frequency_hz / lower.compton_frequency_hz
        mass_ratio = upper.mass_mev / lower.mass_mev
        gaps.append(
            {
                "lower": lower.name,
                "upper": upper.name,
                "frequency_ratio": frequency_ratio,
                "mass_ratio": mass_ratio,
                "log10_gap": math.log10(frequency_ratio),
            }
        )
    return sorted(gaps, key=lambda gap: float(gap["log10_gap"]), reverse=True)


def ranks(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    result = [0.0] * len(values)
    index = 0
    while index < len(indexed):
        end = index + 1
        while end < len(indexed) and indexed[end][1] == indexed[index][1]:
            end += 1
        average_rank = (index + 1 + end) / 2.0
        for original_index, _ in indexed[index:end]:
            result[original_index] = average_rank
        index = end
    return result


def spearman(values_a: list[float], values_b: list[float]) -> float:
    ranked_a = ranks(values_a)
    ranked_b = ranks(values_b)
    mean_a = sum(ranked_a) / len(ranked_a)
    mean_b = sum(ranked_b) / len(ranked_b)
    numerator = sum((a - mean_a) * (b - mean_b) for a, b in zip(ranked_a, ranked_b))
    denominator = math.sqrt(
        sum((a - mean_a) ** 2 for a in ranked_a)
        * sum((b - mean_b) ** 2 for b in ranked_b)
    )
    return numerator / denominator if denominator else math.nan


def hypothesis_metric_rows(particles: list[Particle]) -> list[dict[str, object]]:
    unstable = [particle for particle in particles if particle.lifetime_s is not None]
    log_lifetime = [math.log10(particle.lifetime_s or 0.0) for particle in unstable]
    candidates = {
        "log10_mass_mev": [math.log10(particle.mass_mev) for particle in unstable],
        "log10_compton_frequency_hz": [
            math.log10(particle.compton_frequency_hz) for particle in unstable
        ],
        "decay_channel_count": [float(particle.decay_channel_count) for particle in unstable],
        "log10_compton_cycles": [
            math.log10(particle.compton_cycles or 0.0) for particle in unstable
        ],
    }
    return [
        {
            "candidate": name,
            "target": "log10_lifetime_s",
            "sample_size": len(unstable),
            "spearman_rank_correlation": spearman(values, log_lifetime),
            "interpretation": (
                "descriptive_only; requires a larger catalogue and controlled validation"
            ),
        }
        for name, values in candidates.items()
    ]


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def particle_rows(particles: list[Particle]) -> list[dict[str, object]]:
    return [
        {
            "name": particle.name,
            "family": particle.family,
            "mass_mev": particle.mass_mev,
            "compton_frequency_hz": particle.compton_frequency_hz,
            "lifetime_s": particle.lifetime_s if particle.lifetime_s is not None else "",
            "stability": particle.stability,
            "decay_channel_count": particle.decay_channel_count,
            "compton_cycles": particle.compton_cycles if particle.compton_cycles is not None else "",
            "quality_factor": particle.quality_factor if particle.quality_factor is not None else "",
            "relative_width": particle.relative_width if particle.relative_width is not None else "",
        }
        for particle in particles
    ]


def compact_object_rows() -> list[dict[str, object]]:
    rows = []
    for row in read_csv("compact_objects.csv"):
        phi = compactness(float(row["mass_kg"]), float(row["radius_m"]))
        rows.append({**row, "phi": phi, "two_phi": 2.0 * phi})
    return rows


def exergy_rows() -> list[dict[str, object]]:
    rows = []
    for row in read_csv("exergy_scenarios.csv"):
        input_exergy = float(row["input_exergy_j"])
        entropy_generated = float(row["entropy_generated_j_per_k"])
        ambient_temperature = float(row["ambient_temperature_k"])
        rows.append(
            {
                **row,
                "destroyed_exergy_j": destroyed_exergy(entropy_generated, ambient_temperature),
                "delta_x": exergy_destruction_fraction(
                    input_exergy, entropy_generated, ambient_temperature
                ),
            }
        )
    return rows


def summary_markdown(
    particles: list[Particle],
    gaps: list[dict[str, float | str]],
    objects: list[dict[str, object]],
) -> str:
    unstable = [particle for particle in particles if particle.lifetime_s is not None]
    shortest = min(unstable, key=lambda particle: particle.lifetime_s or math.inf)
    longest = max(unstable, key=lambda particle: particle.lifetime_s or -math.inf)
    largest_gap = gaps[0]
    horizon = next(row for row in objects if row["name"] == "schwarzschild_horizon_reference")
    return f"""# Generated exploratory summary

## Dataset

- Particles: {len(particles)}
- Stable identities represented without an assigned decay lifetime: {len(particles) - len(unstable)}
- Unstable identities: {len(unstable)}

## Immediate checks

- `f_C = mc^2/h` is exactly proportional to mass; its isolated ranking cannot
  improve on mass ranking.
- Largest sampled logarithmic gap: `{largest_gap["lower"]}` to
  `{largest_gap["upper"]}`, ratio `{float(largest_gap["frequency_ratio"]):.6g}`.
- Shortest sampled lifetime: `{shortest.name}`, `{shortest.lifetime_s:.6g} s`.
- Longest finite sampled lifetime: `{longest.name}`, `{longest.lifetime_s:.6g} s`.
- Schwarzschild horizon reference: `Phi = {float(horizon["phi"]):.6g}`;
  equivalently `2 Phi = {float(horizon["two_phi"]):.6g}`.

## Interpretation limits

- Gaps depend on catalogue selection and observational bias.
- Stable particles need lower bounds, not invented finite lifetimes.
- `decay_channel_count` is an exploratory annotation, not a phase-space volume.
- The exergy scenarios are formula demonstrations, not fitted experiments.
"""


def generate() -> None:
    OUTPUT.mkdir(exist_ok=True)
    particles = load_particles()
    gaps = logarithmic_gaps(particles)
    objects = compact_object_rows()
    write_csv(OUTPUT / "particles_analysis.csv", particle_rows(particles))
    write_csv(OUTPUT / "frequency_gaps.csv", gaps)
    write_csv(OUTPUT / "hypothesis_metrics.csv", hypothesis_metric_rows(particles))
    write_csv(OUTPUT / "compactness_analysis.csv", objects)
    write_csv(OUTPUT / "exergy_analysis.csv", exergy_rows())
    (OUTPUT / "summary.md").write_text(
        summary_markdown(particles, gaps, objects), encoding="utf-8"
    )


if __name__ == "__main__":
    generate()
    print(f"Generated reports in {OUTPUT}")
