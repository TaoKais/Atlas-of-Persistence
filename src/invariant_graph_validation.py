"""Invariant graph validation with explicit tolerances.

The routines in this module compare descriptive particle/entity graphs against
randomized controls. They are exploratory validation tools and do not establish
new physics or physical interaction.
"""

from __future__ import annotations

import json
import math
from collections import Counter, deque
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.persistence import compton_frequency

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CONFIG = ROOT / "config"
DOCS = ROOT / "docs"
OUTPUT = ROOT / "outputs" / "invariant_graph_validation"
TWO_PI = 2.0 * math.pi
PHI = (1.0 + math.sqrt(5.0)) / 2.0
FIXED_BASES = {
    "2": 2.0,
    "e": math.e,
    "pi": math.pi,
    "phi": PHI,
    "10": 10.0,
    "2*pi": TWO_PI,
    "sqrt(2*pi)": math.sqrt(TWO_PI),
    "pi^2": math.pi**2,
    "e^pi": math.e**math.pi,
}
CONTROL_TYPES = (
    "shuffle_lifetimes",
    "shuffle_masses",
    "shuffle_N_values",
    "random_logN_same_range",
    "random_phase_same_z",
    "family_label_shuffle",
    "interaction_label_shuffle",
)
REPRESENTATIONS = (
    "unit_circle",
    "cylindrical_helicoid",
    "radial_helicoid",
    "physical_plane",
    "physical_3d",
)


@dataclass(frozen=True)
class GraphBuild:
    adjacency: np.ndarray
    weights: np.ndarray
    mean_distance: np.ndarray
    min_distance: np.ndarray
    max_distance: np.ndarray
    connected_counts: np.ndarray
    base_count: int
    per_base_connected: np.ndarray


def load_tolerances(path: Path = CONFIG / "tolerances.json") -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def base_grid(samples: int = 2000) -> tuple[list[str], np.ndarray]:
    labels = list(FIXED_BASES)
    values = list(FIXED_BASES.values())
    for value in np.linspace(1.2, 50.0, samples):
        labels.append(f"{value:.6f}")
        values.append(float(value))
    return labels, np.asarray(values, dtype=float)


def _stable_tau(mode: str) -> float | None:
    if mode in {"finite_lifetime_only", "stable_excluded", "stable_as_infinity_marker"}:
        return None
    if mode == "stable_truncated_1e25":
        return 1e25
    if mode == "stable_truncated_1e30":
        return 1e30
    if mode == "stable_truncated_1e35":
        return 1e35
    raise ValueError(f"Unknown stable handling mode: {mode}")


def load_entities(
    path: Path = DATA / "particles.csv",
    stable_mode: str = "finite_lifetime_only",
) -> pd.DataFrame:
    particles = pd.read_csv(path).copy()
    particles["mass_mev"] = pd.to_numeric(particles["mass_mev"], errors="coerce")
    particles["lifetime_s"] = pd.to_numeric(particles["lifetime_s"], errors="coerce")
    stable = particles["stability"].fillna("").str.lower().eq("stable")
    stable_tau = _stable_tau(stable_mode)
    if stable_tau is None:
        particles = particles[~stable & particles["lifetime_s"].gt(0.0)].copy()
        particles["stable_handling"] = stable_mode
    else:
        particles.loc[stable, "lifetime_s"] = stable_tau
        particles = particles[particles["lifetime_s"].gt(0.0)].copy()
        particles["stable_handling"] = np.where(stable.loc[particles.index], stable_mode, "finite_value")
    particles = particles[particles["mass_mev"].gt(0.0) & particles["lifetime_s"].gt(0.0)].copy()
    particles["compton_frequency_hz"] = particles["mass_mev"].map(compton_frequency)
    particles["compton_cycles"] = particles.get(
        "compton_cycles", particles["compton_frequency_hz"] * particles["lifetime_s"]
    )
    particles["N"] = particles["compton_cycles"]
    particles["logN"] = np.log10(particles["N"])
    return particles.reset_index(drop=True)


def angular_distance_deg(left: np.ndarray | float, right: np.ndarray | float) -> np.ndarray:
    return np.abs((np.asarray(left) - np.asarray(right) + 180.0) % 360.0 - 180.0)


def phases_deg(entities: pd.DataFrame, bases: np.ndarray) -> np.ndarray:
    log_n = entities["logN"].to_numpy(dtype=float)
    return 360.0 * np.mod(log_n[None, :] / np.log10(bases)[:, None], 1.0)


def representation_coordinates(entities: pd.DataFrame, base: float, representation: str) -> np.ndarray:
    log_n = entities["logN"].to_numpy(dtype=float)
    theta = TWO_PI * np.mod(log_n / math.log10(base), 1.0)
    unit = np.column_stack([np.cos(theta), np.sin(theta)])
    if representation == "unit_circle":
        return unit
    if representation == "cylindrical_helicoid":
        return np.column_stack([unit, log_n])
    if representation == "radial_helicoid":
        span = float(np.ptp(log_n)) or 1.0
        radius = (log_n - float(log_n.min())) / span
        return np.column_stack([radius * np.cos(theta), radius * np.sin(theta), log_n])
    if representation == "physical_plane":
        return np.column_stack([
            np.log10(entities["compton_frequency_hz"].to_numpy(dtype=float)),
            np.log10(entities["lifetime_s"].to_numpy(dtype=float)),
        ])
    if representation == "physical_3d":
        return np.column_stack([
            np.log10(entities["compton_frequency_hz"].to_numpy(dtype=float)),
            np.log10(entities["lifetime_s"].to_numpy(dtype=float)),
            log_n,
        ])
    raise ValueError(f"Unknown representation: {representation}")


def pairwise_distances(coords: np.ndarray) -> np.ndarray:
    delta = coords[:, None, :] - coords[None, :, :]
    return np.sqrt(np.sum(delta * delta, axis=2))


def angular_graph(entities: pd.DataFrame, bases: np.ndarray, tolerance_deg: float) -> GraphBuild:
    theta = phases_deg(entities, bases)
    distances = angular_distance_deg(theta[:, :, None], theta[:, None, :])
    connected = distances < tolerance_deg
    return _build_from_base_connections(connected, distances)


def distance_graph(
    entities: pd.DataFrame, bases: np.ndarray, percentile: float, representation: str = "physical_3d"
) -> GraphBuild:
    if representation in {"physical_plane", "physical_3d"}:
        distances = pairwise_distances(representation_coordinates(entities, float(bases[0]), representation))
        upper = distances[np.triu_indices_from(distances, k=1)]
        threshold = float(np.percentile(upper, percentile)) if len(upper) else 0.0
        connected = distances < threshold
        return _build_from_base_connections(
            np.repeat(connected[None, :, :], len(bases), axis=0),
            np.repeat(distances[None, :, :], len(bases), axis=0),
        )
    distance_layers = []
    connected_layers = []
    for base in bases:
        distances = pairwise_distances(representation_coordinates(entities, float(base), representation))
        upper = distances[np.triu_indices_from(distances, k=1)]
        threshold = float(np.percentile(upper, percentile)) if len(upper) else 0.0
        distance_layers.append(distances)
        connected_layers.append(distances < threshold)
    return _build_from_base_connections(np.asarray(connected_layers), np.asarray(distance_layers))


def knn_graph(entities: pd.DataFrame, bases: np.ndarray, k: int, representation: str = "physical_3d") -> GraphBuild:
    if representation in {"physical_plane", "physical_3d"}:
        distances = pairwise_distances(representation_coordinates(entities, float(bases[0]), representation))
        connected = np.zeros_like(distances, dtype=bool)
        for index in range(len(entities)):
            order = np.argsort(distances[index])
            for neighbor in order[1 : k + 1]:
                connected[index, neighbor] = True
                connected[neighbor, index] = True
        return _build_from_base_connections(
            np.repeat(connected[None, :, :], len(bases), axis=0),
            np.repeat(distances[None, :, :], len(bases), axis=0),
        )
    distance_layers = []
    connected_layers = []
    for base in bases:
        distances = pairwise_distances(representation_coordinates(entities, float(base), representation))
        connected = np.zeros_like(distances, dtype=bool)
        for index in range(len(entities)):
            order = np.argsort(distances[index])
            for neighbor in order[1 : k + 1]:
                connected[index, neighbor] = True
                connected[neighbor, index] = True
        distance_layers.append(distances)
        connected_layers.append(connected)
    return _build_from_base_connections(np.asarray(connected_layers), np.asarray(distance_layers))


def multi_representation_graph(
    entities: pd.DataFrame, bases: np.ndarray, tolerance_deg: float, required_representations: int = 3
) -> GraphBuild:
    layers = []
    distances = []
    for base in bases:
        per_rep = []
        per_dist = []
        for representation in REPRESENTATIONS:
            coords = representation_coordinates(entities, float(base), representation)
            dist = pairwise_distances(coords)
            upper = dist[np.triu_indices_from(dist, k=1)]
            threshold = float(np.percentile(upper, tolerance_deg)) if representation != "unit_circle" else 2.0 * math.sin(math.radians(tolerance_deg) / 2.0)
            per_rep.append(dist < threshold)
            per_dist.append(dist)
        layers.append(np.sum(per_rep, axis=0) >= required_representations)
        distances.append(np.mean(per_dist, axis=0))
    return _build_from_base_connections(np.asarray(layers), np.asarray(distances))


def _build_from_base_connections(connected: np.ndarray, distances: np.ndarray) -> GraphBuild:
    connected = connected.copy()
    for layer in connected:
        np.fill_diagonal(layer, False)
    base_count = int(connected.shape[0])
    counts = connected.sum(axis=0)
    weights = counts / max(base_count, 1)
    mean_distance = distances.mean(axis=0)
    min_distance = distances.min(axis=0)
    max_distance = distances.max(axis=0)
    adjacency = weights > 0.0
    np.fill_diagonal(adjacency, False)
    return GraphBuild(adjacency, weights, mean_distance, min_distance, max_distance, counts, base_count, connected)


def thresholded_adjacency(build: GraphBuild, threshold: float) -> np.ndarray:
    adjacency = build.weights >= threshold
    np.fill_diagonal(adjacency, False)
    return adjacency


def graph_metrics(adjacency: np.ndarray, weights: np.ndarray, labels: Iterable[str]) -> dict[str, float]:
    labels = list(labels)
    n = adjacency.shape[0]
    edges = int(np.triu(adjacency, 1).sum())
    degrees = adjacency.sum(axis=1).astype(float)
    components = connected_components(adjacency)
    largest = max((len(component) for component in components), default=0)
    clustering = clustering_coefficient(adjacency)
    shortest = average_shortest_path(adjacency)
    spectral = spectral_gap(adjacency)
    return {
        "number_of_edges": float(edges),
        "average_edge_weight": float(weights[np.triu(adjacency, 1)].mean()) if edges else 0.0,
        "modularity_by_label": modularity(adjacency, labels),
        "clustering_coefficient": clustering,
        "average_shortest_path": shortest,
        "connected_components": float(len(components)),
        "degree_distribution": float(np.std(degrees)),
        "spectral_gap": spectral,
        "largest_component_size": float(largest),
    }


def connected_components(adjacency: np.ndarray) -> list[list[int]]:
    seen: set[int] = set()
    components: list[list[int]] = []
    for start in range(adjacency.shape[0]):
        if start in seen:
            continue
        queue = deque([start])
        seen.add(start)
        component = []
        while queue:
            node = queue.popleft()
            component.append(node)
            for neighbor in np.flatnonzero(adjacency[node]):
                if int(neighbor) not in seen:
                    seen.add(int(neighbor))
                    queue.append(int(neighbor))
        components.append(component)
    return components


def clustering_coefficient(adjacency: np.ndarray) -> float:
    values = []
    for node in range(adjacency.shape[0]):
        neighbors = np.flatnonzero(adjacency[node])
        if len(neighbors) < 2:
            values.append(0.0)
            continue
        sub_edges = adjacency[np.ix_(neighbors, neighbors)].sum() / 2.0
        values.append(float(sub_edges / (len(neighbors) * (len(neighbors) - 1) / 2.0)))
    return float(np.mean(values)) if values else 0.0


def average_shortest_path(adjacency: np.ndarray) -> float:
    distances = []
    for component in connected_components(adjacency):
        if len(component) < 2:
            continue
        sub = adjacency[np.ix_(component, component)]
        for start in range(len(component)):
            dist = np.full(len(component), np.inf)
            dist[start] = 0.0
            queue = deque([start])
            while queue:
                node = queue.popleft()
                for neighbor in np.flatnonzero(sub[node]):
                    if not np.isfinite(dist[neighbor]):
                        dist[neighbor] = dist[node] + 1.0
                        queue.append(int(neighbor))
            distances.extend(dist[np.isfinite(dist) & (dist > 0.0)].tolist())
    return float(np.mean(distances)) if distances else 0.0


def spectral_gap(adjacency: np.ndarray) -> float:
    if adjacency.shape[0] < 2:
        return 0.0
    degrees = np.diag(adjacency.sum(axis=1))
    laplacian = degrees - adjacency.astype(float)
    values = np.sort(np.linalg.eigvalsh(laplacian))
    return float(values[1]) if len(values) > 1 else 0.0


def modularity(adjacency: np.ndarray, labels: list[str]) -> float:
    m = float(adjacency.sum() / 2.0)
    if m == 0.0:
        return 0.0
    degrees = adjacency.sum(axis=1)
    q = 0.0
    for i in range(adjacency.shape[0]):
        for j in range(adjacency.shape[0]):
            if labels[i] == labels[j]:
                q += adjacency[i, j] - degrees[i] * degrees[j] / (2.0 * m)
    return float(q / (2.0 * m))


def empirical_p_value(observed: float, controls: np.ndarray) -> float:
    if len(controls) == 0:
        return 1.0
    return float((np.sum(controls >= observed) + 1.0) / (len(controls) + 1.0))


def apply_control(entities: pd.DataFrame, control_type: str, rng: np.random.Generator) -> pd.DataFrame:
    controlled = entities.copy()
    if control_type == "shuffle_lifetimes":
        controlled["lifetime_s"] = rng.permutation(controlled["lifetime_s"].to_numpy(dtype=float))
        controlled["N"] = controlled["compton_frequency_hz"] * controlled["lifetime_s"]
        controlled["logN"] = np.log10(controlled["N"])
    elif control_type == "shuffle_masses":
        controlled["mass_mev"] = rng.permutation(controlled["mass_mev"].to_numpy(dtype=float))
        controlled["compton_frequency_hz"] = controlled["mass_mev"].map(compton_frequency)
        controlled["N"] = controlled["compton_frequency_hz"] * controlled["lifetime_s"]
        controlled["logN"] = np.log10(controlled["N"])
    elif control_type == "shuffle_N_values":
        controlled["logN"] = rng.permutation(controlled["logN"].to_numpy(dtype=float))
        controlled["N"] = 10.0 ** controlled["logN"]
    elif control_type == "random_logN_same_range":
        low, high = controlled["logN"].min(), controlled["logN"].max()
        controlled["logN"] = rng.uniform(low, high, len(controlled))
        controlled["N"] = 10.0 ** controlled["logN"]
    elif control_type == "random_phase_same_z":
        base = 10.0
        integer = np.floor(controlled["logN"].to_numpy(dtype=float) / math.log10(base))
        controlled["logN"] = (integer + rng.random(len(controlled))) * math.log10(base)
        controlled["N"] = 10.0 ** controlled["logN"]
    elif control_type == "family_label_shuffle":
        controlled["family"] = rng.permutation(controlled["family"].to_numpy())
    elif control_type == "interaction_label_shuffle":
        controlled["dominant_interaction"] = rng.permutation(controlled["dominant_interaction"].to_numpy())
    else:
        raise ValueError(f"Unknown control type: {control_type}")
    return controlled


def edge_table(
    entities: pd.DataFrame,
    build: GraphBuild,
    tolerance_type: str,
    tolerance_value: float,
    edge_p_values: dict[tuple[str, str], float],
    p_threshold: float,
) -> pd.DataFrame:
    rows = []
    names = entities["name"].tolist()
    for i, j in combinations(range(len(names)), 2):
        if build.weights[i, j] <= 0.0:
            continue
        key = tuple(sorted((names[i], names[j])))
        p_value = edge_p_values.get(key, 1.0)
        rows.append(
            {
                "particle_a": names[i],
                "particle_b": names[j],
                "tolerance_type": tolerance_type,
                "tolerance_value": tolerance_value,
                "base_count": build.base_count,
                "connected_count": int(build.connected_counts[i, j]),
                "persistence_fraction": float(build.weights[i, j]),
                "mean_distance": float(build.mean_distance[i, j]),
                "min_distance": float(build.min_distance[i, j]),
                "max_distance": float(build.max_distance[i, j]),
                "family_match": bool(entities.loc[i, "family"] == entities.loc[j, "family"]),
                "interaction_match": bool(entities.loc[i, "dominant_interaction"] == entities.loc[j, "dominant_interaction"]),
                "empirical_p_value": p_value,
                "passes_threshold": bool(p_value < p_threshold),
            }
        )
    return pd.DataFrame(rows)


def validation_summary(
    entities: pd.DataFrame,
    bases: np.ndarray,
    control_bases: np.ndarray,
    graph_type: str,
    tolerance_type: str,
    tolerance_value: float,
    persistence_threshold: float,
    runs: int,
    p_threshold: float,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, dict[tuple[str, str], float], GraphBuild]:
    observed_build = build_graph(entities, bases, graph_type, tolerance_value)
    observed_adj = thresholded_adjacency(observed_build, persistence_threshold)
    rows = []
    edge_counts: Counter[tuple[str, str]] = Counter()
    names = entities["name"].tolist()
    observed_family = graph_metrics(observed_adj, observed_build.weights, entities["family"])
    observed_interaction = graph_metrics(observed_adj, observed_build.weights, entities["dominant_interaction"])
    observed_metrics = {
        **observed_family,
        "modularity_by_family": observed_family["modularity_by_label"],
        "modularity_by_interaction": observed_interaction["modularity_by_label"],
    }
    observed_metrics.pop("modularity_by_label", None)
    for control_type in CONTROL_TYPES:
        control_metric_rows = []
        for _ in range(runs):
            controlled = apply_control(entities, control_type, rng)
            control_build = build_graph(controlled, control_bases, graph_type, tolerance_value)
            control_adj = thresholded_adjacency(control_build, persistence_threshold)
            if control_type not in {"family_label_shuffle", "interaction_label_shuffle"}:
                for i, j in combinations(range(len(names)), 2):
                    if control_build.weights[i, j] >= persistence_threshold:
                        edge_counts[tuple(sorted((names[i], names[j])))] += 1
            family_metrics = graph_metrics(control_adj, control_build.weights, controlled["family"])
            interaction_metrics = graph_metrics(control_adj, control_build.weights, controlled["dominant_interaction"])
            metric_row = {
                **family_metrics,
                "modularity_by_family": family_metrics["modularity_by_label"],
                "modularity_by_interaction": interaction_metrics["modularity_by_label"],
            }
            metric_row.pop("modularity_by_label", None)
            control_metric_rows.append(metric_row)
        control_frame = pd.DataFrame(control_metric_rows)
        for metric_name, observed in observed_metrics.items():
            controls = control_frame[metric_name].to_numpy(dtype=float)
            mean = float(np.mean(controls))
            std = float(np.std(controls))
            p_value = empirical_p_value(float(observed), controls)
            rows.append(
                {
                    "graph_type": graph_type,
                    "tolerance_type": tolerance_type,
                    "tolerance_value": tolerance_value,
                    "base_count": int(len(bases)),
                    "control_base_count": int(len(control_bases)),
                    "monte_carlo_runs": int(runs),
                    "control_type": control_type,
                    "observed_metric": metric_name,
                    "observed_value": float(observed),
                    "control_mean": mean,
                    "control_std": std,
                    "z_score": float((observed - mean) / std) if std > 0.0 else 0.0,
                    "empirical_p_value": p_value,
                    "passes_p05": bool(p_value < p_threshold),
                    "interpretation": interpretation_label(p_value, observed, mean, tolerance_value),
                }
            )
    denom = max(1, runs * 5)
    edge_p = {key: (count + 1.0) / (denom + 1.0) for key, count in edge_counts.items()}
    return pd.DataFrame(rows), edge_p, observed_build


def build_graph(entities: pd.DataFrame, bases: np.ndarray, graph_type: str, tolerance_value: float) -> GraphBuild:
    if graph_type == "angular_closeness":
        return angular_graph(entities, bases, tolerance_value)
    if graph_type == "3d_distance":
        return distance_graph(entities, bases, tolerance_value, "physical_3d")
    if graph_type == "knn":
        return knn_graph(entities, bases, int(tolerance_value), "physical_3d")
    if graph_type == "multi_representation":
        return multi_representation_graph(entities, bases, tolerance_value, 3)
    if graph_type == "multi_base_persistent":
        return angular_graph(entities, bases, tolerance_value)
    raise ValueError(f"Unknown graph type: {graph_type}")


def interpretation_label(p_value: float, observed: float, control_mean: float, tolerance_value: float) -> str:
    if not math.isfinite(observed):
        return "Insufficient data"
    if p_value < 0.05 and observed > control_mean:
        return "Robust candidate"
    if p_value >= 0.05 and abs(observed - control_mean) <= max(abs(control_mean) * 0.05, 1e-12):
        return "Likely artifact"
    if tolerance_value in {5, 30, 1, 5}:
        return "Tolerance-sensitive"
    return "Weak candidate"


def robust_clusters(
    entities: pd.DataFrame,
    build: GraphBuild,
    graph_type: str,
    tolerance_type: str,
    tolerance_value: float,
    persistence_threshold: float,
    edge_p_values: dict[tuple[str, str], float],
) -> pd.DataFrame:
    adjacency = thresholded_adjacency(build, persistence_threshold)
    names = entities["name"].tolist()
    rows = []
    for cluster_id, component in enumerate(connected_components(adjacency), start=1):
        if len(component) < 2:
            continue
        sub_weights = build.weights[np.ix_(component, component)]
        upper = sub_weights[np.triu_indices(len(component), 1)]
        member_names = [names[index] for index in component]
        families = entities.loc[component, "family"].tolist()
        interactions = entities.loc[component, "dominant_interaction"].tolist()
        p_values = [
            edge_p_values.get(tuple(sorted((names[i], names[j]))), 1.0)
            for i, j in combinations(component, 2)
            if adjacency[i, j]
        ]
        rows.append(
            {
                "cluster_id": cluster_id,
                "members": ";".join(member_names),
                "graph_type": graph_type,
                "tolerance_type": tolerance_type,
                "tolerance_value": tolerance_value,
                "persistence_fraction": float(np.mean(upper)) if len(upper) else 0.0,
                "family_purity": max(Counter(families).values()) / len(families),
                "interaction_purity": max(Counter(interactions).values()) / len(interactions),
                "empirical_p_value": float(np.mean(p_values)) if p_values else 1.0,
            }
        )
    return pd.DataFrame(rows)


def run_validation(stable_mode: str = "finite_lifetime_only", samples: int = 2000, seed: int = 1729) -> None:
    DATA.mkdir(exist_ok=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(exist_ok=True)
    config = load_tolerances()
    runs = int(config["monte_carlo_runs"])
    p_threshold = float(config["empirical_p_value_threshold"])
    persistence_threshold = float(config["persistence_fraction_threshold"][1])
    _, bases = base_grid(samples)
    control_indices = np.unique(np.linspace(0, len(bases) - 1, 25, dtype=int))
    control_bases = bases[control_indices]
    entities = load_entities(stable_mode=stable_mode)
    rng = np.random.default_rng(seed)

    summary_frames = []
    edge_frames = []
    cluster_frames = []
    graph_specs = [
        ("angular_closeness", "angular_close_deg", float(config["angular_close_deg"][1])),
        ("3d_distance", "euclidean_close_percentile", float(config["euclidean_close_percentile"][1])),
        ("knn", "neighbor_rank_k", float(config["neighbor_rank_k"][1])),
        ("multi_representation", "angular_close_deg", float(config["angular_close_deg"][1])),
        ("multi_base_persistent", "angular_close_deg", float(config["angular_close_deg"][1])),
    ]
    primary_build = None
    primary_edge_p: dict[tuple[str, str], float] = {}
    for graph_type, tolerance_type, tolerance_value in graph_specs:
        summary, edge_p, build = validation_summary(
            entities, bases, control_bases, graph_type, tolerance_type, tolerance_value,
            persistence_threshold, runs, p_threshold, rng
        )
        summary_frames.append(summary)
        edge_frames.append(edge_table(entities, build, tolerance_type, tolerance_value, edge_p, p_threshold))
        cluster_frames.append(
            robust_clusters(entities, build, graph_type, tolerance_type, tolerance_value, persistence_threshold, edge_p)
        )
        if graph_type == "angular_closeness":
            primary_build = build
            primary_edge_p = edge_p

    assert primary_build is not None
    invariant_edges = pd.concat(edge_frames, ignore_index=True)
    summary = pd.concat(summary_frames, ignore_index=True)
    clusters = pd.concat(cluster_frames, ignore_index=True) if cluster_frames else pd.DataFrame()
    invariant_edges.to_csv(DATA / "invariant_edges.csv", index=False)
    summary.to_csv(DATA / "graph_validation_summary.csv", index=False)
    clusters.to_csv(DATA / "robust_clusters.csv", index=False)
    public_data = ROOT / "archipelago-explorer" / "public" / "data"
    if public_data.exists():
        invariant_edges.to_csv(public_data / "invariant_edges.csv", index=False)
        summary.to_csv(public_data / "graph_validation_summary.csv", index=False)
        clusters.to_csv(public_data / "robust_clusters.csv", index=False)
    _save_figures(entities, bases, primary_build, invariant_edges, summary, clusters, primary_edge_p, config)
    _write_report(entities, invariant_edges, summary, clusters, config, len(bases), runs)


def _save_figures(
    entities: pd.DataFrame,
    bases: np.ndarray,
    build: GraphBuild,
    edges: pd.DataFrame,
    summary: pd.DataFrame,
    clusters: pd.DataFrame,
    edge_p: dict[tuple[str, str], float],
    config: dict[str, object],
) -> None:
    tolerance = float(config["angular_close_deg"][1])
    threshold = float(config["persistence_fraction_threshold"][1])
    _plot_network(entities, build, threshold, tolerance, OUTPUT / "01_persistent_neighbor_graph.png")
    _plot_tolerance_sweep(entities, bases, config, OUTPUT / "02_tolerance_sweep_heatmap.png")
    _plot_base_heatmap(entities, bases[:300], tolerance, OUTPUT / "03_base_robustness_heatmap.png")
    _plot_metric_controls(summary, OUTPUT / "04_real_vs_randomized_graph_metric_plots.png")
    _plot_modularity(summary, "modularity_by_family", OUTPUT / "05_modularity_by_family_vs_controls.png")
    _plot_modularity(summary, "modularity_by_interaction", OUTPUT / "06_modularity_by_interaction_vs_controls.png")
    _plot_modularity(summary, "spectral_gap", OUTPUT / "07_spectral_graph_comparison.png")
    _plot_clusters(clusters, OUTPUT / "08_robust_cluster_diagram.png")
    _plot_survival(edges, OUTPUT / "09_edge_survival_curve.png")
    _plot_helicoid_edges(entities, build, threshold, tolerance, OUTPUT / "10_3d_helicoid_with_edges.png")


def _pair_labels(names: list[str]) -> list[str]:
    return [f"{a}/{b}" for a, b in combinations(names, 2)]


def _plot_network(entities: pd.DataFrame, build: GraphBuild, threshold: float, tolerance: float, path: Path) -> None:
    theta = np.linspace(0.0, TWO_PI, len(entities), endpoint=False)
    pos = np.column_stack([np.cos(theta), np.sin(theta)])
    families = {family: index for index, family in enumerate(sorted(entities["family"].unique()))}
    fig, ax = plt.subplots(figsize=(9, 9))
    for i, j in combinations(range(len(entities)), 2):
        if build.weights[i, j] >= threshold:
            ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]], color="0.25", alpha=build.weights[i, j], linewidth=1 + 3 * build.weights[i, j])
            midpoint = (pos[i] + pos[j]) / 2.0
            ax.text(midpoint[0], midpoint[1], f"{tolerance:g}deg/{build.weights[i,j]:.2f}", fontsize=6)
    for index, row in entities.iterrows():
        ax.scatter(pos[index, 0], pos[index, 1], s=90, c=[families[row["family"]]], cmap="tab10", vmin=0, vmax=max(families.values()) or 1)
        ax.text(pos[index, 0] * 1.08, pos[index, 1] * 1.08, row["name"], fontsize=7, ha="center", va="center")
    ax.set_title("Persistent neighbor graph: tolerance 10deg, threshold 0.7")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _plot_tolerance_sweep(entities: pd.DataFrame, bases: np.ndarray, config: dict[str, object], path: Path) -> None:
    names = entities["name"].tolist()
    labels = _pair_labels(names)
    values = []
    for tolerance in config["angular_close_deg"]:
        build = angular_graph(entities, bases, float(tolerance))
        values.append([build.weights[i, j] for i, j in combinations(range(len(names)), 2)])
    matrix = np.asarray(values).T
    fig, ax = plt.subplots(figsize=(8, max(5, len(labels) * 0.18)))
    image = ax.imshow(matrix, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(config["angular_close_deg"])), labels=[str(x) for x in config["angular_close_deg"]])
    ax.set_yticks(range(len(labels)), labels=labels, fontsize=5)
    ax.set_xlabel("angular_close_deg")
    ax.set_title("Tolerance sweep heatmap: persistence_fraction")
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _plot_base_heatmap(entities: pd.DataFrame, bases: np.ndarray, tolerance: float, path: Path) -> None:
    names = entities["name"].tolist()
    build = angular_graph(entities, bases, tolerance)
    matrix = np.asarray([
        build.per_base_connected[:, i, j].astype(int) for i, j in combinations(range(len(names)), 2)
    ])
    fig, ax = plt.subplots(figsize=(10, max(5, matrix.shape[0] * 0.18)))
    image = ax.imshow(matrix, aspect="auto", cmap="Greys", vmin=0, vmax=1)
    ax.set_xlabel("base index")
    ax.set_yticks(range(len(_pair_labels(names))), labels=_pair_labels(names), fontsize=5)
    ax.set_title(f"Base robustness heatmap: angular_close_deg={tolerance:g}")
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _plot_metric_controls(summary: pd.DataFrame, path: Path) -> None:
    frame = summary[summary["graph_type"].eq("angular_closeness")].head(80)
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.errorbar(range(len(frame)), frame["control_mean"], yerr=frame["control_std"], fmt="o", label="controls")
    ax.scatter(range(len(frame)), frame["observed_value"], marker="x", color="black", label="observed")
    ax.set_xticks(range(len(frame)), frame["observed_metric"], rotation=90, fontsize=6)
    ax.set_title("Real vs randomized graph metrics: angular_close_deg=10")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _plot_modularity(summary: pd.DataFrame, metric: str, path: Path) -> None:
    frame = summary[summary["observed_metric"].eq(metric)]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(range(len(frame)), frame["observed_value"], label="observed")
    ax.errorbar(range(len(frame)), frame["control_mean"], yerr=frame["control_std"], fmt=".", color="black", label="control mean")
    ax.set_xticks(range(len(frame)), frame["control_type"], rotation=45, ha="right", fontsize=7)
    ax.set_title(metric.replace("_", " "))
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _plot_clusters(clusters: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    if clusters.empty:
        ax.text(0.5, 0.5, "No robust clusters at selected tolerance", ha="center", va="center")
    else:
        frame = clusters.sort_values("persistence_fraction", ascending=False).head(12)
        ax.barh(frame["members"], frame["persistence_fraction"], color="tab:green")
        ax.set_xlim(0, 1)
    ax.set_title("Robust cluster diagram: threshold 0.7")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _plot_survival(edges: pd.DataFrame, path: Path) -> None:
    values = np.sort(edges["persistence_fraction"].to_numpy(dtype=float))[::-1]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(np.arange(1, len(values) + 1), values)
    ax.axhline(0.7, color="0.4", linestyle="--", label="threshold 0.7")
    ax.set(xlabel="edge rank", ylabel="persistence_fraction", title="Edge survival curve")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _plot_helicoid_edges(entities: pd.DataFrame, build: GraphBuild, threshold: float, tolerance: float, path: Path) -> None:
    coords = representation_coordinates(entities, 10.0, "cylindrical_helicoid")
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    for i, j in combinations(range(len(entities)), 2):
        if build.weights[i, j] >= threshold:
            ax.plot(coords[[i, j], 0], coords[[i, j], 1], coords[[i, j], 2], color="0.2", alpha=build.weights[i, j])
    for index, row in entities.iterrows():
        ax.scatter(coords[index, 0], coords[index, 1], coords[index, 2], s=35)
        ax.text(coords[index, 0], coords[index, 1], coords[index, 2], row["name"], fontsize=6)
    ax.set_title(f"3D helicoid with persistent edges: angular_close_deg={tolerance:g}")
    ax.set_xlabel("cos(theta)")
    ax.set_ylabel("sin(theta)")
    ax.set_zlabel("log10(N)")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _write_report(
    entities: pd.DataFrame,
    edges: pd.DataFrame,
    summary: pd.DataFrame,
    clusters: pd.DataFrame,
    config: dict[str, object],
    base_count: int,
    runs: int,
) -> None:
    top_edges = edges.sort_values(["passes_threshold", "persistence_fraction"], ascending=[False, False]).head(8)
    significant = summary[summary["passes_p05"]]
    sensitive = edges.groupby("tolerance_value")["persistence_fraction"].mean().to_dict()
    text = [
        "# Invariant Graph Validation",
        "",
        "Persistent graph structure does not imply physical interaction. It indicates only that selected entities remain close under the chosen mathematical representations and tolerances.",
        "",
        "All outputs are exploratory unless they pass randomized controls and tolerance checks. No new physics or discovery is claimed.",
        "",
        "## Validation Settings",
        f"- Particles/entities tested: {len(entities)}",
        f"- Bases tested: {base_count}",
        "- Randomized controls use the generated representative control-base grid recorded in `graph_validation_summary.csv`.",
        f"- Monte Carlo runs per control: {runs}",
        f"- Primary angular tolerance: {config['angular_close_deg'][1]} degrees",
        f"- Primary persistence threshold: {config['persistence_fraction_threshold'][1]}",
        f"- Empirical p-value threshold: {config['empirical_p_value_threshold']}",
        "",
        "## 1. Which relationships survive changes of base?",
    ]
    if top_edges.empty:
        text.append("No persistent edges were found at the reported tolerances.")
    else:
        for row in top_edges.itertuples():
            text.append(
                f"- {row.particle_a} / {row.particle_b}: tolerance {row.tolerance_type}={row.tolerance_value}, "
                f"{row.connected_count}/{row.base_count} bases, persistence {row.persistence_fraction:.3f}, empirical p={row.empirical_p_value:.4f}."
            )
    text.extend(
        [
            "",
            "## 2. Which relationships survive changes of representation?",
            "The multi-representation graph requires closeness in at least three representations. Relationships absent from that graph are labeled representation-sensitive rather than robust.",
            "",
            "## 3. Which relationships survive randomized controls?",
            f"{len(significant)} graph metric rows pass p < 0.05. Rows that do not pass should be treated as weak candidates or likely artifacts.",
            "",
            "## 4. Which relationships are tolerance-sensitive?",
            f"Mean persistence by reported tolerance: {json.dumps({str(k): round(float(v), 4) for k, v in sensitive.items()})}. Strong changes across this map indicate tolerance-sensitive structure.",
            "",
            "## 5. Which graph metrics are stronger than random?",
        ]
    )
    for row in significant.head(12).itertuples():
        text.append(
            f"- {row.graph_type}, {row.observed_metric}, tolerance {row.tolerance_type}={row.tolerance_value}: "
            f"observed {row.observed_value:.4g}, null mean {row.control_mean:.4g}, null std {row.control_std:.4g}, empirical p={row.empirical_p_value:.4f}."
        )
    if significant.empty:
        text.append("- No metric is stronger than random at p < 0.05 in the generated summary.")
    text.extend(
        [
            "",
            "## 6. Which clusters are most robust?",
        ]
    )
    if clusters.empty:
        text.append("No multi-node robust clusters were produced at the selected persistence threshold.")
    else:
        for row in clusters.sort_values("persistence_fraction", ascending=False).head(8).itertuples():
            text.append(
                f"- Cluster {row.cluster_id}: {row.members}; tolerance {row.tolerance_type}={row.tolerance_value}; "
                f"persistence {row.persistence_fraction:.3f}; family purity {row.family_purity:.3f}; interaction purity {row.interaction_purity:.3f}; p={row.empirical_p_value:.4f}."
            )
    text.extend(
        [
            "",
            "## 7. Are family or interaction labels predictive?",
            "Use `modularity_by_family` and `modularity_by_interaction` in `data/graph_validation_summary.csv`. A label is treated as predictive only when observed modularity exceeds the matching randomized-label controls with p < 0.05.",
            "",
            "## 8. Is the observed structure mainly driven by N?",
            "Controls that shuffle or randomize N directly test this. If structure disappears under `shuffle_N_values`, `random_logN_same_range`, or `random_phase_same_z`, the reported graph is mainly driven by the chosen N projection.",
            "",
            "## 9. Is there evidence beyond visualization?",
            "Only rows with explicit tolerance, base count, Monte Carlo count, observed metric, null mean, null std, and empirical p-value provide evidence beyond visualization. Non-significant rows are reported as non-significant.",
            "",
            "## 10. What should not be claimed?",
            "Do not claim discovery, new interaction, preferred mathematical base, or physical causation. Persistent closeness is a property of selected representations, tolerances, and data preprocessing.",
        ]
    )
    (DOCS / "invariant_graph_validation.md").write_text("\n".join(text) + "\n", encoding="utf-8")


if __name__ == "__main__":
    run_validation()
