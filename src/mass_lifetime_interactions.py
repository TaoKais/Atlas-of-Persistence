"""Validate mass-lifetime scaling against interaction labels.

The analysis deliberately separates measured lifetime from the persistence
cycle definition N = f_C * tau, because N reuses mass through f_C.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.atlas_persistence import DATA, EV_J, H, HBAR, MEV_J, ROOT, load_particles

CONFIG = ROOT / "config" / "mass_lifetime_validation_tolerances.json"
OUTPUT = ROOT / "outputs" / "mass_lifetime_interactions"
DOC = ROOT / "docs" / "mass_lifetime_interaction_validation.md"
SEED = 20260603


def load_tolerances() -> dict[str, float | int]:
    return json.loads(CONFIG.read_text(encoding="ascii"))


def safe_log10(value: float, label: str = "value") -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{label} must be positive and finite for log10")
    return math.log10(value)


def finite_unstable_particles() -> pd.DataFrame:
    rows = []
    for particle in load_particles():
        if (
            particle.stability == "stable"
            or particle.lifetime_s is None
            or not math.isfinite(particle.lifetime_s)
            or particle.lifetime_s <= 0.0
        ):
            continue
        mass_e_v = particle.mass_mev * 1.0e6
        gamma_e_v = HBAR / particle.lifetime_s / EV_J
        compton_frequency = particle.mass_mev * MEV_J / H
        n_cycles = compton_frequency * particle.lifetime_s
        rows.append(
            {
                "name": particle.name,
                "family": particle.family,
                "mass_mev": particle.mass_mev,
                "lifetime_s": particle.lifetime_s,
                "stability": particle.stability,
                "dominant_interaction": particle.dominant_interaction,
                "representative_decay_mode_count": particle.representative_decay_mode_count,
                "log_mass": safe_log10(particle.mass_mev, "mass_mev"),
                "log_lifetime": safe_log10(particle.lifetime_s, "lifetime_s"),
                "compton_frequency_hz": compton_frequency,
                "N": n_cycles,
                "Q": n_cycles,
                "logN": safe_log10(n_cycles, "N"),
                "Gamma_eV": gamma_e_v,
                "log_width": safe_log10(gamma_e_v, "Gamma_eV"),
                "relative_width": gamma_e_v / mass_e_v,
            }
        )
    return pd.DataFrame(rows)


def _rank(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    index = 0
    while index < len(values):
        end = index + 1
        while end < len(values) and values[order[end]] == values[order[index]]:
            end += 1
        ranks[order[index:end]] = (index + 1 + end) / 2.0
        index = end
    return ranks


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(_rank(x), _rank(y))


def kendall_tau(x: np.ndarray, y: np.ndarray) -> float:
    concordant = discordant = 0
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            sign = np.sign((x[i] - x[j]) * (y[i] - y[j]))
            if sign > 0:
                concordant += 1
            elif sign < 0:
                discordant += 1
    denominator = concordant + discordant
    return (concordant - discordant) / denominator if denominator else math.nan


def permutation_p_value(
    x: np.ndarray,
    y: np.ndarray,
    metric: Callable[[np.ndarray, np.ndarray], float],
    runs: int,
    rng: np.random.Generator,
) -> float:
    observed = abs(metric(x, y))
    count = 0
    for _ in range(runs):
        if abs(metric(x, rng.permutation(y))) >= observed:
            count += 1
    return (count + 1.0) / (runs + 1.0)


def bootstrap_ci(
    x: np.ndarray,
    y: np.ndarray,
    metric: Callable[[np.ndarray, np.ndarray], float],
    runs: int,
    confidence: float,
    rng: np.random.Generator,
) -> tuple[float, float]:
    values = []
    for _ in range(runs):
        sample = rng.integers(0, len(x), len(x))
        value = metric(x[sample], y[sample])
        if math.isfinite(value):
            values.append(value)
    alpha = (1.0 - confidence) / 2.0
    return (
        float(np.quantile(values, alpha)),
        float(np.quantile(values, 1.0 - alpha)),
    )


def design_matrix(frame: pd.DataFrame, model: str) -> tuple[np.ndarray, list[str]]:
    n = len(frame)
    columns = [np.ones(n)]
    names = ["intercept"]
    interactions = sorted(frame["dominant_interaction"].unique())
    baseline = interactions[0] if interactions else ""
    if model in {"mass", "mass_interaction", "mass_x_interaction"}:
        columns.append(frame["log_mass"].to_numpy())
        names.append("log_mass")
    if model in {"interaction", "mass_interaction", "mass_x_interaction"}:
        for interaction in interactions:
            if interaction == baseline:
                continue
            values = (frame["dominant_interaction"] == interaction).astype(float).to_numpy()
            columns.append(values)
            names.append(f"interaction_{interaction}")
            if model == "mass_x_interaction":
                columns.append(values * frame["log_mass"].to_numpy())
                names.append(f"log_mass_x_{interaction}")
    return np.column_stack(columns), names


def fit_ols(x: np.ndarray, y: np.ndarray) -> dict[str, object]:
    beta = np.linalg.pinv(x) @ y
    fitted = x @ beta
    residuals = y - fitted
    n, p = x.shape
    sse = float(np.sum(residuals**2))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - sse / tss if tss else 0.0
    adjusted = 1.0 - (1.0 - r2) * (n - 1) / (n - p) if n > p else math.nan
    sigma2 = max(sse / n, 1e-300)
    return {
        "beta": beta,
        "fitted": fitted,
        "residuals": residuals,
        "sse": sse,
        "r2": r2,
        "adjusted_r2": adjusted,
        "aic": float(n * math.log(sigma2) + 2 * p),
        "bic": float(n * math.log(sigma2) + p * math.log(n)),
        "n": n,
        "p": p,
    }


def loo_cv_r2(frame: pd.DataFrame, model: str) -> float:
    y = frame["log_lifetime"].to_numpy()
    predictions = np.zeros(len(frame))
    for index in range(len(frame)):
        train = frame.drop(frame.index[index]).reset_index(drop=True)
        test = frame.iloc[[index]].reset_index(drop=True)
        x_train, names = design_matrix(train, model)
        beta = fit_ols(x_train, train["log_lifetime"].to_numpy())["beta"]
        x_test, test_names = design_matrix(test, model)
        aligned = np.zeros((1, len(names)))
        for name_index, name in enumerate(names):
            if name in test_names:
                aligned[0, name_index] = x_test[0, test_names.index(name)]
        predictions[index] = float((aligned @ beta)[0])
    tss = float(np.sum((y - np.mean(y)) ** 2))
    sse = float(np.sum((y - predictions) ** 2))
    return 1.0 - sse / tss if tss else math.nan


def slope_metric(x: np.ndarray, y: np.ndarray) -> float:
    design = np.column_stack([np.ones(len(x)), x])
    return float(fit_ols(design, y)["beta"][1])


def global_correlation(frame: pd.DataFrame, tolerances: dict[str, float | int]) -> pd.DataFrame:
    x = frame["log_mass"].to_numpy()
    y = frame["log_lifetime"].to_numpy()
    runs = int(tolerances["permutation_runs"])
    confidence = float(tolerances["bootstrap_confidence"])
    rng = np.random.default_rng(SEED)
    rows = []
    for name, metric, tolerance in [
        ("pearson", pearson, tolerances["pearson_rho_tolerance"]),
        ("spearman", spearman, tolerances["spearman_rho_tolerance"]),
        ("kendall_tau", kendall_tau, tolerances["spearman_rho_tolerance"]),
    ]:
        low, high = bootstrap_ci(x, y, metric, runs, confidence, rng)
        rows.append(
            {
                "metric": name,
                "coefficient": metric(x, y),
                "p_value": permutation_p_value(x, y, metric, runs, rng),
                "n": len(frame),
                "ci_low": low,
                "ci_high": high,
                "bootstrap_confidence": confidence,
                "coefficient_tolerance": tolerance,
                "p_value_threshold": tolerances["p_value_threshold"],
            }
        )
    return pd.DataFrame(rows)


def interaction_regression(frame: pd.DataFrame, tolerances: dict[str, float | int]) -> pd.DataFrame:
    rows = []
    rng = np.random.default_rng(SEED + 1)
    min_n = int(tolerances["min_group_n_for_regression"])
    for interaction, group in frame.groupby("dominant_interaction"):
        x = group["log_mass"].to_numpy()
        y = group["log_lifetime"].to_numpy()
        if len(group) < min_n:
            rows.append(
                {
                    "interaction": interaction,
                    "status": "insufficient_data",
                    "n": len(group),
                    "alpha": math.nan,
                    "intercept": math.nan,
                    "r2": math.nan,
                    "adjusted_r2": math.nan,
                    "pearson_rho": math.nan,
                    "spearman_rho": math.nan,
                    "p_value": math.nan,
                    "alpha_ci_low": math.nan,
                    "alpha_ci_high": math.nan,
                    "bootstrap_confidence": tolerances["bootstrap_confidence"],
                    "min_group_n_for_regression": min_n,
                }
            )
            continue
        fit = fit_ols(np.column_stack([np.ones(len(x)), x]), y)
        low, high = bootstrap_ci(
            x,
            y,
            slope_metric,
            int(tolerances["permutation_runs"]),
            float(tolerances["bootstrap_confidence"]),
            rng,
        )
        rows.append(
            {
                "interaction": interaction,
                "status": "fit",
                "n": len(group),
                "alpha": float(fit["beta"][1]),
                "intercept": float(fit["beta"][0]),
                "r2": fit["r2"],
                "adjusted_r2": fit["adjusted_r2"],
                "pearson_rho": pearson(x, y),
                "spearman_rho": spearman(x, y),
                "p_value": permutation_p_value(x, y, pearson, int(tolerances["permutation_runs"]), rng),
                "alpha_ci_low": low,
                "alpha_ci_high": high,
                "bootstrap_confidence": tolerances["bootstrap_confidence"],
                "min_group_n_for_regression": min_n,
            }
        )
    return pd.DataFrame(rows).sort_values("interaction")


def expectation_comparison(regression: pd.DataFrame, tolerances: dict[str, float | int]) -> pd.DataFrame:
    rows = []
    for _, row in regression.iterrows():
        expected = float(tolerances["weak_exponent_expected"]) if row["interaction"] == "weak" else math.nan
        tolerance = float(tolerances["weak_exponent_tolerance"]) if row["interaction"] == "weak" else math.nan
        error = abs(row["alpha"] - expected) if math.isfinite(expected) and math.isfinite(row["alpha"]) else math.nan
        rows.append(
            {
                "interaction": row["interaction"],
                "fitted_alpha": row["alpha"],
                "expected_alpha": expected,
                "tolerance": tolerance,
                "absolute_error": error,
                "relative_error": error / abs(expected) if math.isfinite(error) and expected else math.nan,
                "compatible_with_expectation": bool(error <= tolerance) if math.isfinite(error) else False,
                "caution_note": (
                    "Weak -5 scaling is approximate and muon-like; it is not universal."
                    if row["interaction"] == "weak"
                    else "No fixed exponent was imposed for this interaction."
                ),
            }
        )
    return pd.DataFrame(rows)


def model_rows(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, model in [
        ("Model 1: log_lifetime ~ log_mass", "mass"),
        ("Model 2: log_lifetime ~ interaction", "interaction"),
        ("Model 3: log_lifetime ~ log_mass + interaction", "mass_interaction"),
        ("Model 4: log_lifetime ~ log_mass * interaction", "mass_x_interaction"),
    ]:
        x, _ = design_matrix(frame, model)
        fit = fit_ols(x, frame["log_lifetime"].to_numpy())
        rows.append(
            {
                "test": "model_comparison",
                "model": label,
                "aic": fit["aic"],
                "bic": fit["bic"],
                "r2": fit["r2"],
                "adjusted_r2": fit["adjusted_r2"],
                "cross_validated_r2": loo_cv_r2(frame, model),
                "n": len(frame),
            }
        )
    return pd.DataFrame(rows)


def tautology_decomposition(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    x_mass = frame["log_mass"].to_numpy()
    y_tau = frame["log_lifetime"].to_numpy()
    y_n = frame["logN"].to_numpy()
    mass_tau = fit_ols(np.column_stack([np.ones(len(frame)), x_mass]), y_tau)
    mass_n = fit_ols(np.column_stack([np.ones(len(frame)), x_mass]), y_n)
    rows.extend(
        [
            {
                "test": "mass_vs_lifetime",
                "model": "log_lifetime ~ log_mass",
                "coefficient": pearson(x_mass, y_tau),
                "r2": mass_tau["r2"],
                "adjusted_r2": mass_tau["adjusted_r2"],
                "aic": mass_tau["aic"],
                "bic": mass_tau["bic"],
                "cross_validated_r2": loo_cv_r2(frame, "mass"),
                "answer": "Measured lifetime is not defined from mass, but this model tests only mass association.",
                "n": len(frame),
            },
            {
                "test": "mass_vs_N",
                "model": "logN ~ log_mass",
                "coefficient": pearson(x_mass, y_n),
                "r2": mass_n["r2"],
                "adjusted_r2": mass_n["adjusted_r2"],
                "aic": mass_n["aic"],
                "bic": mass_n["bic"],
                "cross_validated_r2": math.nan,
                "answer": "Partly tautological because logN includes log_mass plus log_lifetime and a constant.",
                "n": len(frame),
            },
            {
                "test": "mass_vs_Q",
                "model": "logQ ~ log_mass",
                "coefficient": pearson(x_mass, y_n),
                "r2": mass_n["r2"],
                "adjusted_r2": mass_n["adjusted_r2"],
                "aic": mass_n["aic"],
                "bic": mass_n["bic"],
                "cross_validated_r2": math.nan,
                "answer": "Q is treated here as N by task definition, so it has the same tautological mass term.",
                "n": len(frame),
            },
        ]
    )
    rows.extend(model_rows(frame).to_dict("records"))
    residual_frame = frame.copy()
    residual_frame["mass_residual"] = mass_tau["residuals"]
    interaction_resid = fit_ols(design_matrix(frame, "interaction")[0], mass_tau["residuals"])
    rows.append(
        {
            "test": "interaction_labels_explain_mass_residuals",
            "model": "residual(log_lifetime ~ log_mass) ~ interaction",
            "coefficient": math.nan,
            "r2": interaction_resid["r2"],
            "adjusted_r2": interaction_resid["adjusted_r2"],
            "aic": interaction_resid["aic"],
            "bic": interaction_resid["bic"],
            "cross_validated_r2": math.nan,
            "answer": "Interaction-label residual structure is non-tautological if above random-label controls.",
            "n": len(frame),
        }
    )
    return pd.DataFrame(rows)


def interaction_significance(frame: pd.DataFrame, tolerances: dict[str, float | int]) -> pd.DataFrame:
    y = frame["log_lifetime"].to_numpy()
    fit1 = fit_ols(design_matrix(frame, "mass")[0], y)
    fit3 = fit_ols(design_matrix(frame, "mass_interaction")[0], y)
    observed = float(fit3["r2"] - fit1["r2"])
    rng = np.random.default_rng(SEED + 2)
    runs = int(tolerances["permutation_runs"])
    count = 0
    for _ in range(runs):
        shuffled = frame.copy()
        shuffled["dominant_interaction"] = rng.permutation(shuffled["dominant_interaction"].to_numpy())
        perm_fit = fit_ols(design_matrix(shuffled, "mass_interaction")[0], y)
        if float(perm_fit["r2"] - fit1["r2"]) >= observed:
            count += 1
    p_value = (count + 1.0) / (runs + 1.0)
    return pd.DataFrame(
        [
            {
                "model_comparison": "Model 3 vs Model 1",
                "observed_delta_R2": observed,
                "observed_delta_AIC": float(fit3["aic"] - fit1["aic"]),
                "observed_delta_BIC": float(fit3["bic"] - fit1["bic"]),
                "permutation_p_value": p_value,
                "n_permutations": runs,
                "passes_p05": bool(p_value < float(tolerances["p_value_threshold"])),
                "interpretation": (
                    "Interaction labels explain lifetime beyond mass at p<0.05."
                    if p_value < float(tolerances["p_value_threshold"])
                    else "Interaction-label improvement is not significant at p<0.05."
                ),
                "p_value_threshold": tolerances["p_value_threshold"],
                "n": len(frame),
            }
        ]
    )


def control_rows(frame: pd.DataFrame, tolerances: dict[str, float | int]) -> pd.DataFrame:
    rng = np.random.default_rng(SEED + 3)
    observed_rho = pearson(frame["log_mass"].to_numpy(), frame["log_lifetime"].to_numpy())
    observed_fit1 = fit_ols(design_matrix(frame, "mass")[0], frame["log_lifetime"].to_numpy())
    observed_fit3 = fit_ols(design_matrix(frame, "mass_interaction")[0], frame["log_lifetime"].to_numpy())
    observed_improvement = float(observed_fit3["r2"] - observed_fit1["r2"])
    rows = []
    runs = int(tolerances["permutation_runs"])
    for control in [
        "shuffle_lifetimes",
        "shuffle_masses",
        "shuffle_interaction_labels",
        "shuffle_family_labels",
        "random_log_lifetime_same_range",
        "random_log_mass_same_range",
    ]:
        rhos = []
        slopes = []
        improvements = []
        for _ in range(runs):
            controlled = frame.copy()
            if control == "shuffle_lifetimes":
                controlled["log_lifetime"] = rng.permutation(controlled["log_lifetime"].to_numpy())
            elif control == "shuffle_masses":
                controlled["log_mass"] = rng.permutation(controlled["log_mass"].to_numpy())
            elif control == "shuffle_interaction_labels":
                controlled["dominant_interaction"] = rng.permutation(controlled["dominant_interaction"].to_numpy())
            elif control == "shuffle_family_labels":
                controlled["family"] = rng.permutation(controlled["family"].to_numpy())
            elif control == "random_log_lifetime_same_range":
                low, high = controlled["log_lifetime"].min(), controlled["log_lifetime"].max()
                controlled["log_lifetime"] = rng.uniform(low, high, len(controlled))
            elif control == "random_log_mass_same_range":
                low, high = controlled["log_mass"].min(), controlled["log_mass"].max()
                controlled["log_mass"] = rng.uniform(low, high, len(controlled))
            x = controlled["log_mass"].to_numpy()
            y = controlled["log_lifetime"].to_numpy()
            fit1 = fit_ols(design_matrix(controlled, "mass")[0], y)
            fit3 = fit_ols(design_matrix(controlled, "mass_interaction")[0], y)
            rhos.append(pearson(x, y))
            slopes.append(slope_metric(x, y))
            improvements.append(float(fit3["r2"] - fit1["r2"]))
        rows.append(
            {
                "control": control,
                "rho": float(np.mean(rhos)),
                "rho_ci_low": float(np.quantile(rhos, 0.025)),
                "rho_ci_high": float(np.quantile(rhos, 0.975)),
                "regression_alpha": float(np.mean(slopes)),
                "alpha_ci_low": float(np.quantile(slopes, 0.025)),
                "alpha_ci_high": float(np.quantile(slopes, 0.975)),
                "model3_delta_r2": float(np.mean(improvements)),
                "model3_delta_r2_ci_low": float(np.quantile(improvements, 0.025)),
                "model3_delta_r2_ci_high": float(np.quantile(improvements, 0.975)),
                "empirical_p_value": (sum(abs(value) >= abs(observed_rho) for value in rhos) + 1.0) / (runs + 1.0),
                "model3_improvement_empirical_p_value": (sum(value >= observed_improvement for value in improvements) + 1.0) / (runs + 1.0),
                "observed_rho_reference": observed_rho,
                "observed_model3_delta_r2_reference": observed_improvement,
                "n_resamples": runs,
                "n": len(frame),
                "p_value_threshold": tolerances["p_value_threshold"],
            }
        )
    return pd.DataFrame(rows)


def add_tolerance_columns(frame: pd.DataFrame, tolerances: dict[str, float | int]) -> pd.DataFrame:
    frame = frame.copy()
    frame["relative_numeric_tolerance"] = tolerances["relative_numeric_tolerance"]
    frame["absolute_numeric_tolerance"] = tolerances["absolute_numeric_tolerance"]
    return frame


def generate_figures(
    frame: pd.DataFrame,
    global_corr: pd.DataFrame,
    by_interaction: pd.DataFrame,
    decomposition: pd.DataFrame,
    significance: pd.DataFrame,
    controls: pd.DataFrame,
) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    colors = {name: color for name, color in zip(sorted(frame["dominant_interaction"].unique()), plt.cm.tab10.colors)}

    def scatter_by_interaction(y_col: str, ylabel: str, path: str) -> None:
        fig, ax = plt.subplots(figsize=(8, 5))
        for interaction, group in frame.groupby("dominant_interaction"):
            ax.scatter(group["log_mass"], group[y_col], label=interaction, color=colors[interaction])
            if len(group) >= 3 and y_col == "log_lifetime":
                alpha = by_interaction.loc[by_interaction["interaction"] == interaction, "alpha"].iloc[0]
                intercept = by_interaction.loc[by_interaction["interaction"] == interaction, "intercept"].iloc[0]
                xs = np.linspace(group["log_mass"].min(), group["log_mass"].max(), 40)
                ax.plot(xs, intercept + alpha * xs, color=colors[interaction], linestyle="--")
        if y_col == "log_lifetime":
            fit = fit_ols(np.column_stack([np.ones(len(frame)), frame["log_mass"]]), frame["log_lifetime"].to_numpy())
            xs = np.linspace(frame["log_mass"].min(), frame["log_mass"].max(), 80)
            ax.plot(xs, float(fit["beta"][0]) + float(fit["beta"][1]) * xs, color="black", linewidth=2, label="global")
        ax.set_xlabel("log10(mass / MeV)")
        ax.set_ylabel(ylabel)
        ax.legend()
        fig.tight_layout()
        fig.savefig(OUTPUT / path, dpi=160)
        plt.close(fig)

    scatter_by_interaction("log_lifetime", "log10(lifetime / s)", "01_log_mass_vs_log_lifetime.png")
    scatter_by_interaction("log_width", "log10(width / eV)", "02_log_mass_vs_log_width.png")

    fit_mass = fit_ols(np.column_stack([np.ones(len(frame)), frame["log_mass"]]), frame["log_lifetime"].to_numpy())
    residuals = fit_mass["residuals"]
    fig, ax = plt.subplots(figsize=(8, 5))
    for interaction, group in frame.assign(residual=residuals).groupby("dominant_interaction"):
        ax.scatter(group["log_mass"], group["residual"], label=interaction, color=colors[interaction])
    ax.axhline(0.0, color="black", linewidth=1)
    ax.set_xlabel("log10(mass / MeV)")
    ax.set_ylabel("lifetime residual after global mass fit")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT / "03_lifetime_residuals_after_mass_fit.png", dpi=160)
    plt.close(fig)

    fitted = by_interaction[by_interaction["status"] == "fit"]
    fig, ax = plt.subplots(figsize=(8, 5))
    yerr = np.vstack([fitted["alpha"] - fitted["alpha_ci_low"], fitted["alpha_ci_high"] - fitted["alpha"]])
    ax.errorbar(fitted["interaction"], fitted["alpha"], yerr=yerr, fmt="o", capsize=4)
    ax.axhline(-5.0, color="gray", linestyle=":", label="weak reference -5")
    ax.set_ylabel("fitted alpha in log10(tau) = a + alpha log10(mass)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT / "04_alpha_by_interaction.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(["observed"] + controls["control"].tolist(), [global_corr.loc[0, "coefficient"]] + controls["rho"].tolist())
    ax.tick_params(axis="x", rotation=45)
    ax.set_ylabel("Pearson rho")
    fig.tight_layout()
    fig.savefig(OUTPUT / "05_observed_rho_vs_controls.png", dpi=160)
    plt.close(fig)

    model_rows_only = decomposition[decomposition["test"] == "model_comparison"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, col in zip(axes, ["adjusted_r2", "aic", "bic"]):
        ax.bar(model_rows_only["model"].str.extract(r"(Model \\d)")[0], model_rows_only[col])
        ax.set_title(col)
    fig.tight_layout()
    fig.savefig(OUTPUT / "06_model_comparison.png", dpi=160)
    plt.close(fig)

    rng = np.random.default_rng(SEED + 4)
    fit1 = fit_ols(design_matrix(frame, "mass")[0], frame["log_lifetime"].to_numpy())
    deltas = []
    for _ in range(1000):
        shuffled = frame.copy()
        shuffled["dominant_interaction"] = rng.permutation(shuffled["dominant_interaction"].to_numpy())
        fit3 = fit_ols(design_matrix(shuffled, "mass_interaction")[0], shuffled["log_lifetime"].to_numpy())
        deltas.append(float(fit3["r2"] - fit1["r2"]))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(deltas, bins=30, color="#4c78a8")
    ax.axvline(significance.loc[0, "observed_delta_R2"], color="black", label="observed")
    ax.set_xlabel("delta R2 from shuffled interaction labels")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT / "07_interaction_permutation_null.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    constant = math.log10(MEV_J / H)
    ax.scatter(frame["log_mass"] + frame["log_lifetime"] + constant, frame["logN"])
    low, high = frame["logN"].min(), frame["logN"].max()
    ax.plot([low, high], [low, high], color="black")
    ax.set_xlabel("log_mass + log_lifetime + log10(MeV_J / h)")
    ax.set_ylabel("logN")
    fig.tight_layout()
    fig.savefig(OUTPUT / "08_logN_decomposition.png", dpi=160)
    plt.close(fig)

    variables = ["log_mass", "log_lifetime", "logN", "log_width", "relative_width"]
    matrix = frame.set_index("name")[variables]
    normalized = (matrix - matrix.mean()) / matrix.std(ddof=0)
    fig, ax = plt.subplots(figsize=(8, 8))
    image = ax.imshow(normalized, aspect="auto", cmap="coolwarm")
    ax.set_xticks(range(len(variables)), variables, rotation=45, ha="right")
    ax.set_yticks(range(len(matrix.index)), matrix.index)
    fig.colorbar(image, ax=ax, label="z score")
    fig.tight_layout()
    fig.savefig(OUTPUT / "09_particle_variable_heatmap.png", dpi=160)
    plt.close(fig)

    residual_frame = frame.assign(residual=residuals, abs_residual=np.abs(residuals)).sort_values("abs_residual", ascending=False)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(residual_frame["log_mass"], residual_frame["residual"], color="#f58518")
    for _, row in residual_frame.head(8).iterrows():
        ax.annotate(row["name"], (row["log_mass"], row["residual"]), fontsize=8)
    ax.axhline(0.0, color="black", linewidth=1)
    ax.set_xlabel("log10(mass / MeV)")
    ax.set_ylabel("global mass-fit residual")
    fig.tight_layout()
    fig.savefig(OUTPUT / "10_largest_residual_outliers.png", dpi=160)
    plt.close(fig)


def report(
    frame: pd.DataFrame,
    global_corr: pd.DataFrame,
    by_interaction: pd.DataFrame,
    decomposition: pd.DataFrame,
    significance: pd.DataFrame,
    tolerances: dict[str, float | int],
) -> str:
    pearson_row = global_corr[global_corr["metric"] == "pearson"].iloc[0]
    spearman_row = global_corr[global_corr["metric"] == "spearman"].iloc[0]
    mass_model = decomposition[decomposition["model"] == "log_lifetime ~ log_mass"].iloc[0]
    model3 = decomposition[decomposition["model"] == "Model 3: log_lifetime ~ log_mass + interaction"].iloc[0]
    conclusion_supported = (
        pearson_row["p_value"] < float(tolerances["p_value_threshold"])
        and significance.loc[0, "passes_p05"]
    )
    slope = slope_metric(frame["log_mass"].to_numpy(), frame["log_lifetime"].to_numpy())
    interaction_lines = "\n".join(
        f"- {row['interaction']}: status={row['status']}, n={int(row['n'])}, alpha={row['alpha']:.6g}, "
        f"95% CI=[{row['alpha_ci_low']:.6g}, {row['alpha_ci_high']:.6g}]"
        if row["status"] == "fit"
        else f"- {row['interaction']}: status=insufficient_data, n={int(row['n'])}"
        for _, row in by_interaction.iterrows()
    )
    valid_conclusion = (
        "After removing the definitional contribution of Compton frequency to N, the measured lifetime still "
        "shows a statistically meaningful relationship with mass and interaction class. This supports the "
        "interpretation that lifetime is governed by energy scale, phase space, coupling strength, and decay "
        "topology, not by mass alone."
        if conclusion_supported
        else "The current dataset is too small or too biased to support the stronger conclusion that interaction class adds statistically robust explanatory power beyond mass."
    )
    return f"""# Mass-Lifetime Interaction Validation

Tolerances used: coefficient tolerance +/-{tolerances['pearson_rho_tolerance']}, p-value threshold {tolerances['p_value_threshold']}, bootstrap confidence {tolerances['bootstrap_confidence']}, weak exponent reference {tolerances['weak_exponent_expected']} +/-{tolerances['weak_exponent_tolerance']}, minimum regression group n={tolerances['min_group_n_for_regression']}, permutation runs={tolerances['permutation_runs']}.

## Section 1 - What was tested

This validation tests the mass-lifetime correlation for finite-lifetime unstable particles only, with n={len(frame)}. Stable particles, missing lifetimes, infinite lifetimes, and manually truncated stable values are excluded by default.

The Compton frequency is `f_C = m c^2 / h`, and persistence cycles are `N = f_C tau`. Because `f_C` is proportional to mass, correlations involving `N` and mass are partly tautological. Lifetime itself is not tautological because `tau` is measured independently of the definition of `f_C`.

Interaction type matters because decay width follows `Gamma = hbar / tau` and, schematically, `Gamma` depends on matrix elements and phase space. Weak three-body beta-like decays can show an approximate `Q^5` or muon-like `m^-5` lifetime scaling, but that exponent is not universal.

## Section 2 - Results

- Global Pearson rho: {pearson_row['coefficient']:.6g}, p={pearson_row['p_value']:.6g}, n={int(pearson_row['n'])}, 95% bootstrap CI=[{pearson_row['ci_low']:.6g}, {pearson_row['ci_high']:.6g}], tolerance +/-{pearson_row['coefficient_tolerance']}.
- Global Spearman rho: {spearman_row['coefficient']:.6g}, p={spearman_row['p_value']:.6g}, n={int(spearman_row['n'])}, 95% bootstrap CI=[{spearman_row['ci_low']:.6g}, {spearman_row['ci_high']:.6g}], tolerance +/-{spearman_row['coefficient_tolerance']}.
- Global regression slope alpha: {slope:.6g}, n={len(frame)}.
- Mass-only adjusted R2: {mass_model['adjusted_r2']:.6g}, AIC={mass_model['aic']:.6g}, BIC={mass_model['bic']:.6g}.
- Mass plus interaction adjusted R2: {model3['adjusted_r2']:.6g}, AIC={model3['aic']:.6g}, BIC={model3['bic']:.6g}.
- Interaction permutation p-value for Model 3 vs Model 1: {significance.loc[0, 'permutation_p_value']:.6g}, n={len(frame)}, permutation runs={int(significance.loc[0, 'n_permutations'])}.

Per-interaction slopes:

{interaction_lines}

## Section 3 - Physical interpretation

Mass tends to open phase space and often increases decay width, reducing lifetime. But lifetime is not determined by mass alone. It also depends on coupling constants, available decay channels, phase space, conservation laws, selection rules, and interaction type.

The weak exponent comparison uses alpha approximately -5 only as a loose reference motivated by muon-like weak decay scaling. It should not be applied as a universal weak-decay law.

## Section 4 - Tautology audit

Tautological or partly tautological:

- `N` includes mass through `f_C`.
- `logN = log_tau + log_mass + log10(MeV_J / h)`.
- Correlations between mass and `N` reuse mass by construction.

Not tautological:

- Correlation between mass and measured lifetime.
- Interaction-dependent residual structure after fitting lifetime against mass.
- Improvement of models when interaction labels are included, subject to permutation controls.

In this run, mass explains R2={mass_model['r2']:.6g} of `log_lifetime`, while mass explains the reported `logN` variance partly through the definitional mass term. Interaction labels explain additional variance only if the permutation p-value is below {tolerances['p_value_threshold']}.

## Section 5 - Valid conclusion

{valid_conclusion}

This conclusion is evaluated with n={len(frame)}, coefficient tolerance +/-{tolerances['pearson_rho_tolerance']}, p-value threshold {tolerances['p_value_threshold']}, and bootstrap confidence {tolerances['bootstrap_confidence']}.

## Section 6 - Warnings

- Small dataset warning: per-interaction groups are small, and groups with n < {tolerances['min_group_n_for_regression']} are marked insufficient.
- PDG/catalog selection bias warning: the sample is curated and not a complete particle catalogue.
- Stable particle exclusion warning: stable particles are excluded from default regression to avoid artificial truncation.
- Mixed interaction ambiguity warning: `mixed` labels combine multiple mechanisms and should not be overinterpreted.
- No new physics claim: these are descriptive validation checks on known measured quantities.

## Reproduction

```powershell
python -m src.mass_lifetime_interactions
```
"""


def generate() -> dict[str, pd.DataFrame]:
    tolerances = load_tolerances()
    DATA.mkdir(parents=True, exist_ok=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    frame = finite_unstable_particles()
    global_corr = add_tolerance_columns(global_correlation(frame, tolerances), tolerances)
    by_interaction = add_tolerance_columns(interaction_regression(frame, tolerances), tolerances)
    expectation = add_tolerance_columns(expectation_comparison(by_interaction, tolerances), tolerances)
    decomposition = add_tolerance_columns(tautology_decomposition(frame), tolerances)
    significance = add_tolerance_columns(interaction_significance(frame, tolerances), tolerances)
    controls = add_tolerance_columns(control_rows(frame, tolerances), tolerances)
    global_corr.to_csv(DATA / "mass_lifetime_global_correlation.csv", index=False)
    by_interaction.to_csv(DATA / "mass_lifetime_by_interaction.csv", index=False)
    expectation.to_csv(DATA / "scaling_expectation_comparison.csv", index=False)
    decomposition.to_csv(DATA / "tautology_decomposition.csv", index=False)
    significance.to_csv(DATA / "interaction_significance.csv", index=False)
    controls.to_csv(DATA / "mass_lifetime_controls.csv", index=False)
    generate_figures(frame, global_corr, by_interaction, decomposition, significance, controls)
    DOC.write_text(report(frame, global_corr, by_interaction, decomposition, significance, tolerances), encoding="ascii")
    return {
        "frame": frame,
        "global": global_corr,
        "by_interaction": by_interaction,
        "expectation": expectation,
        "decomposition": decomposition,
        "significance": significance,
        "controls": controls,
    }


if __name__ == "__main__":
    result = generate()
    print(f"Generated mass-lifetime validation for {len(result['frame'])} finite-lifetime particles in {OUTPUT}")
