"""Main script to run all experiments for Chapter 5.

This script runs MMCP and NW models for all 5 time periods,
computes LOOCV metrics, and generates all figures.

Usage:
    python -m bayesian_statistics.experiments.run_all --help
    python -m bayesian_statistics.experiments.run_all --run-models
    python -m bayesian_statistics.experiments.run_all --run-loocv
    python -m bayesian_statistics.experiments.run_all --generate-figures
    python -m bayesian_statistics.experiments.run_all --all
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any, Dict

import numpy as np

from .config import ExperimentConfig
from .evaluation import MMCPLOOCVEvaluator, NWLOOCVEvaluator
from .models.mmcp_runner import MMCPRunner
from .models.nw_runner import NWRunner
from .output import ProgressManager, set_progress_manager


def _compute_grid_probs_std(mcmc_results, n_samples: int = 50) -> np.ndarray:
    """Compute posterior standard deviation of grid probabilities.

    Uses conditional sampling to get multiple probability predictions,
    then computes the standard deviation across samples.

    Parameters
    ----------
    mcmc_results : MarkedPointProcessResults
        MCMC results containing beta samples.
    n_samples : int
        Number of conditional samples to draw.

    Returns
    -------
    np.ndarray
        Posterior standard deviation (K, n_grid).
    """
    # Get multiple samples using conditional sampling
    probs_samples = []
    for i in range(n_samples):
        probs = mcmc_results.predict_probabilities(
            location="grid", sample_conditional=True, seed=42 + i
        )
        probs_samples.append(probs)

    probs_array = np.stack(probs_samples, axis=0)  # (n_samples, K, n_grid)
    return probs_array.std(axis=0)  # (K, n_grid)


def _compute_p0_from_dataset(dataset) -> np.ndarray | None:
    """Compute distance prior p0 from dataset's distance features.

    The distance_features_grid contains log-ratios: log(p_k) - log(p_K)
    where p_K is the baseline category. We convert these back to probabilities
    using softmax transformation.

    Parameters
    ----------
    dataset : MarkedPointProcessDataset
        Dataset containing distance_features_grid.

    Returns
    -------
    np.ndarray or None
        Prior probabilities (n_grid, K), or None if no distance features.
        Format matches plot_distance_prior expectations.
    """
    if dataset.distance_features_grid is None:
        return None

    # distance_features_grid has shape (n_grid, K-1)
    # Contains g_k = log(p_k / p_K) for k = 1, ..., K-1
    g = dataset.distance_features_grid  # (n_grid, K-1)

    # Softmax inversion:
    # p_k = exp(g_k) / (1 + sum(exp(g_k))) for k = 1, ..., K-1
    # p_K = 1 / (1 + sum(exp(g_k)))
    exp_g = np.exp(g)  # (n_grid, K-1)
    sum_exp_g = exp_g.sum(axis=1)  # (n_grid,)
    denom = 1.0 + sum_exp_g  # (n_grid,)

    p_k = exp_g / denom[:, np.newaxis]  # (n_grid, K-1)
    p_K = 1.0 / denom  # (n_grid,)

    # Combine: p0 = [p_1, p_2, ..., p_{K-1}, p_K] with shape (n_grid, K)
    p0 = np.column_stack([p_k, p_K])  # (n_grid, K)

    return p0


def _check_grid_sizes_match(mmcp_probs: np.ndarray, nw_probs: np.ndarray) -> bool:
    """Check if MMCP and NW grid sizes match.

    Parameters
    ----------
    mmcp_probs : np.ndarray
        MMCP grid probabilities (K, n_grid).
    nw_probs : np.ndarray
        NW grid probabilities (K, n_grid).

    Returns
    -------
    bool
        True if grid sizes match, False otherwise.
    """
    return mmcp_probs.shape[1] == nw_probs.shape[1]


def load_preprocessor(config: ExperimentConfig):
    """Load and prepare the data preprocessor.

    Parameters
    ----------
    config : ExperimentConfig
        Experiment configuration.

    Returns
    -------
    ObsidianDataPreprocessor
        Loaded preprocessor.
    """
    from bayesian_statistics.models.preprocessing.data_preprocessor import (
        ObsidianDataPreprocessor,
    )

    preprocessor = ObsidianDataPreprocessor(str(config.data_dir), scale_variables=True)
    preprocessor.load_data()
    return preprocessor


def run_mmcp_all_periods(
    config: ExperimentConfig,
    preprocessor,
    output_dir: Path,
    progress: ProgressManager,
) -> Dict[int, Dict[str, Any]]:
    """Run MMCP model for all periods and save results.

    Parameters
    ----------
    config : ExperimentConfig
        Experiment configuration.
    preprocessor : ObsidianDataPreprocessor
        Data preprocessor.
    output_dir : Path
        Directory to save results.
    progress : ProgressManager
        Progress manager for output control.

    Returns
    -------
    dict
        Results for all periods.
    """
    progress.section("Running MMCP model for all periods")

    runner = MMCPRunner(config=config, progress=progress)
    results = runner.run_all_periods(preprocessor)

    # Save results
    results_dir = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    for period, result in results.items():
        period_dir = results_dir / f"mmcp_period_{period}"
        period_dir.mkdir(exist_ok=True)

        np.save(period_dir / "site_probs.npy", result["site_probs"])
        np.save(period_dir / "grid_probs.npy", result["grid_probs"])

        # Save effects
        for effect_name, effect_data in result["effects"].items():
            np.save(period_dir / f"effect_{effect_name}.npy", effect_data)

        # Save lambda_star samples
        np.save(
            period_dir / "lambda_star.npy",
            result["results"].lambda_star_samples,
        )

        # Compute and save posterior standard deviation at grid
        # Using sample_conditional=True for proper uncertainty quantification
        mcmc_results = result["results"]
        grid_probs_std = _compute_grid_probs_std(mcmc_results, n_samples=50)
        np.save(period_dir / "grid_probs_std.npy", grid_probs_std)

        # Compute and save site existence probability q(s) = sigmoid(η_int)
        # intensity = λ* × q(s), so q(s) = intensity / λ*
        intensity = mcmc_results.predict_intensity(location="grid")
        lambda_star_mean = mcmc_results.lambda_star_samples.mean()
        site_probability = intensity / lambda_star_mean
        np.save(period_dir / "site_probability.npy", site_probability)

    progress.info(f"\nMMCP results saved to {results_dir}")
    return results


def run_nw_all_periods(
    config: ExperimentConfig,
    preprocessor,
    output_dir: Path,
    progress: ProgressManager,
) -> Dict[int, Dict[str, Any]]:
    """Run NW model for all periods and save results.

    Parameters
    ----------
    config : ExperimentConfig
        Experiment configuration.
    preprocessor : ObsidianDataPreprocessor
        Data preprocessor.
    output_dir : Path
        Directory to save results.
    progress : ProgressManager
        Progress manager for output control.

    Returns
    -------
    dict
        Results for all periods.
    """
    progress.section("Running NW model for all periods")

    runner = NWRunner(config=config, progress=progress)
    results = runner.run_all_periods(preprocessor)

    # Save results
    results_dir = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    for period, result in results.items():
        period_dir = results_dir / f"nw_period_{period}"
        period_dir.mkdir(exist_ok=True)

        np.save(period_dir / "site_probs.npy", result["site_probs"])
        np.save(period_dir / "grid_probs.npy", result["grid_probs"])

    progress.info(f"\nNW results saved to {results_dir}")
    return results


def run_loocv_evaluation(
    config: ExperimentConfig,
    preprocessor,
    output_dir: Path,
    progress: ProgressManager,
) -> Dict[str, Dict[int, Dict[str, Any]]]:
    """Run LOOCV evaluation for MMCP and NW models.

    Parameters
    ----------
    config : ExperimentConfig
        Experiment configuration.
    preprocessor : ObsidianDataPreprocessor
        Data preprocessor.
    output_dir : Path
        Directory to save results.
    progress : ProgressManager
        Progress manager for output control.

    Returns
    -------
    dict
        LOOCV results for both models.
    """
    progress.section("Running LOOCV evaluation")

    results = {}

    # NW LOOCV (faster, do first)
    progress.subsection("NW LOOCV")
    nw_evaluator = NWLOOCVEvaluator(config, progress=progress)
    results["nw"] = nw_evaluator.evaluate_all_periods(preprocessor)

    # MMCP LOOCV (slower due to MCMC)
    progress.subsection("MMCP LOOCV")
    mmcp_runner = MMCPRunner(config, progress=progress)
    mmcp_evaluator = MMCPLOOCVEvaluator(config, progress=progress)
    results["mmcp"] = mmcp_evaluator.evaluate_all_periods(preprocessor, mmcp_runner)

    # Save summary table (Table 5.2)
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    save_loocv_table(results, config, tables_dir / "table_5_2_loocv.csv", progress)

    return results


def save_loocv_table(
    results: Dict[str, Dict[int, Dict[str, Any]]],
    config: ExperimentConfig,
    output_path: Path,
    progress: ProgressManager | None = None,
) -> None:
    """Save LOOCV results as CSV table.

    Parameters
    ----------
    results : dict
        LOOCV results for both models.
    config : ExperimentConfig
        Experiment configuration.
    output_path : Path
        Path to save CSV file.
    progress : ProgressManager, optional
        Progress manager for output control.
    """
    rows = []
    header = ["時期", "MMCP平均", "MMCP標準偏差", "NW平均", "NW標準偏差"]

    for period in config.time_periods.keys():
        period_name = config.time_periods[period]

        mmcp_result = results.get("mmcp", {}).get(period, {})
        nw_result = results.get("nw", {}).get(period, {})

        row = [
            period_name,
            f"{mmcp_result.get('mean_aitchison_distance', float('nan')):.4f}",
            f"{mmcp_result.get('std_aitchison_distance', float('nan')):.4f}",
            f"{nw_result.get('mean_aitchison_distance', float('nan')):.4f}",
            f"{nw_result.get('std_aitchison_distance', float('nan')):.4f}",
        ]
        rows.append(row)

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

    if progress:
        progress.info(f"\nLOOCV table saved to {output_path}")


def save_data_summary_table(
    preprocessor,
    config: ExperimentConfig,
    output_dir: Path,
    progress: ProgressManager | None = None,
) -> None:
    """Save data summary table (Table 5.1).

    Parameters
    ----------
    preprocessor : ObsidianDataPreprocessor
        Data preprocessor.
    config : ExperimentConfig
        Experiment configuration.
    output_dir : Path
        Directory to save tables.
    progress : ProgressManager, optional
        Progress manager for output control.
    """
    import polars as pl

    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for period, period_name in config.time_periods.items():
        period_data = preprocessor.df_obsidian.filter(pl.col("時期") == period)
        site_count = period_data["遺跡ID"].n_unique()
        total_count = period_data.height
        rows.append({
            "時期": period_name,
            "遺跡数": site_count,
            "総出土数": total_count,
        })

    output_path = tables_dir / "table_5_1_data_summary.csv"
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["時期", "遺跡数", "総出土数"])
        writer.writeheader()
        writer.writerows(rows)

    if progress:
        progress.info(f"Data summary table saved to {output_path}")


def save_lambda_posterior_table(
    config: ExperimentConfig,
    output_dir: Path,
    progress: ProgressManager | None = None,
) -> None:
    """Save lambda* posterior summary table (Table 5.3).

    Parameters
    ----------
    config : ExperimentConfig
        Experiment configuration.
    output_dir : Path
        Directory containing results and tables.
    progress : ProgressManager, optional
        Progress manager for output control.
    """
    results_dir = output_dir / "results"
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    # Collect lambda_star samples from all periods
    all_lambda_star = []
    for period in config.time_periods.keys():
        lambda_path = results_dir / f"mmcp_period_{period}" / "lambda_star.npy"
        if lambda_path.exists():
            samples = np.load(lambda_path)
            all_lambda_star.extend(samples)

    if not all_lambda_star:
        if progress:
            progress.info("Warning: No lambda_star samples found. Skipping Table 5.3.")
        return

    arr = np.array(all_lambda_star)
    stats = {
        "統計量": ["平均", "中央値", "95%CI下限", "95%CI上限"],
        "値": [
            f"{np.mean(arr):.4f}",
            f"{np.median(arr):.4f}",
            f"{np.percentile(arr, 2.5):.4f}",
            f"{np.percentile(arr, 97.5):.4f}",
        ],
    }

    output_path = tables_dir / "table_5_3_lambda_posterior.csv"
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["統計量", "値"])
        for stat, val in zip(stats["統計量"], stats["値"]):
            writer.writerow([stat, val])

    if progress:
        progress.info(f"Lambda posterior table saved to {output_path}")


def _get_boundary_coords(preprocessor) -> np.ndarray:
    """Get boundary coordinates from preprocessor.

    Boundary = cells where elevation is null AND not sea.
    Based on Notebook 28 implementation.

    Parameters
    ----------
    preprocessor : ObsidianDataPreprocessor
        Data preprocessor.

    Returns
    -------
    np.ndarray
        Boundary coordinates (n_points, 2).
    """
    import polars as pl

    boundary = (
        preprocessor.df_elevation.filter(
            pl.col("average_elevation").is_null(), ~pl.col("is_sea")
        )
        .select(["x", "y"])
        .to_numpy()
    )
    return boundary


def generate_figures(
    config: ExperimentConfig,
    output_dir: Path,
    progress: ProgressManager,
) -> None:
    """Generate all figures for Chapter 5.

    Uses saved .npy results and recreates datasets for metadata.

    Parameters
    ----------
    config : ExperimentConfig
        Experiment configuration.
    output_dir : Path
        Directory containing results and figures.
    progress : ProgressManager
        Progress manager for output control.
    """
    from bayesian_statistics.nngp.model.marked_point_process import (
        prepare_marked_point_process_dataset,
    )

    from .visualization import (
        MapPlotter,
        apply_japanese_font,
        plot_all_origins_all_periods,
        plot_distance_prior,
        plot_effect_by_periods_and_origins,
        plot_estimated_vs_observed,
        plot_lambda_star_diagnostics,
        plot_loocv_comparison,
        plot_model_comparison_map,
        plot_site_probability_map,
        plot_uncertainty_map,
    )

    progress.section("Generating figures")
    apply_japanese_font()

    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    results_dir = output_dir / "results"

    if not results_dir.exists():
        progress.info("Error: No results found. Run --run-models first.")
        return

    # 1. Load preprocessor
    preprocessor = load_preprocessor(config)

    # 2. Load saved results and recreate datasets
    all_grid_probs: Dict[int, np.ndarray] = {}
    all_effects: Dict[int, Dict[str, np.ndarray]] = {}
    all_datasets: Dict[int, Any] = {}
    all_lambda_star: Dict[int, np.ndarray] = {}

    for period in config.time_periods.keys():
        period_dir = results_dir / f"mmcp_period_{period}"

        if not period_dir.exists():
            progress.info(f"Warning: Results for period {period} not found, skipping.")
            continue

        # Load saved numpy arrays
        all_grid_probs[period] = np.load(period_dir / "grid_probs.npy")
        all_lambda_star[period] = np.load(period_dir / "lambda_star.npy")
        all_effects[period] = {
            "distance": np.load(period_dir / "effect_distance.npy"),
            "intercept": np.load(period_dir / "effect_intercept.npy"),
            "intercept_adjustment": np.load(
                period_dir / "effect_intercept_adjustment.npy"
            ),
            "full": np.load(period_dir / "effect_full.npy"),
        }

        # Recreate dataset for visualization metadata
        all_datasets[period] = prepare_marked_point_process_dataset(
            preprocessor=preprocessor,
            period=period,
            origins=config.origins,
            grid_subsample_ratio=config.grid_subsample_ratio,
            distance_column_names=config.distance_column_names,
            source_weights=config.source_weights,
            tau=config.tau,
            alpha=config.alpha,
        )

    if not all_datasets:
        progress.info("Error: No valid results found.")
        return

    # 3. Create MapPlotter with boundary and land mask
    first_dataset = list(all_datasets.values())[0]
    boundary = _get_boundary_coords(preprocessor)

    plotter = MapPlotter(
        boundary=boundary,
        is_land=first_dataset.valid_grids,
    )
    grid_coords = first_dataset.grid_coords

    # 4. Generate Figure 5.1: All origins x all periods (4x5 grid)
    progress.info("Generating Figure 5.1: All origins x all periods")
    plot_all_origins_all_periods(
        plotter=plotter,
        grid_coords=grid_coords,
        all_grid_probs=all_grid_probs,
        all_datasets=all_datasets,
        origins=config.origins,
        time_periods=config.time_periods,
        output_path=figures_dir / "fig_5_1_all_origins_periods.png",
    )

    # 5. Generate Figure 5.2: Distance effect decomposition
    progress.info("Generating Figure 5.2: Effect decomposition (distance)")
    plot_effect_by_periods_and_origins(
        plotter=plotter,
        grid_coords=grid_coords,
        all_effects=all_effects,
        all_datasets=all_datasets,
        effect_name="distance",
        periods=list(config.time_periods.keys()),
        origin_indices=[0, 1, 2, 3],  # First 4 origins
        origins=config.origins,
        time_periods=config.time_periods,
        output_path=figures_dir / "fig_5_2_effect_distance.png",
    )

    # 5b. Generate Figure 5.2b: Intercept adjustment effect
    progress.info("Generating Figure 5.2b: Intercept adjustment (data correction)")
    plot_effect_by_periods_and_origins(
        plotter=plotter,
        grid_coords=grid_coords,
        all_effects=all_effects,
        all_datasets=all_datasets,
        effect_name="intercept_adjustment",
        periods=list(config.time_periods.keys()),
        origin_indices=[0, 1, 2, 3],
        origins=config.origins,
        time_periods=config.time_periods,
        scatter=True,
        vcenter=0.2,
        output_path=figures_dir / "fig_5_2b_effect_intercept_adjustment.png",
    )

    # 6. Generate Figure 5.3: Estimated vs Observed scatter
    progress.info("Generating Figure 5.3: Estimated vs Observed")
    all_estimated = []
    all_observed = []
    for period, dataset in all_datasets.items():
        site_probs = np.load(results_dir / f"mmcp_period_{period}" / "site_probs.npy")
        counts = dataset.counts
        total = counts.sum(axis=1, keepdims=True)
        observed = np.where(total > 0, counts / total, 0).T  # (K, n_sites)
        all_estimated.append(site_probs)
        all_observed.append(observed)

    estimated_concat = np.concatenate(all_estimated, axis=1)
    observed_concat = np.concatenate(all_observed, axis=1)

    plot_estimated_vs_observed(
        estimated=estimated_concat,
        observed=observed_concat,
        origins=config.origins,
        output_path=figures_dir / "fig_5_3_estimated_vs_observed.png",
    )

    # 7. Generate Figure 5.4: Lambda star diagnostics
    progress.info("Generating Figure 5.4: Lambda star diagnostics")
    first_period = list(all_lambda_star.keys())[0]
    plot_lambda_star_diagnostics(
        lambda_star_samples=all_lambda_star[first_period],
        output_path=figures_dir / "fig_5_4_lambda_diagnostics.png",
    )

    # 8. Generate Figure 5.5: Distance prior p0
    progress.info("Generating Figure 5.5: Distance prior p0")
    first_dataset = all_datasets[first_period]
    p0_grid = _compute_p0_from_dataset(first_dataset)

    if p0_grid is not None:
        plot_distance_prior(
            plotter=plotter,
            grid_coords=grid_coords,
            p0_grid=p0_grid,
            origins=config.origins,
            output_path=figures_dir / "fig_5_5_distance_prior.png",
        )
    else:
        progress.info("Skipping Figure 5.5: No distance features in dataset")

    # 9. Generate Figure 5.6: Model comparison (NW vs MMCP)
    # Check if NW results exist
    nw_period_dir = results_dir / f"nw_period_{first_period}"
    if nw_period_dir.exists():
        nw_grid_probs = np.load(nw_period_dir / "grid_probs.npy")

        # Check if grid sizes match (NW might use full grid, MMCP uses subsampled)
        if _check_grid_sizes_match(all_grid_probs[first_period], nw_grid_probs):
            progress.info("Generating Figure 5.6: Model comparison (NW vs MMCP)")
            first_dataset = all_datasets[first_period]
            true_ratio = first_dataset.counts / first_dataset.counts.sum(
                axis=1, keepdims=True
            )
            true_ratio = np.nan_to_num(true_ratio).T

            plot_model_comparison_map(
                plotter=plotter,
                grid_coords=grid_coords,
                mmcp_probs=all_grid_probs[first_period],
                nw_probs=nw_grid_probs,
                site_coords=first_dataset.site_coords,
                true_ratio=true_ratio,
                origins=config.origins,
                period_name=config.time_periods[first_period],
                output_path=figures_dir / "fig_5_6_model_comparison.png",
            )
        else:
            progress.info(
                "Skipping Figure 5.6: Grid sizes don't match "
                f"(MMCP: {all_grid_probs[first_period].shape[1]}, "
                f"NW: {nw_grid_probs.shape[1]})"
            )
    else:
        progress.info("Skipping Figure 5.6: NW results not found")

    # 10. Generate Figure 5.7: Uncertainty map
    progress.info("Generating Figure 5.7: Uncertainty map")
    # Load posterior standard deviation if available
    first_period_dir = results_dir / f"mmcp_period_{first_period}"
    grid_probs_std_path = first_period_dir / "grid_probs_std.npy"

    if grid_probs_std_path.exists():
        posterior_std = np.load(grid_probs_std_path)
    else:
        progress.info("  Warning: grid_probs_std.npy not found, using approximation")
        first_effects = all_effects[first_period]
        posterior_std = np.zeros_like(all_grid_probs[first_period])
        for k in range(4):
            adj = first_effects["intercept_adjustment"][k]
            posterior_std[k] = np.abs(adj - np.mean(adj[first_dataset.valid_grids]))

    plot_uncertainty_map(
        plotter=plotter,
        grid_coords=grid_coords,
        posterior_std=posterior_std,
        site_coords=first_dataset.site_coords,
        origins=config.origins,
        period_name=config.time_periods[first_period],
        output_path=figures_dir / "fig_5_7_uncertainty.png",
    )

    # 11. Generate Figure 5.8: LOOCV comparison
    # Check if LOOCV results exist
    loocv_table_path = output_dir / "tables" / "table_5_2_loocv.csv"
    if loocv_table_path.exists():
        progress.info("Generating Figure 5.8: LOOCV comparison")
        import csv

        loocv_results: Dict[str, Dict[int, Dict[str, float]]] = {"mmcp": {}, "nw": {}}
        with open(loocv_table_path, encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row_idx, row in enumerate(reader):
                period = row_idx
                try:
                    loocv_results["mmcp"][period] = {
                        "mean_aitchison_distance": float(row["MMCP平均"]),
                        "std_aitchison_distance": float(row["MMCP標準偏差"]),
                    }
                    loocv_results["nw"][period] = {
                        "mean_aitchison_distance": float(row["NW平均"]),
                        "std_aitchison_distance": float(row["NW標準偏差"]),
                    }
                except (ValueError, KeyError):
                    pass

        if loocv_results["mmcp"] or loocv_results["nw"]:
            plot_loocv_comparison(
                loocv_results=loocv_results,
                time_periods=config.time_periods,
                output_path=figures_dir / "fig_5_8_loocv_comparison.png",
            )
    else:
        progress.info("Skipping Figure 5.8: LOOCV results not found")

    # 12. Generate Figure 5.9: Site existence probability map
    progress.info("Generating Figure 5.9: Site existence probability")
    site_prob_path = first_period_dir / "site_probability.npy"

    if site_prob_path.exists():
        site_probability = np.load(site_prob_path)
        first_dataset = all_datasets[first_period]

        plot_site_probability_map(
            plotter=plotter,
            grid_coords=grid_coords,
            site_probability=site_probability,
            site_coords=first_dataset.site_coords,
            period_name=config.time_periods[first_period],
            output_path=figures_dir / "fig_5_9_site_probability.png",
        )
    else:
        progress.info("Skipping Figure 5.9: site_probability.npy not found")

    progress.info(f"\nAll figures saved to {figures_dir}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Run experiments for master thesis Chapter 5"
    )
    parser.add_argument(
        "--run-models",
        action="store_true",
        help="Run MMCP and NW models for all periods",
    )
    parser.add_argument(
        "--mmcp-only",
        action="store_true",
        help="Run only MMCP model (skip NW)",
    )
    parser.add_argument(
        "--run-loocv",
        action="store_true",
        help="Run LOOCV evaluation for both models",
    )
    parser.add_argument(
        "--generate-figures",
        action="store_true",
        help="Generate all figures",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run models, LOOCV, and generate figures",
    )
    parser.add_argument(
        "--n-iter",
        type=int,
        default=1000,
        help="Number of MCMC iterations (default: 1000)",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Quiet mode: suppress progress bars and most output",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Verbose mode: show all output including nested progress bars",
    )

    args = parser.parse_args()

    # Default: show help
    if not any(
        [args.run_models, args.mmcp_only, args.run_loocv, args.generate_figures, args.all]
    ):
        parser.print_help()
        return

    # Determine verbosity level
    if args.quiet:
        verbosity = ProgressManager.QUIET
    elif args.verbose:
        verbosity = ProgressManager.VERBOSE
    else:
        verbosity = ProgressManager.NORMAL

    # Create progress manager and set globally
    progress = ProgressManager(verbosity=verbosity)
    set_progress_manager(progress)

    # Create config
    config = ExperimentConfig(n_iter=args.n_iter, verbosity=verbosity)

    # Auto-adjust for testing if n_iter is small
    if args.n_iter < config.burn_in + config.thinning * 10:
        progress.info(
            f"Note: n_iter={args.n_iter} is small. "
            f"Auto-adjusting burn_in/thinning for testing."
        )
        config.adjust_for_testing(args.n_iter)
        progress.info(
            f"  Adjusted: burn_in={config.burn_in}, "
            f"thinning={config.thinning}, n_saved={config.n_saved()}"
        )

    # Validate configuration
    warnings = config.validate()
    for warning in warnings:
        progress.info(f"WARNING: {warning}")

    if config.n_saved() == 0 and (
        args.run_models or args.mmcp_only or args.run_loocv or args.all
    ):
        progress.info("ERROR: n_saved=0. Cannot run models without saving samples.")
        progress.info(f"  Current: n_iter={config.n_iter}, burn_in={config.burn_in}")
        progress.info(f"  Fix: Use --n-iter {config.burn_in + config.thinning * 10} or higher")
        return

    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    preprocessor = None
    if args.run_models or args.mmcp_only or args.run_loocv or args.all:
        preprocessor = load_preprocessor(config)

    if args.run_models or args.mmcp_only or args.all:
        run_mmcp_all_periods(config, preprocessor, output_dir, progress)
        if not args.mmcp_only:
            run_nw_all_periods(config, preprocessor, output_dir, progress)

        # Generate summary tables after model runs
        save_data_summary_table(preprocessor, config, output_dir, progress)
        save_lambda_posterior_table(config, output_dir, progress)

    if args.run_loocv or args.all:
        run_loocv_evaluation(config, preprocessor, output_dir, progress)

    if args.generate_figures or args.all:
        generate_figures(config, output_dir, progress)

    progress.section("All experiments completed!")


if __name__ == "__main__":
    main()
