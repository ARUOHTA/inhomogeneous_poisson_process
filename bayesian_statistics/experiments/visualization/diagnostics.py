"""MCMC diagnostics visualization."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np


def plot_lambda_star_diagnostics(
    lambda_star_samples: np.ndarray,
    output_path: Optional[Union[str, Path]] = None,
) -> tuple:
    """Plot lambda star trace plot and posterior distribution.

    Parameters
    ----------
    lambda_star_samples : np.ndarray
        Array of lambda star samples.
    output_path : str or Path, optional
        Path to save the figure.

    Returns
    -------
    fig : Figure
        Matplotlib figure.
    axes : np.ndarray
        Array of matplotlib axes.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Check for empty or invalid data
    if lambda_star_samples.size == 0:
        for ax in axes:
            ax.text(
                0.5,
                0.5,
                "データなし\n(MCMCサンプルが保存されていません)",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=12,
            )
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
        axes[0].set_title("λ* のトレースプロット")
        axes[1].set_title("λ* の事後分布")
        plt.tight_layout()
        if output_path:
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
        return fig, axes

    # Trace plot
    ax = axes[0]
    ax.plot(lambda_star_samples)
    ax.set_xlabel("サンプル番号")
    ax.set_ylabel("λ*")
    ax.set_title("λ* のトレースプロット")

    # Posterior distribution
    ax = axes[1]
    ax.hist(lambda_star_samples, bins=30, density=True, alpha=0.7)
    ax.axvline(
        lambda_star_samples.mean(),
        color="red",
        linestyle="--",
        label=f"平均: {lambda_star_samples.mean():.2f}",
    )
    ax.set_xlabel("λ*")
    ax.set_ylabel("密度")
    ax.set_title("λ* の事後分布")
    ax.legend()

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")

    return fig, axes


def print_lambda_star_summary(lambda_star_samples: np.ndarray) -> None:
    """Print lambda star summary statistics.

    Parameters
    ----------
    lambda_star_samples : np.ndarray
        Array of lambda star samples.
    """
    print("=" * 50)
    print("λ* の事後サマリー")
    print("=" * 50)
    print(f"平均:     {lambda_star_samples.mean():.2f}")
    print(f"中央値:   {np.median(lambda_star_samples):.2f}")
    print(f"標準偏差: {lambda_star_samples.std():.2f}")
    print(
        f"95%CI:    [{np.percentile(lambda_star_samples, 2.5):.2f}, "
        f"{np.percentile(lambda_star_samples, 97.5):.2f}]"
    )


def print_pseudo_absence_summary(n_pseudo_absence: list) -> None:
    """Print pseudo absence statistics.

    Parameters
    ----------
    n_pseudo_absence : list
        List of pseudo absence counts per iteration.
    """
    n_U_array = np.array(n_pseudo_absence)
    print("=" * 50)
    print("偽不在点の統計")
    print("=" * 50)
    print(f"平均: {n_U_array.mean():.1f}")
    print(f"最小: {n_U_array.min()}")
    print(f"最大: {n_U_array.max()}")


def plot_beta_trace(
    beta_samples: np.ndarray,
    category_idx: int = 0,
    feature_idx: int = 0,
    site_idx: int = 0,
    output_path: Optional[Union[str, Path]] = None,
) -> tuple:
    """Plot beta coefficient trace plot.

    Parameters
    ----------
    beta_samples : np.ndarray
        Array of beta samples (n_samples, K-1, p+1, n_sites).
    category_idx : int
        Category index to plot.
    feature_idx : int
        Feature index to plot.
    site_idx : int
        Site index to plot.
    output_path : str or Path, optional
        Path to save the figure.

    Returns
    -------
    fig : Figure
        Matplotlib figure.
    ax : Axes
        Matplotlib axes.
    """
    fig, ax = plt.subplots(figsize=(10, 4))

    samples = beta_samples[:, category_idx, feature_idx, site_idx]
    ax.plot(samples)
    ax.set_xlabel("サンプル番号")
    ax.set_ylabel(f"β[{category_idx},{feature_idx},{site_idx}]")
    ax.set_title("βのトレースプロット")

    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")

    return fig, ax
