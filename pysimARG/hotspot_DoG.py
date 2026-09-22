"""Difference-of-Gaussians hotspot detection for irregularly spaced genes.

The detector treats gene-level log rates as a one-dimensional genomic signal.
It subtracts a broad Gaussian smoother (background) from a narrow Gaussian
smoother (local signal), then propagates NPE uncertainty by applying the
filter to posterior draws.

This is an exploratory posterior-stability method, not a formally coherent
chromosome-level Bayesian model. Calibrate detection thresholds using
end-to-end chromosomes simulated under a no-hotspot model.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import ArrayLike


_EPS = 1e-12


def _weighted_gaussian_matrix(
    source_positions: np.ndarray,
    target_positions: np.ndarray,
    bandwidth_bp: float,
    precision: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return normalized uncertainty-weighted Gaussian smoothing matrix."""
    distance = target_positions[:, None] - source_positions[None, :]
    kernel = np.exp(-0.5 * (distance / bandwidth_bp) ** 2)
    weighted = kernel * precision[None, :]
    support = weighted.sum(axis=1)
    matrix = weighted / np.maximum(support[:, None], _EPS)

    # Kish effective sample size, useful for diagnosing unsupported locations.
    effective_n = support**2 / np.maximum((weighted**2).sum(axis=1), _EPS)
    return matrix, effective_n


@dataclass
class DoGHotspotResult:
    """Results from uncertainty-aware Difference-of-Gaussians filtering."""

    positions_bp: np.ndarray
    posterior_rate_median: np.ndarray
    posterior_rate_lower: np.ndarray
    posterior_rate_upper: np.ndarray
    narrow_rate_median: np.ndarray
    broad_rate_median: np.ndarray
    dog_log_enrichment_median: np.ndarray
    dog_log_enrichment_lower: np.ndarray
    dog_log_enrichment_upper: np.ndarray
    hotspot_probability: np.ndarray
    hotspot_call: np.ndarray
    narrow_effective_n: np.ndarray
    broad_effective_n: np.ndarray
    narrow_bandwidth_bp: float
    broad_bandwidth_bp: float
    min_fold_enrichment: float
    probability_threshold: float
    sort_order: np.ndarray

    @property
    def local_fold_enrichment_median(self) -> np.ndarray:
        """Median narrow-scale/broad-scale fold enrichment."""
        return np.exp(self.dog_log_enrichment_median)

    def segments(self, max_gap_bp: Optional[float] = None) -> list[dict]:
        """Return contiguous called regions, optionally split across large gaps."""
        if max_gap_bp is None:
            max_gap_bp = 2.0 * self.narrow_bandwidth_bp
        if max_gap_bp <= 0:
            raise ValueError("max_gap_bp must be positive")

        segments: list[dict] = []
        start: Optional[int] = None
        n = len(self.positions_bp)

        for i in range(n):
            starts_new = self.hotspot_call[i] and (
                i == 0
                or not self.hotspot_call[i - 1]
                or self.positions_bp[i] - self.positions_bp[i - 1] > max_gap_bp
            )
            if starts_new:
                start = i

            ends_here = self.hotspot_call[i] and (
                i == n - 1
                or not self.hotspot_call[i + 1]
                or self.positions_bp[i + 1] - self.positions_bp[i] > max_gap_bp
            )
            if ends_here and start is not None:
                sl = slice(start, i + 1)
                segments.append(
                    {
                        "start_index": start,
                        "end_index": i,
                        "start_position_bp": float(self.positions_bp[start]),
                        "end_position_bp": float(self.positions_bp[i]),
                        "n_genes": i - start + 1,
                        "peak_probability": float(self.hotspot_probability[sl].max()),
                        "median_probability": float(np.median(self.hotspot_probability[sl])),
                        "peak_fold_enrichment": float(
                            self.local_fold_enrichment_median[sl].max()
                        ),
                    }
                )
                start = None
        return segments

    def plot(
        self,
        rate_label: str = "Recombination rate",
        position_unit_bp: float = 1e6,
        position_unit_name: str = "Mb",
        max_gap_bp: Optional[float] = None,
        figsize: tuple[float, float] = (12.0, 9.0),
        save_path: Optional[str] = None,
        dpi: int = 300,
    ):
        """Plot gene rates, DoG enrichment, and hotspot probability."""
        import matplotlib.pyplot as plt

        if position_unit_bp <= 0:
            raise ValueError("position_unit_bp must be positive")

        x = self.positions_bp / position_unit_bp
        fig, axes = plt.subplots(
            3,
            1,
            figsize=figsize,
            sharex=True,
            gridspec_kw={"height_ratios": [3.0, 2.1, 1.8], "hspace": 0.08},
        )
        ax_rate, ax_dog, ax_prob = axes

        # Shade called regions on all panels.
        for segment in self.segments(max_gap_bp=max_gap_bp):
            i0, i1 = segment["start_index"], segment["end_index"]
            left = x[i0] if i0 == 0 else 0.5 * (x[i0 - 1] + x[i0])
            right = x[i1] if i1 == len(x) - 1 else 0.5 * (x[i1] + x[i1 + 1])
            for ax in axes:
                ax.axvspan(left, right, color="#d62728", alpha=0.10, lw=0)

        colors = np.where(self.hotspot_call, "#d62728", "#4c78a8")
        ax_rate.vlines(
            x,
            self.posterior_rate_lower,
            self.posterior_rate_upper,
            color=colors,
            alpha=0.28,
            linewidth=0.8,
        )
        ax_rate.scatter(
            x,
            self.posterior_rate_median,
            c=colors,
            s=18,
            edgecolor="white",
            linewidth=0.3,
            zorder=3,
            label="Gene posterior median",
        )
        ax_rate.plot(
            x,
            self.narrow_rate_median,
            color="#e68613",
            linewidth=1.7,
            label=f"Narrow smoother ({self.narrow_bandwidth_bp / 1e3:g} kb)",
        )
        ax_rate.plot(
            x,
            self.broad_rate_median,
            color="#222222",
            linewidth=1.7,
            label=f"Broad background ({self.broad_bandwidth_bp / 1e3:g} kb)",
        )
        ax_rate.set_yscale("linear")
        ax_rate.set_ylabel(rate_label)
        ax_rate.set_title("Uncertainty-aware Difference-of-Gaussians hotspot detection")
        ax_rate.grid(alpha=0.2)
        ax_rate.legend(frameon=False, ncol=3, fontsize=9)

        fold_median = np.exp(self.dog_log_enrichment_median)
        fold_lower = np.exp(self.dog_log_enrichment_lower)
        fold_upper = np.exp(self.dog_log_enrichment_upper)
        ax_dog.fill_between(
            x, fold_lower, fold_upper, color="#7f7f7f", alpha=0.22, label="95% interval"
        )
        ax_dog.plot(x, fold_median, color="#7b3294", linewidth=1.6)
        ax_dog.axhline(1.0, color="#555555", linestyle=":", linewidth=1.0)
        ax_dog.axhline(
            self.min_fold_enrichment,
            color="#d62728",
            linestyle="--",
            linewidth=1.2,
            label=f"Minimum enrichment = {self.min_fold_enrichment:g}x",
        )
        ax_dog.set_yscale("log")
        ax_dog.set_ylabel("Local / broad\nfold enrichment")
        ax_dog.grid(alpha=0.2)
        ax_dog.legend(frameon=False, loc="upper right")

        ax_prob.plot(x, self.hotspot_probability, color="#222222", linewidth=1.4)
        ax_prob.scatter(x, self.hotspot_probability, c=colors, s=18, zorder=3)
        ax_prob.fill_between(
            x,
            self.probability_threshold,
            self.hotspot_probability,
            where=self.hotspot_probability >= self.probability_threshold,
            color="#d62728",
            alpha=0.28,
            interpolate=True,
        )
        ax_prob.axhline(
            self.probability_threshold,
            color="#d62728",
            linestyle="--",
            linewidth=1.2,
            label=f"Calling threshold = {self.probability_threshold:g}",
        )
        ax_prob.set_ylim(-0.03, 1.03)
        ax_prob.set_ylabel("P(enrichment\nexceeds threshold)")
        ax_prob.set_xlabel(f"Chromosome position ({position_unit_name})")
        ax_prob.grid(alpha=0.2)
        ax_prob.legend(frameon=False, loc="upper right")

        fig.subplots_adjust(left=0.10, right=0.98, bottom=0.08, top=0.94, hspace=0.08)
        if save_path is not None:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        return fig, axes


class DoGHotspotDetector:
    """Detect local rate elevations relative to a broad genomic background.

    Parameters
    ----------
    narrow_bandwidth_bp
        Gaussian bandwidth for the local signal. Choose this near the
        expected hotspot radius, not its full width.
    broad_bandwidth_bp
        Gaussian bandwidth for the background. It must be larger than the
        narrow bandwidth, commonly by a factor of 4 to 10.
    min_fold_enrichment
        A location is considered enriched in a posterior draw when the
        narrow smoother exceeds the broad smoother by at least this factor.
    probability_threshold
        Minimum posterior-draw frequency used to call a hotspot.
    residual_log_sd
        Extra log-rate variability added to each gene's posterior variance
        when constructing smoothing weights. This prevents very precise genes
        from dominating their neighborhoods.
    min_effective_genes
        Locations with fewer than this effective number of genes in either
        smoother are not called, although their scores remain available.
    max_posterior_draws
        Optional random subsample size for computational efficiency.
    random_state
        Seed used only when posterior draws are subsampled.
    """

    def __init__(
        self,
        narrow_bandwidth_bp: float,
        broad_bandwidth_bp: float,
        min_fold_enrichment: float = 1.5,
        probability_threshold: float = 0.90,
        residual_log_sd: float = 0.20,
        min_effective_genes: float = 2.0,
        max_posterior_draws: Optional[int] = 5000,
        random_state: Optional[int] = 0,
    ):
        if narrow_bandwidth_bp <= 0:
            raise ValueError("narrow_bandwidth_bp must be positive")
        if broad_bandwidth_bp <= narrow_bandwidth_bp:
            raise ValueError("broad_bandwidth_bp must exceed narrow_bandwidth_bp")
        if min_fold_enrichment <= 1:
            raise ValueError("min_fold_enrichment must exceed 1")
        if not 0 < probability_threshold < 1:
            raise ValueError("probability_threshold must be between 0 and 1")
        if residual_log_sd < 0:
            raise ValueError("residual_log_sd cannot be negative")
        if min_effective_genes <= 0:
            raise ValueError("min_effective_genes must be positive")
        if max_posterior_draws is not None and max_posterior_draws < 2:
            raise ValueError("max_posterior_draws must be at least 2")

        self.narrow_bandwidth_bp = float(narrow_bandwidth_bp)
        self.broad_bandwidth_bp = float(broad_bandwidth_bp)
        self.min_fold_enrichment = float(min_fold_enrichment)
        self.probability_threshold = float(probability_threshold)
        self.residual_log_sd = float(residual_log_sd)
        self.min_effective_genes = float(min_effective_genes)
        self.max_posterior_draws = max_posterior_draws
        self.random_state = random_state

    def fit(
        self,
        positions_bp: ArrayLike,
        posterior_draws: ArrayLike,
        rates_are_log: bool = False,
    ) -> DoGHotspotResult:
        """Apply the DoG filter to gene-level posterior draws.

        Parameters
        ----------
        positions_bp
            One genomic coordinate per gene, typically gene midpoint.
        posterior_draws
            Array with shape (n_draws, n_genes). Pass joint NPE draws for one
            parameter at a time, preserving each gene's posterior uncertainty.
        rates_are_log
            Set True if posterior_draws already contain natural-log rates.

        Returns
        -------
        DoGHotspotResult
            Arrays are sorted by genomic position.
        """
        positions = np.asarray(positions_bp, dtype=float)
        draws = np.asarray(posterior_draws, dtype=float)

        if positions.ndim != 1:
            raise ValueError("positions_bp must be one-dimensional")
        if draws.ndim != 2:
            raise ValueError("posterior_draws must have shape (n_draws, n_genes)")
        if draws.shape[1] != len(positions):
            raise ValueError("posterior_draws columns must match positions_bp")
        if draws.shape[0] < 2:
            raise ValueError("At least two posterior draws are required")
        if not np.all(np.isfinite(positions)) or not np.all(np.isfinite(draws)):
            raise ValueError("positions_bp and posterior_draws must be finite")
        if len(np.unique(positions)) != len(positions):
            raise ValueError("Gene positions must be unique")
        if not rates_are_log and np.any(draws <= 0):
            raise ValueError("Rates must be positive when rates_are_log=False")

        order = np.argsort(positions)
        positions = positions[order]
        draws = draws[:, order]

        if (
            self.max_posterior_draws is not None
            and draws.shape[0] > self.max_posterior_draws
        ):
            rng = np.random.default_rng(self.random_state)
            selected = rng.choice(
                draws.shape[0], size=self.max_posterior_draws, replace=False
            )
            draws = draws[selected]

        log_draws = draws if rates_are_log else np.log(draws)
        posterior_log_variance = np.var(log_draws, axis=0, ddof=1)
        precision = 1.0 / (
            posterior_log_variance + self.residual_log_sd**2 + _EPS
        )

        narrow_matrix, narrow_effective_n = _weighted_gaussian_matrix(
            positions, positions, self.narrow_bandwidth_bp, precision
        )
        broad_matrix, broad_effective_n = _weighted_gaussian_matrix(
            positions, positions, self.broad_bandwidth_bp, precision
        )

        narrow_log_draws = log_draws @ narrow_matrix.T
        broad_log_draws = log_draws @ broad_matrix.T
        dog_draws = narrow_log_draws - broad_log_draws

        threshold = np.log(self.min_fold_enrichment)
        hotspot_probability = np.mean(dog_draws > threshold, axis=0)
        adequate_support = (
            (narrow_effective_n >= self.min_effective_genes)
            & (broad_effective_n >= self.min_effective_genes)
        )
        hotspot_call = (
            hotspot_probability >= self.probability_threshold
        ) & adequate_support

        rate_draws = np.exp(log_draws)
        rate_q = np.quantile(rate_draws, [0.025, 0.5, 0.975], axis=0)
        dog_q = np.quantile(dog_draws, [0.025, 0.5, 0.975], axis=0)

        return DoGHotspotResult(
            positions_bp=positions,
            posterior_rate_median=rate_q[1],
            posterior_rate_lower=rate_q[0],
            posterior_rate_upper=rate_q[2],
            narrow_rate_median=np.exp(np.median(narrow_log_draws, axis=0)),
            broad_rate_median=np.exp(np.median(broad_log_draws, axis=0)),
            dog_log_enrichment_median=dog_q[1],
            dog_log_enrichment_lower=dog_q[0],
            dog_log_enrichment_upper=dog_q[2],
            hotspot_probability=hotspot_probability,
            hotspot_call=hotspot_call,
            narrow_effective_n=narrow_effective_n,
            broad_effective_n=broad_effective_n,
            narrow_bandwidth_bp=self.narrow_bandwidth_bp,
            broad_bandwidth_bp=self.broad_bandwidth_bp,
            min_fold_enrichment=self.min_fold_enrichment,
            probability_threshold=self.probability_threshold,
            sort_order=order,
        )


def simulate_example(
    n_genes: int = 260,
    n_posterior_draws: int = 3000,
    chromosome_length_bp: float = 30e6,
    random_state: int = 7,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Simulate irregular gene positions and uncertain gene-level rate posteriors."""
    rng = np.random.default_rng(random_state)
    positions = np.sort(rng.uniform(0, chromosome_length_bp, n_genes))

    # Slowly varying background plus two local elevations of different widths.
    broad_background = (
        np.log(1.0)
        + 0.18 * np.sin(2.0 * np.pi * positions / chromosome_length_bp)
        + 0.10 * np.cos(4.0 * np.pi * positions / chromosome_length_bp)
    )
    hotspot_1 = np.log(3.2) * np.exp(
        -0.5 * ((positions - 8.5e6) / 0.40e6) ** 2
    )
    hotspot_2 = np.log(2.3) * np.exp(
        -0.5 * ((positions - 21.5e6) / 0.85e6) ** 2
    )
    true_log_rate = broad_background + hotspot_1 + hotspot_2

    # Mimic heterogeneous NPE uncertainty and slight gene-level estimation error.
    posterior_log_sd = rng.uniform(0.12, 0.38, n_genes)
    posterior_center = true_log_rate + rng.normal(0.0, 0.10, n_genes)
    log_draws = rng.normal(
        posterior_center[None, :],
        posterior_log_sd[None, :],
        size=(n_posterior_draws, n_genes),
    )
    return positions, np.exp(log_draws), np.exp(true_log_rate)

