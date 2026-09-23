"""Difference-of-Gaussians hotspot detection for irregularly spaced genes.

The detector treats gene-level log rates as a one-dimensional signal on a
circular chromosome. It subtracts a broad Gaussian smoother (background) from
a narrow Gaussian smoother (local signal), then propagates NPE uncertainty by
applying the filter to posterior draws. Distances are shortest-path distances
around the chromosome, so genes near the coordinate origin borrow information
from genes near the chromosome end.
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
    chromosome_length_bp: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a circular, uncertainty-weighted Gaussian smoothing matrix.

    Distances use the shorter of the clockwise and anticlockwise paths:
    d_circular(x, y) = min(|x-y|, L-|x-y|).
    """
    linear_distance = np.abs(
        target_positions[:, None] - source_positions[None, :]
    )
    distance = np.minimum(
        linear_distance, chromosome_length_bp - linear_distance
    )
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
    chromosome_length_bp: float
    min_fold_enrichment: float
    probability_threshold: float
    sort_order: np.ndarray

    @property
    def local_fold_enrichment_median(self) -> np.ndarray:
        """Median narrow-scale/broad-scale fold enrichment."""
        return np.exp(self.dog_log_enrichment_median)

    def segments(self, max_gap_bp: Optional[float] = None) -> list[dict]:
        """Return called regions, merging a region that crosses the origin.

        Consecutive hotspot-called genes are grouped while their circular
        inter-gene distance does not exceed ``max_gap_bp``. If the final and
        first genes are both called and are sufficiently close across the
        chromosome origin, their two linear runs are returned as one wrapped
        segment with ``wraps_origin=True``.
        """
        if max_gap_bp is None:
            max_gap_bp = 2.0 * self.narrow_bandwidth_bp
        if max_gap_bp <= 0:
            raise ValueError("max_gap_bp must be positive")

        called = self.hotspot_call
        n = len(called)
        if n == 0 or not np.any(called):
            return []

        # First construct ordinary runs in sorted coordinate order.
        runs: list[list[int]] = []
        current: list[int] = []
        for i in range(n):
            if not called[i]:
                if current:
                    runs.append(current)
                    current = []
                continue

            if current:
                gap = self.positions_bp[i] - self.positions_bp[current[-1]]
                if gap > max_gap_bp:
                    runs.append(current)
                    current = []
            current.append(i)

        if current:
            runs.append(current)

        # The last and first genes are adjacent on a circular chromosome.
        boundary_gap = (
            self.positions_bp[0]
            + self.chromosome_length_bp
            - self.positions_bp[-1]
        )
        if (
            len(runs) >= 2
            and called[0]
            and called[-1]
            and boundary_gap <= max_gap_bp
        ):
            wrapped = runs[-1] + runs[0]
            runs = [wrapped] + runs[1:-1]

        segments: list[dict] = []
        for indices in runs:
            idx = np.asarray(indices, dtype=int)
            wraps_origin = bool(np.any(np.diff(idx) < 0))
            start_index = int(idx[0])
            end_index = int(idx[-1])
            segments.append(
                {
                    "start_index": start_index,
                    "end_index": end_index,
                    "start_position_bp": float(
                        self.positions_bp[start_index]
                    ),
                    "end_position_bp": float(self.positions_bp[end_index]),
                    "n_genes": int(len(idx)),
                    "wraps_origin": wraps_origin,
                    "peak_probability": float(
                        self.hotspot_probability[idx].max()
                    ),
                    "median_probability": float(
                        np.median(self.hotspot_probability[idx])
                    ),
                    "peak_fold_enrichment": float(
                        self.local_fold_enrichment_median[idx].max()
                    ),
                }
            )
        return segments

    def plot(
        self,
        rate_label: str = "Recombination rate",
        position_unit_bp: float = 1e6,
        position_unit_name: str = "Mb",
        max_gap_bp: Optional[float] = None,
        rate_axis_log: bool = False,
        figsize: tuple[float, float] = (14.0, 9.0),
        save_path: Optional[str] = None,
        format: str = "pdf",
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

        # Shade called regions on all panels, including regions crossing
        # the circular chromosome origin.
        chromosome_end = self.chromosome_length_bp / position_unit_bp
        for segment in self.segments(max_gap_bp=max_gap_bp):
            i0, i1 = segment["start_index"], segment["end_index"]
            if segment["wraps_origin"]:
                left = 0.5 * (x[i0 - 1] + x[i0])
                right = 0.5 * (x[i1] + x[i1 + 1])
                spans = [(left, chromosome_end), (0.0, right)]
            else:
                left = x[i0] if i0 == 0 else 0.5 * (x[i0 - 1] + x[i0])
                right = (
                    x[i1]
                    if i1 == len(x) - 1
                    else 0.5 * (x[i1] + x[i1 + 1])
                )
                spans = [(left, right)]

            for span_left, span_right in spans:
                for ax in axes:
                    ax.axvspan(
                        span_left,
                        span_right,
                        color="#d62728",
                        alpha=0.10,
                        lw=0,
                    )

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
        if rate_axis_log:
            ax_rate.set_yscale("log")
        ax_rate.set_ylabel(rate_label)
        ax_rate.set_title("Difference-of-Gaussians Hotspot Detection")
        ax_rate.grid(alpha=0.2)
        ax_rate.legend(
            frameon=False,
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            borderaxespad=0.0,
            fontsize=9,
        )

        fold_median = np.exp(self.dog_log_enrichment_median)
        fold_lower = np.exp(self.dog_log_enrichment_lower)
        fold_upper = np.exp(self.dog_log_enrichment_upper)
        ax_dog.fill_between(
            x,
            fold_lower,
            fold_upper,
            color="#7f7f7f",
            alpha=0.22,
            label="95% interval",
        )
        ax_dog.plot(
            x,
            fold_median,
            color="#7b3294",
            linewidth=1.6,
            label="Median fold enrichment",
        )
        ax_dog.axhline(
            1.0,
            color="#555555",
            linestyle=":",
            linewidth=1.0,
            label="No enrichment",
        )
        ax_dog.axhline(
            self.min_fold_enrichment,
            color="#d62728",
            linestyle="--",
            linewidth=1.2,
            label=f"Minimum enrichment = {self.min_fold_enrichment:g}x",
        )
        ax_dog.set_yscale("linear")
        ax_dog.set_ylabel("Local / broad\nfold enrichment")
        ax_dog.grid(alpha=0.2)
        ax_dog.legend(
            frameon=False,
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            borderaxespad=0.0,
            fontsize=9,
        )

        ax_prob.plot(
            x,
            self.hotspot_probability,
            color="#222222",
            linewidth=1.4,
        )
        ax_prob.scatter(
            x,
            self.hotspot_probability,
            c=colors,
            s=18,
            zorder=3,
        )
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
            label=f"Threshold = {self.probability_threshold:g}",
        )
        ax_prob.set_ylim(-0.03, 1.03)
        ax_prob.set_ylabel("P(enrichment\nexceeds threshold)")
        ax_prob.set_xlabel(f"Chromosome position ({position_unit_name})")
        ax_prob.grid(alpha=0.2)
        ax_prob.legend(
            frameon=False,
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            borderaxespad=0.0,
            fontsize=9,
        )
        ax_prob.set_xlim(0.0, chromosome_end)

        # Reserve the right side of the figure for subplot-specific legends.
        fig.subplots_adjust(
            left=0.10,
            right=0.76,
            bottom=0.08,
            top=0.94,
            hspace=0.08,
        )
        fig.align_ylabels(axes)

        if save_path is not None:
            fig.savefig(save_path, bbox_inches="tight", format=format)

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
    chromosome_length_bp
        Total circular chromosome length. This is required to calculate the
        shortest distance across the coordinate origin.
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
        chromosome_length_bp: float,
        min_fold_enrichment: float = 1.5,
        probability_threshold: float = 0.90,
        residual_log_sd: float = 0.20,
        min_effective_genes: float = 1.0,
        max_posterior_draws: Optional[int] = 5000,
        random_state: Optional[int] = 0,
    ):
        if chromosome_length_bp <= 0:
            raise ValueError("chromosome_length_bp must be positive")
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
        self.chromosome_length_bp = float(chromosome_length_bp)
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
        if np.any(positions < 0) or np.any(
            positions > self.chromosome_length_bp
        ):
            raise ValueError(
                "positions_bp must lie in [0, chromosome_length_bp]"
            )
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
            positions,
            positions,
            self.narrow_bandwidth_bp,
            precision,
            self.chromosome_length_bp,
        )
        broad_matrix, broad_effective_n = _weighted_gaussian_matrix(
            positions,
            positions,
            self.broad_bandwidth_bp,
            precision,
            self.chromosome_length_bp,
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
            chromosome_length_bp=self.chromosome_length_bp,
            min_fold_enrichment=self.min_fold_enrichment,
            probability_threshold=self.probability_threshold,
            sort_order=order,
        )

