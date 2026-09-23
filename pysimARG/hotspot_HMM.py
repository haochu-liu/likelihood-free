"""Distance-aware two-state HMM using NPE moments on the raw-rate scale.

Model
-----
For ordered genes g = 1, ..., G:

    z_g in {0=background, 1=hotspot}
    y_g | z_g=k ~ Normal(mu_k, sigma_k^2 + v_g)

where y_g and v_g are the posterior mean and variance of the gene's NPE
posterior for the rate itself (not the log rate). State changes follow a
continuous-time Markov chain along genomic distance:

    Q = [[-alpha, alpha],
         [ beta, -beta ]]

and the transition matrix over distance d is exp(Q*d).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import norm


_EPS = 1e-12


def _softplus(x: float | np.ndarray) -> float | np.ndarray:
    return np.logaddexp(0.0, x)


def _inv_softplus(y: float | np.ndarray) -> float | np.ndarray:
    y = np.maximum(y, 1e-10)
    return y + np.log(-np.expm1(-y))


@dataclass
class HotspotHMMResult:
    """Fitted raw-rate HMM results, sorted by genomic position."""

    positions: np.ndarray
    rate_mean: np.ndarray
    rate_variance: np.ndarray
    hotspot_probability: np.ndarray
    viterbi_state: np.ndarray
    background_mean: float
    hotspot_mean: float
    background_sd: float
    hotspot_sd: float
    background_to_hotspot_rate: float
    hotspot_to_background_rate: float
    expected_background_length: float
    expected_hotspot_length: float
    log_likelihood: float
    distance_unit_bp: float
    optimizer_success: bool
    optimizer_message: str
    rate_scale: float

    @property
    def fold_enrichment(self) -> float:
        """Estimated hotspot/background ratio of raw-rate state means."""
        return float(self.hotspot_mean / self.background_mean)

    def hotspot_calls(self, probability_threshold: float = 0.9) -> np.ndarray:
        if not 0.0 <= probability_threshold <= 1.0:
            raise ValueError("probability_threshold must be between 0 and 1")
        return self.hotspot_probability >= probability_threshold

    def segments(self, probability_threshold: float = 0.9) -> list[dict]:
        """Return contiguous called segments; boundaries are gene positions."""
        called = self.hotspot_calls(probability_threshold)
        segments: list[dict] = []
        start: Optional[int] = None

        for i, is_hot in enumerate(called):
            if is_hot and start is None:
                start = i
            if start is not None and (not is_hot or i == len(called) - 1):
                end = i if is_hot and i == len(called) - 1 else i - 1
                p = self.hotspot_probability[start : end + 1]
                segments.append(
                    {
                        "start_index": start,
                        "end_index": end,
                        "start_position": float(self.positions[start]),
                        "end_position": float(self.positions[end]),
                        "n_genes": end - start + 1,
                        "mean_hotspot_probability": float(p.mean()),
                        "min_hotspot_probability": float(p.min()),
                    }
                )
                start = None
        return segments

    def plot_states(
        self,
        position_unit_bp: float = 1e6,
        position_unit_name: str = "Mb",
        rate_label: str = "Recombination rate",
        probability_threshold: float = None,
        figsize: tuple[float, float] = (12.0, 8.0),
        save_path: Optional[str] = None,
        format: str = "pdf",
    ):
        """Plot raw-rate posterior moments, hotspot probabilities, and states."""
        if probability_threshold is not None and not 0.0 <= probability_threshold <= 1.0:
            raise ValueError("probability_threshold must be between 0 and 1")
        if position_unit_bp <= 0:
            raise ValueError("position_unit_bp must be positive")

        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        x = self.positions / position_unit_bp
        state = self.viterbi_state.astype(int)
        p_hot = self.hotspot_probability
        rate_sd = np.sqrt(self.rate_variance)
        lower = np.maximum(self.rate_mean - 1.96 * rate_sd, np.finfo(float).tiny)
        upper = self.rate_mean + 1.96 * rate_sd

        colors = np.where(state == 1, "#d62728", "#4c78a8")
        fig, axes = plt.subplots(
            3,
            1,
            figsize=figsize,
            sharex=True,
            gridspec_kw={"height_ratios": [3.2, 2.0, 0.8], "hspace": 0.08},
        )
        ax_rate, ax_prob, ax_state = axes

        edges = np.empty(len(x) + 1)
        edges[1:-1] = (x[:-1] + x[1:]) / 2.0
        edges[0] = x[0] - (x[1] - x[0]) / 2.0
        edges[-1] = x[-1] + (x[-1] - x[-2]) / 2.0
        for i, is_hot in enumerate(state):
            if is_hot:
                for ax in axes:
                    ax.axvspan(edges[i], edges[i + 1], color="#d62728", alpha=0.10, lw=0)

        ax_rate.vlines(x, lower, upper, color=colors, alpha=0.38, linewidth=1.0)
        ax_rate.scatter(
            x,
            self.rate_mean,
            c=colors,
            s=25,
            edgecolor="white",
            linewidth=0.4,
            zorder=3,
        )
        ax_rate.axhline(
            self.background_mean,
            color="#4c78a8",
            linestyle="--",
            linewidth=1.4,
            label="Estimated background mean",
        )
        ax_rate.axhline(
            self.hotspot_mean,
            color="#d62728",
            linestyle="--",
            linewidth=1.4,
            label="Estimated hotspot mean",
        )
        ax_rate.set_ylabel(rate_label)
        ax_rate.set_title(
            f"Hidden Markov Model Hotspot Detection"
            # f"({self.fold_enrichment:.2f}x estimated enrichment)"
        )
        ax_rate.grid(alpha=0.2)
        ax_rate.legend(frameon=False, loc="best")

        ax_prob.plot(x, p_hot, color="#222222", linewidth=1.2)
        ax_prob.scatter(x, p_hot, c=colors, s=22, zorder=3)
        if probability_threshold is not None:
            ax_prob.fill_between(
                x,
                probability_threshold,
                p_hot,
                where=p_hot >= probability_threshold,
                color="#d62728",
                alpha=0.25,
                interpolate=True,
            )
            ax_prob.axhline(
                probability_threshold,
                color="#d62728",
                linestyle="--",
                linewidth=1.1,
                label=f"Calling threshold = {probability_threshold:.2f}",
            )
        ax_prob.set_ylim(-0.03, 1.03)
        ax_prob.set_ylabel("P(hotspot)")
        ax_prob.grid(alpha=0.2)
        if probability_threshold is not None:
            ax_prob.legend(frameon=False, loc="best")

        ax_state.step(x, state, where="mid", color="#333333", linewidth=1.3)
        ax_state.scatter(x, state, c=colors, s=22, zorder=3)
        ax_state.set_yticks([0, 1], labels=["Background", "Hotspot"])
        ax_state.set_ylim(-0.45, 1.45)
        ax_state.set_ylabel("Viterbi")
        ax_state.set_xlabel(f"Chromosome position ({position_unit_name})")
        ax_state.grid(axis="x", alpha=0.2)

        state_legend = [
            Line2D(
                [0], [0], marker="o", color="none", markerfacecolor="#4c78a8",
                markeredgecolor="white", markersize=7, label="Background gene"
            ),
            Line2D(
                [0], [0], marker="o", color="none", markerfacecolor="#d62728",
                markeredgecolor="white", markersize=7, label="Hotspot gene"
            ),
        ]
        ax_state.legend(handles=state_legend, frameon=False, loc="upper right", ncol=2)

        fig.align_ylabels(axes)
        if save_path is not None:
            fig.savefig(save_path, bbox_inches="tight", format=format)
        return fig, axes


class DistanceAwareHotspotHMM:
    """Two-state continuous-distance HMM with raw-rate Gaussian emissions.

    Parameters
    ----------
    distance_unit_bp
        Unit used for transition rates. With 1e6, alpha and beta are rates
        per Mb, and expected segment lengths are reported in base pairs.
    min_fold_enrichment
        Minimum ratio between the hotspot and background raw-rate means.
        A value of 1.0 only enforces hotspot_mean > background_mean.
    n_starts
        Number of numerical optimization restarts.
    random_state
        Seed controlling randomized optimizer starting values.
    """

    def __init__(
        self,
        distance_unit_bp: float = 1e6,
        min_fold_enrichment: float = 1.0,
        n_starts: int = 8,
        random_state: Optional[int] = 1,
    ) -> None:
        if distance_unit_bp <= 0:
            raise ValueError("distance_unit_bp must be positive")
        if min_fold_enrichment < 1.0:
            raise ValueError("min_fold_enrichment must be at least 1")
        if n_starts < 1:
            raise ValueError("n_starts must be at least 1")

        self.distance_unit_bp = float(distance_unit_bp)
        self.min_fold_enrichment = float(min_fold_enrichment)
        self.n_starts = int(n_starts)
        self.rng = np.random.default_rng(random_state)
        self.result_: Optional[HotspotHMMResult] = None

    @staticmethod
    def summarize_npe_samples(
        posterior_samples: ArrayLike,
        variance_floor: float = 0.0,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return per-gene means and variances of raw-rate NPE samples.

        posterior_samples must have shape (n_genes, n_draws). If joint NPE
        samples contain recombination and mutation rates, pass one parameter
        slice at a time, for example samples[:, :, 0].
        """
        samples = np.asarray(posterior_samples, dtype=float)
        if samples.ndim != 2:
            raise ValueError("posterior_samples must have shape (genes, draws)")
        if samples.shape[1] < 2:
            raise ValueError("At least two posterior draws per gene are required")
        if not np.all(np.isfinite(samples)):
            raise ValueError("posterior_samples contains non-finite values")
        if np.any(samples < 0):
            raise ValueError("Rate samples cannot be negative")
        if variance_floor < 0:
            raise ValueError("variance_floor cannot be negative")

        means = samples.mean(axis=1)
        variances = samples.var(axis=1, ddof=1)
        return means, np.maximum(variances, variance_floor)

    @staticmethod
    def _validate_and_sort(
        positions: ArrayLike,
        rate_mean: ArrayLike,
        rate_variance: ArrayLike,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        x = np.asarray(positions, dtype=float).reshape(-1)
        y = np.asarray(rate_mean, dtype=float).reshape(-1)
        v = np.asarray(rate_variance, dtype=float).reshape(-1)

        if not (len(x) == len(y) == len(v)):
            raise ValueError("positions, means, and variances must have equal length")
        if len(x) < 3:
            raise ValueError("At least three genes are required")
        if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y)) and np.all(np.isfinite(v))):
            raise ValueError("Inputs contain non-finite values")
        if np.any(y <= 0):
            raise ValueError("Posterior mean rates must be positive")
        if np.any(v < 0):
            raise ValueError("Posterior variances cannot be negative")

        order = np.argsort(x)
        x, y, v = x[order], y[order], v[order]
        if np.any(np.diff(x) <= 0):
            raise ValueError("Gene positions must be unique")
        return x, y, v

    @staticmethod
    def _transition_matrix(distance: float, alpha: float, beta: float) -> np.ndarray:
        """Closed-form exp(Q*d) for a two-state continuous-time chain."""
        total = alpha + beta
        stationary_hot = alpha / total
        stationary_bg = beta / total
        decay = np.exp(-total * distance)
        return np.array(
            [
                [stationary_bg + stationary_hot * decay, stationary_hot * (1.0 - decay)],
                [stationary_bg * (1.0 - decay), stationary_hot + stationary_bg * decay],
            ]
        )

    def _unpack(
        self, raw: np.ndarray
    ) -> tuple[float, float, float, float, float, float]:
        # Parameters here are on the internally scaled raw-rate axis.
        mu0 = _softplus(raw[0]) + 1e-10
        mu1 = self.min_fold_enrichment * mu0 + _softplus(raw[1])
        sigma0 = np.exp(raw[2])
        sigma1 = np.exp(raw[3])
        alpha = np.exp(raw[4])
        beta = np.exp(raw[5])
        return float(mu0), float(mu1), float(sigma0), float(sigma1), float(alpha), float(beta)

    def _components(
        self,
        raw: np.ndarray,
        y: np.ndarray,
        v: np.ndarray,
        distances: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        mu0, mu1, sigma0, sigma1, alpha, beta = self._unpack(raw)
        initial = np.array([beta / (alpha + beta), alpha / (alpha + beta)])

        log_emission = np.column_stack(
            [
                norm.logpdf(y, loc=mu0, scale=np.sqrt(sigma0**2 + v)),
                norm.logpdf(y, loc=mu1, scale=np.sqrt(sigma1**2 + v)),
            ]
        )
        transitions = np.stack(
            [self._transition_matrix(d, alpha, beta) for d in distances], axis=0
        )
        return (
            np.log(np.maximum(initial, _EPS)),
            log_emission,
            np.log(np.maximum(transitions, _EPS)),
        )

    @staticmethod
    def _forward(
        log_initial: np.ndarray,
        log_emission: np.ndarray,
        log_transition: np.ndarray,
    ) -> tuple[float, np.ndarray]:
        n = len(log_emission)
        forward = np.empty((n, 2))
        forward[0] = log_initial + log_emission[0]
        for g in range(1, n):
            forward[g] = log_emission[g] + logsumexp(
                forward[g - 1][:, None] + log_transition[g - 1], axis=0
            )
        return float(logsumexp(forward[-1])), forward

    @staticmethod
    def _backward(log_emission: np.ndarray, log_transition: np.ndarray) -> np.ndarray:
        n = len(log_emission)
        backward = np.zeros((n, 2))
        for g in range(n - 2, -1, -1):
            backward[g] = logsumexp(
                log_transition[g]
                + log_emission[g + 1][None, :]
                + backward[g + 1][None, :],
                axis=1,
            )
        return backward

    @staticmethod
    def _viterbi(
        log_initial: np.ndarray,
        log_emission: np.ndarray,
        log_transition: np.ndarray,
    ) -> np.ndarray:
        n = len(log_emission)
        score = np.empty((n, 2))
        pointer = np.zeros((n, 2), dtype=int)
        score[0] = log_initial + log_emission[0]

        for g in range(1, n):
            candidate = score[g - 1][:, None] + log_transition[g - 1]
            pointer[g] = np.argmax(candidate, axis=0)
            score[g] = log_emission[g] + candidate[pointer[g], np.arange(2)]

        states = np.empty(n, dtype=int)
        states[-1] = int(np.argmax(score[-1]))
        for g in range(n - 2, -1, -1):
            states[g] = pointer[g + 1, states[g + 1]]
        return states

    def _initial_parameters(self, y: np.ndarray, distances: np.ndarray) -> np.ndarray:
        q30, q80 = np.quantile(y, [0.30, 0.80])
        spread = max(float(np.std(y)), 0.05)

        mu0_init = max(float(q30), 1e-4)
        minimum_hotspot_mean = self.min_fold_enrichment * mu0_init
        mu1_init = max(float(q80), minimum_hotspot_mean + 0.05 * spread)
        excess_hotspot_mean = max(mu1_init - minimum_hotspot_mean, 1e-4)

        median_distance = max(float(np.median(distances)), 1e-4)
        expected_bg = max(10.0 * median_distance, 0.1)
        expected_hot = max(3.0 * median_distance, 0.05)

        return np.array(
            [
                _inv_softplus(mu0_init),
                _inv_softplus(excess_hotspot_mean),
                np.log(max(0.5 * spread, 1e-4)),
                np.log(max(0.7 * spread, 1e-4)),
                np.log(1.0 / expected_bg),
                np.log(1.0 / expected_hot),
            ],
            dtype=float,
        )

    def fit(
        self,
        positions: ArrayLike,
        posterior_samples: Optional[ArrayLike] = None,
        *,
        rate_mean: Optional[ArrayLike] = None,
        rate_variance: Optional[ArrayLike] = None,
    ) -> HotspotHMMResult:
        """Fit the raw-rate HMM and return smoothed hotspot probabilities.

        Supply either posterior_samples with shape (genes, draws), or both
        rate_mean and rate_variance. Positions should be in base pairs.
        Returned arrays are sorted by position.
        """
        if posterior_samples is not None:
            if rate_mean is not None or rate_variance is not None:
                raise ValueError("Use posterior_samples or summary moments, not both")
            y, v = self.summarize_npe_samples(posterior_samples)
        else:
            if rate_mean is None or rate_variance is None:
                raise ValueError("Supply posterior_samples or both posterior moments")
            y = np.asarray(rate_mean, dtype=float)
            v = np.asarray(rate_variance, dtype=float)

        x, y, v = self._validate_and_sort(positions, y, v)
        distances = np.diff(x) / self.distance_unit_bp

        # Optimize on a dimensionless scale to remain stable for rates such as 1e-8.
        rate_scale = float(np.median(y[y > 0]))
        if not np.isfinite(rate_scale) or rate_scale <= 0:
            rate_scale = float(np.max(y))
        y_scaled = y / rate_scale
        v_scaled = np.maximum(v / rate_scale**2, 1e-12)

        initial_raw = self._initial_parameters(y_scaled, distances)

        def objective(raw: np.ndarray) -> float:
            try:
                log_pi, log_e, log_t = self._components(
                    raw, y_scaled, v_scaled, distances
                )
                log_likelihood, _ = self._forward(log_pi, log_e, log_t)
                if not np.isfinite(log_likelihood):
                    return 1e100
                return -log_likelihood
            except (FloatingPointError, ValueError, OverflowError):
                return 1e100

        candidates = []
        for start in range(self.n_starts):
            raw0 = initial_raw.copy()
            if start:
                raw0 += self.rng.normal(
                    0.0, [0.35, 0.5, 0.35, 0.35, 0.8, 0.8]
                )
            opt = minimize(
                objective,
                raw0,
                method="L-BFGS-B",
                options={"maxiter": 3000},
            )
            candidates.append(opt)

        best = min(candidates, key=lambda z: z.fun)
        mu0_s, mu1_s, sigma0_s, sigma1_s, alpha, beta = self._unpack(best.x)
        log_pi, log_e, log_t = self._components(
            best.x, y_scaled, v_scaled, distances
        )
        log_likelihood_scaled, forward = self._forward(log_pi, log_e, log_t)
        backward = self._backward(log_e, log_t)
        log_gamma = forward + backward - log_likelihood_scaled
        gamma = np.exp(log_gamma - logsumexp(log_gamma, axis=1, keepdims=True))
        states = self._viterbi(log_pi, log_e, log_t)

        # The scaled density differs by 1/rate_scale for each gene.
        log_likelihood_original_units = (
            log_likelihood_scaled - len(y) * np.log(rate_scale)
        )

        self.result_ = HotspotHMMResult(
            positions=x,
            rate_mean=y,
            rate_variance=v,
            hotspot_probability=gamma[:, 1],
            viterbi_state=states,
            background_mean=float(mu0_s * rate_scale),
            hotspot_mean=float(mu1_s * rate_scale),
            background_sd=float(sigma0_s * rate_scale),
            hotspot_sd=float(sigma1_s * rate_scale),
            background_to_hotspot_rate=float(alpha),
            hotspot_to_background_rate=float(beta),
            expected_background_length=float(self.distance_unit_bp / alpha),
            expected_hotspot_length=float(self.distance_unit_bp / beta),
            log_likelihood=float(log_likelihood_original_units),
            distance_unit_bp=self.distance_unit_bp,
            optimizer_success=bool(best.success),
            optimizer_message=str(best.message),
            rate_scale=rate_scale,
        )
        return self.result_

