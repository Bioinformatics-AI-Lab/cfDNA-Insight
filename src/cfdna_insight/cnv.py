"""Coverage-based copy-number features.

This is deliberately *not* a reimplementation of ichorCNA. There is no HMM, no
tumour-fraction estimate, no purity/ploidy solution. What it does is compute
GC-corrected, healthy-referenced coverage z-scores at 5 Mb resolution and
summarise them per chromosome arm, which is enough to give the classifier
large-scale aneuploidy signal alongside the fragmentomics.

Two summary statistics are worth more than the raw arm profile:

  - `cnv_arm_absz_mean`: mean |z| across arms, a scalar aneuploidy burden
  - `cnv_n_arms_altered`: arms exceeding |z| > 3

Both are close in spirit to what a PA (plasma aneuploidy) score measures, and
both behave sensibly at the 1-2x coverage typical of these datasets.

The healthy reference must be built from *training-fold healthy samples only*.
Building it from all samples leaks label information into every feature; see
`modeling.build_reference_within_fold`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .bins import aggregate_to_windows
from .fragmentomics import gc_correct


def coverage_windows(binned: pd.DataFrame, window_factor: int = 50,
                     gc: np.ndarray | None = None) -> pd.DataFrame:
    """5 Mb GC-corrected, depth-normalised coverage."""
    b = binned.copy()
    if gc is not None:
        b["n_total"] = gc_correct(b["n_total"].to_numpy(), gc)
    win = aggregate_to_windows(b, factor=window_factor, value_cols=["n_total"])
    win = win[win["frac_centromeric"] < 0.5].copy()
    total = win["n_total"].sum()
    win["cov"] = win["n_total"] / total if total > 0 else np.nan
    return win


class HealthyReference:
    """Per-window median and MAD coverage across a set of healthy samples."""

    def __init__(self, median: pd.Series, mad: pd.Series):
        self.median = median
        self.mad = mad

    @classmethod
    def fit(cls, coverage_matrix: pd.DataFrame) -> "HealthyReference":
        """coverage_matrix: samples x windows of normalised coverage."""
        med = coverage_matrix.median(axis=0)
        mad = (coverage_matrix - med).abs().median(axis=0) * 1.4826
        # windows with zero dispersion would give infinite z; floor the MAD
        floor = mad[mad > 0].median() if (mad > 0).any() else 1e-9
        mad = mad.clip(lower=floor * 0.1)
        return cls(med, mad)

    def zscore(self, coverage: pd.Series) -> pd.Series:
        common = self.median.index.intersection(coverage.index)
        return (coverage[common] - self.median[common]) / self.mad[common]


def arm_summary(z: pd.Series, window_arms: pd.Series) -> pd.Series:
    """Collapse window z-scores to per-arm medians plus scalar burden stats."""
    df = pd.DataFrame({"z": z, "arm": window_arms.reindex(z.index)})
    df = df[np.isfinite(df["z"])]
    if df.empty:
        return pd.Series(dtype=float)

    arm_z = df.groupby("arm")["z"].median()
    # per-arm z of a median of k windows shrinks by ~sqrt(k); rescale so that
    # the |z| > 3 threshold means the same thing regardless of arm length
    n_per_arm = df.groupby("arm")["z"].size()
    arm_z = arm_z * np.sqrt(n_per_arm)

    out = arm_z.add_prefix("cnv_arm_")
    out["cnv_arm_absz_mean"] = float(arm_z.abs().mean())
    out["cnv_arm_absz_max"] = float(arm_z.abs().max())
    out["cnv_n_arms_altered"] = float((arm_z.abs() > 3).sum())
    out["cnv_window_mad"] = float(
        (df["z"] - df["z"].median()).abs().median() * 1.4826)
    return out


def technical_covariates(frags: pd.DataFrame,
                         binned: pd.DataFrame) -> pd.Series:
    """Covariates that describe the *library*, not the biology.

    These are extracted so they can be used as an explicit negative control:
    a model trained on these alone should not separate cancer from healthy.
    If it does, the main model's AUC is measuring batch, not disease.
    """
    n = len(frags)
    cov = binned["n_total"].to_numpy(dtype=float)
    nonzero = cov[cov > 0]
    return pd.Series({
        "tech_n_fragments": float(n),
        "tech_log10_n_fragments": float(np.log10(max(n, 1))),
        "tech_median_mapq": float(frags["mapq"].median())
        if "mapq" in frags.columns else np.nan,
        "tech_frac_bins_zero": float(np.mean(cov == 0)),
        "tech_coverage_cv": float(nonzero.std() / nonzero.mean())
        if nonzero.size and nonzero.mean() > 0 else np.nan,
        "tech_chr_entropy": _chrom_entropy(frags),
    })


def _chrom_entropy(frags: pd.DataFrame) -> float:
    counts = frags["chrom"].value_counts(normalize=True)
    p = counts[counts > 0].to_numpy()
    return float(-(p * np.log(p)).sum() / np.log(len(p))) if p.size > 1 else np.nan
