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
from sklearn.base import BaseEstimator, TransformerMixin

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
    """Covariates that describe the *library*, and nothing else.

    These exist to be a negative control: a model trained on these alone
    should not separate cancer from healthy. If it does, the main model's AUC
    is measuring batch rather than disease.

    For that to be a real test, the block must contain nothing biological.
    An earlier version of this function also returned chromosome-distribution
    entropy and bin-coverage dispersion. Both are aneuploidy readouts —
    per-chromosome fragment proportions and coverage variance are precisely
    what large-scale copy-number change perturbs — so including them put a
    cancer signal inside the negative control and made it fire by
    construction. On the CRAG cohort chromosome entropy alone reached AUC
    0.879, which inflated the "technical" block to 0.886 and made the design
    look confounded when it was not. They now live in
    `aneuploidy_covariates`.

    The rule for anything added here: if a tumour could change it, it is not
    a technical covariate.
    """
    n = len(frags)
    cov = binned["n_total"].to_numpy(dtype=float)
    return pd.Series({
        "tech_n_fragments": float(n),
        "tech_log10_n_fragments": float(np.log10(max(n, 1))),
        "tech_median_mapq": float(frags["mapq"].median())
        if "mapq" in frags.columns else np.nan,
        "tech_frac_bins_zero": float(np.mean(cov == 0)),
    })


def aneuploidy_covariates(frags: pd.DataFrame,
                          binned: pd.DataFrame) -> pd.Series:
    """Cheap whole-genome aneuploidy proxies, needing no healthy reference.

    Both describe how unevenly fragments are spread across the genome, which
    is what arm- and chromosome-scale copy-number change does to coverage.
    Unlike the arm-level z-scores in `arm_summary`, neither needs a reference
    panel, so neither can leak label information across CV folds.

    `an_chr_entropy` is the normalised Shannon entropy of the per-chromosome
    fragment proportions. In the CRAG cohort it separated cases from their
    own matched controls in 91% of 33 pairs (Wilcoxon p = 3e-07) — matched on
    age, sex, smoking and alcohol, and run in the same batch, so the paired
    design rules out those explanations.
    """
    cov = binned["n_total"].to_numpy(dtype=float)
    nonzero = cov[cov > 0]
    return pd.Series({
        "an_chr_entropy": _chrom_entropy(frags),
        "an_coverage_cv": float(nonzero.std() / nonzero.mean())
        if nonzero.size and nonzero.mean() > 0 else np.nan,
    })


def _chrom_entropy(frags: pd.DataFrame) -> float:
    counts = frags["chrom"].value_counts(normalize=True)
    p = counts[counts > 0].to_numpy()
    return float(-(p * np.log(p)).sum() / np.log(len(p))) if p.size > 1 else np.nan


# --------------------------------------------------------------------------
# Fold-safe CNV features
# --------------------------------------------------------------------------

class HealthyReferenceTransformer(BaseEstimator, TransformerMixin):
    """Build the healthy coverage reference inside the CV fold.

    `build_cnv_block` in 02_train_evaluate.py fits the reference once, over
    every healthy sample in the cohort. Each test-fold healthy sample
    therefore contributes to the median and MAD it is later scored against,
    and with only a few dozen healthy samples that contribution is not small:
    at n healthy, each one is roughly 1/n of its own reference. The effect is
    to pull test healthy samples toward z = 0 and inflate CNV performance.

    Wrapping the reference in a transformer hands the problem to sklearn.
    Pipeline.fit passes the training fold's X and y, so the reference is built
    from training healthy samples only, and the test fold is genuinely held
    out.

    Input X is the samples x windows coverage matrix. Output is the per-arm
    summary that `arm_summary` produces.
    """

    def __init__(self, window_arms: pd.Series | None = None,
                 altered_z: float = 3.0):
        self.window_arms = window_arms
        self.altered_z = altered_z

    def fit(self, X, y=None):
        X = pd.DataFrame(X)
        if y is None:
            raise ValueError("HealthyReferenceTransformer needs y to identify "
                             "healthy samples")
        y = np.asarray(y)
        healthy = X[y == 0]
        if len(healthy) < 3:
            raise ValueError(f"need >=3 healthy samples in the training fold "
                             f"to build a reference, got {len(healthy)}")
        self.reference_ = HealthyReference.fit(healthy)
        self.columns_ = list(X.columns)

        # Build the summaries for the whole training fold, not just one probe
        # row: an arm can be well defined for one sample and empty for the
        # next. Columns with no finite value anywhere in training are dropped
        # here rather than left for the imputer to skip with a warning.
        # In practice these are the acrocentric p-arms (chr13p, chr14p,
        # chr15p, chr21p, chr22p), which are unassembled N-sequence in hg19
        # and carry no fragments in any sample.
        summaries = pd.DataFrame(
            [arm_summary(self.reference_.zscore(X.iloc[i]), self.window_arms)
             for i in range(len(X))],
            index=X.index,
        )
        usable = summaries.columns[summaries.notna().any(axis=0)]
        self.dropped_ = [c for c in summaries.columns if c not in set(usable)]
        self.feature_names_ = list(usable)
        if not self.feature_names_:
            raise ValueError("no usable CNV arm features in this fold")
        return self

    def transform(self, X):
        X = pd.DataFrame(X)
        rows = [arm_summary(self.reference_.zscore(X.iloc[i]),
                            self.window_arms)
                for i in range(len(X))]
        out = pd.DataFrame(rows, index=X.index)
        return out.reindex(columns=self.feature_names_)

    def get_feature_names_out(self, input_features=None):
        return np.asarray(self.feature_names_, dtype=object)
