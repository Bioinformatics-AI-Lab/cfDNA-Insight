"""In silico dilution series.

The question a screening assay is actually judged on is not "can you separate
late-stage cancer from healthy" but "at what tumour fraction does the signal
disappear". Fragment files make that directly testable without any new
sequencing: build synthetic samples by mixing fragments from a cancer sample
and a healthy sample at a known ratio, hold total depth constant, recompute
features, and score them with a model that never saw the mixture.

Caveats that belong in the write-up, not buried here:

  - The mixing fraction is a *fragment* fraction, not a tumour fraction. The
    cancer donor's plasma is itself mostly haematopoietic cfDNA, so a 10%
    mixture of a sample with 5% ctDNA is roughly 0.5% ctDNA. The x-axis is
    "fraction of cancer-patient plasma", and the resulting LOD is an upper
    bound on the true ctDNA LOD.
  - Mixing is done at the fragment level with replacement, so at high
    depths the same fragment can be drawn twice. This under-represents the
    true diversity of a real library at that depth.
  - Real dilution also changes library complexity, duplication rate and error
    profile. None of that is captured.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import cnv, fragmentomics as fx
from .io import read_fragments


def mix_fragments(cancer: pd.DataFrame, healthy: pd.DataFrame,
                  fraction: float, total: int,
                  rng: np.random.Generator) -> pd.DataFrame:
    """Draw `total` fragments, `fraction` of them from the cancer sample."""
    n_cancer = int(round(total * fraction))
    n_healthy = total - n_cancer
    parts = []
    if n_cancer > 0:
        idx = rng.choice(len(cancer), n_cancer, replace=n_cancer > len(cancer))
        parts.append(cancer.iloc[idx])
    if n_healthy > 0:
        idx = rng.choice(len(healthy), n_healthy,
                         replace=n_healthy > len(healthy))
        parts.append(healthy.iloc[idx])
    mixed = pd.concat(parts, ignore_index=True)
    return mixed.sort_values(["chrom", "start"], ignore_index=True)


def dilution_series(
    cancer_path: str,
    healthy_path: str,
    model,
    feature_columns: list[str],
    fractions=(0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 1.0),
    total_fragments: int = 5_000_000,
    n_replicates: int = 5,
    gc_track: pd.DataFrame | None = None,
    bin_size: int = 100_000,
    window_factor: int = 50,
    genome: str = "hg38",
    seed: int = 0,
) -> pd.DataFrame:
    """Score synthetic mixtures across a range of cancer-plasma fractions."""
    rng = np.random.default_rng(seed)
    cancer = read_fragments(cancer_path)
    healthy = read_fragments(healthy_path)
    gc = gc_track["gc"].to_numpy() if gc_track is not None else None

    rows = []
    for frac in fractions:
        for rep in range(n_replicates):
            mixed = mix_fragments(cancer, healthy, frac, total_fragments, rng)
            feats = mixture_features(mixed, gc=gc, bin_size=bin_size,
                                     window_factor=window_factor,
                                     genome=genome,
                                     need=set(feature_columns))
            X = feats.reindex(feature_columns).to_frame().T
            score = float(model.predict_proba(X)[:, 1][0])
            row = {"fraction": frac, "replicate": rep, "score": score}
            # carry the raw aneuploidy scalars through: when the model is the
            # two-feature one, the curve in the feature itself is easier to
            # read than the curve in the model's output
            for k in ("an_chr_entropy", "an_coverage_cv"):
                if k in feats.index:
                    row[k] = float(feats[k])
            rows.append(row)

    return pd.DataFrame(rows)


def mixture_features(mixed: pd.DataFrame, gc=None, bin_size: int = 100_000,
                     window_factor: int = 50, genome: str = "hg38",
                     need: set[str] | None = None) -> pd.Series:
    """Every feature recomputable from a mixture without a reference panel.

    This is the constraint that decides what can appear in a dilution series
    at all. Arm-level CNV z-scores need a healthy reference panel, and a panel
    built from the same donors that produced the mixture would leak; the
    fragmentomic blocks and the two aneuploidy scalars need nothing but the
    fragments in hand.

    `an_chr_entropy` in particular is just the evenness of the per-chromosome
    fragment distribution, so it is exactly as meaningful on a synthetic
    mixture as on a real library.
    """
    binned = fx.bin_short_long(mixed, bin_size=bin_size, genome=genome)
    parts = [cnv.aneuploidy_covariates(mixed, binned), fx.size_summary(mixed)]
    # the 5 Mb DELFI profile is by far the most expensive piece, and a series
    # run on the aneuploidy block alone does not need it
    if need is None or any(str(n).startswith(("delfi_", "size_")) for n in need):
        parts.append(fx.size_distribution(mixed))
        parts.append(fx.delfi_profile(binned, window_factor=window_factor,
                                      gc=gc))
    return pd.concat(parts)


def estimate_lod(series: pd.DataFrame, threshold: float,
                 min_detection_rate: float = 0.95) -> float:
    """Lowest tested fraction at which detection rate meets the target.

    Reported as "LOD95 <= x" rather than "LOD = x": with a handful of
    replicates per level the estimate is coarse, and the answer is bounded by
    the grid of fractions actually tested.
    """
    rates = (series.assign(detected=series["score"] >= threshold)
             .groupby("fraction")["detected"].mean().sort_index())
    hits = rates[rates >= min_detection_rate]
    return float(hits.index.min()) if len(hits) else float("nan")


def threshold_at_specificity(scores_healthy: np.ndarray,
                             specificity: float = 0.95) -> float:
    """Score cutoff giving the requested specificity in known-healthy samples.

    Screening assays are specificity-anchored: at population prevalence, a
    95%-specific test generates far more false than true positives. The
    threshold is therefore set on healthy samples first, and sensitivity is
    read off afterwards, never the reverse.
    """
    return float(np.quantile(scores_healthy, specificity))


def dose_response(series: pd.DataFrame, feature: str,
                  n_sd: float = 3.0) -> dict:
    """Characterise one feature's response to the mixing fraction.

    A population-threshold LOD answers "would a screening assay flag this
    sample", and on a feature with large between-person variation the answer
    can be no even when the feature tracks tumour content perfectly. Two
    further numbers separate those cases:

    `r2` and `slope` say whether the feature is quantitative. A feature linear
    in mixing fraction is measuring tumour content, not merely correlated
    with disease status.

    `lod_own_baseline` is the smallest fraction distinguishable from the same
    donor's own zero point, at `n_sd` times the replicate standard deviation.
    That is the relevant limit when a prior sample from the same person
    exists -- monitoring and MRD -- and it can be far lower than the
    population-threshold LOD, because it never pays the between-person
    variance.

    Both assume the replicate spread seen here, which captures sampling noise
    in the mixture but not library prep, extraction or batch, so the
    own-baseline figure is optimistic for a real longitudinal assay.
    """
    if feature not in series.columns:
        raise KeyError(f"{feature} not in the series; columns: "
                       f"{list(series.columns)}")
    g = series.groupby("fraction")[feature]
    frac = np.array(sorted(series["fraction"].unique()), dtype=float)
    mean = g.mean().reindex(frac).to_numpy()
    sd = g.std(ddof=0).reindex(frac).to_numpy()

    slope, intercept = np.polyfit(frac, mean, 1)
    pred = intercept + slope * frac
    ss_res = float(((mean - pred) ** 2).sum())
    ss_tot = float(((mean - mean.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    sd0 = float(sd[0]) if np.isfinite(sd[0]) and sd[0] > 0 else float(np.nanmean(sd))
    span = float(mean[-1] - mean[0])
    lod_own = abs(n_sd * sd0 / slope) if slope else float("nan")

    return {
        "feature": feature,
        "slope_per_unit_fraction": float(slope),
        "direction": "increases with cancer fraction" if slope > 0
                     else "decreases with cancer fraction",
        "r2_linear": float(r2),
        "span_over_series": span,
        "replicate_sd_at_zero": sd0,
        "replicate_sd_as_pct_of_span": (abs(sd0 / span) * 100
                                        if span else float("nan")),
        "n_sd": n_sd,
        "lod_own_baseline_fraction": float(lod_own),
    }


def check_direction(dose: dict, cohort_case_higher: bool | None) -> str | None:
    """Flag a donor whose dose-response runs against the cohort direction.

    A dilution series is one donor pair, and a feature can be real at cohort
    level while this particular case sits on the wrong side of it. The series
    will still be beautifully linear, because mixing two samples that differ
    always produces a linear response in anything computed from their
    proportions; linearity says the feature resolves the mixture, not that it
    resolves cancer.

    `cohort_case_higher` is whether cases exceed controls for this feature
    across the cohort. When the donor's slope disagrees, the LOD is a
    statement about measurement resolution only, and quoting it as a
    detection limit would get the sign wrong.
    """
    if cohort_case_higher is None:
        return None
    donor_up = dose["slope_per_unit_fraction"] > 0
    if donor_up == cohort_case_higher:
        return None
    return (f"{dose['feature']}: this donor's case is on the opposite side of "
            f"the cohort trend (slope {dose['slope_per_unit_fraction']:+.3e}, "
            f"cohort has cases {'higher' if cohort_case_higher else 'lower'}). "
            f"The LOD measures resolution, not detection, for this pair.")
