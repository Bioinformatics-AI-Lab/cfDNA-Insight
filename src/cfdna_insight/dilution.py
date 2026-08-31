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
            binned = fx.bin_short_long(mixed, bin_size=bin_size, genome=genome)

            feats = pd.concat([
                fx.size_summary(mixed),
                fx.size_distribution(mixed),
                fx.delfi_profile(binned, window_factor=window_factor, gc=gc),
            ])
            X = feats.reindex(feature_columns).to_frame().T
            score = float(model.predict_proba(X)[:, 1][0])
            rows.append({"fraction": frac, "replicate": rep, "score": score})

    return pd.DataFrame(rows)


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
