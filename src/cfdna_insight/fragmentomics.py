"""Fragmentomic feature extraction.

Three feature blocks, each capturing a different physical property:

1. Global size distribution
   Normal cfDNA is dominated by mononucleosomal fragments (~167 bp, with 10 bp
   periodicity from the helical pitch of DNA on the histone octamer). Tumour
   cfDNA is shifted shorter. This block is one histogram per sample.

2. DELFI short/long ratio profile (Cristiano et al. 2019)
   The ratio of short (100-150 bp) to long (151-220 bp) fragments in 5 Mb
   windows. What matters is the *shape* of this profile across the genome, not
   its mean, so profiles are z-scored within sample. That within-sample
   centring is deliberate: it removes global size shifts that also correlate
   with library prep and storage, leaving positional structure.

3. End-motif frequencies and motif diversity score (Jiang et al. 2020)
   The 4-mer at the 5' end of each fragment reflects the nuclease that cut it.
   Cancer plasma shows altered nuclease activity (notably DNASE1L3 loss), which
   shows up as a shifted 4-mer spectrum and reduced motif diversity.

Block 3 requires the reference FASTA; blocks 1 and 2 do not.
"""

from __future__ import annotations

from collections import Counter
from itertools import product

import numpy as np
import pandas as pd

from .bins import aggregate_to_windows, assign_bins, make_bins

SHORT_RANGE = (100, 150)
LONG_RANGE = (151, 220)
SIZE_HIST_RANGE = (50, 400)

FOUR_MERS = ["".join(p) for p in product("ACGT", repeat=4)]
_MER_INDEX = {m: i for i, m in enumerate(FOUR_MERS)}


# --------------------------------------------------------------------------
# Block 1: size distribution
# --------------------------------------------------------------------------

def size_distribution(frags: pd.DataFrame,
                      lo: int = SIZE_HIST_RANGE[0],
                      hi: int = SIZE_HIST_RANGE[1]) -> pd.Series:
    """Density-normalised fragment length histogram, 1 bp resolution."""
    lengths = frags["length"].to_numpy()
    counts = np.bincount(np.clip(lengths, lo, hi) - lo, minlength=hi - lo + 1)
    dens = counts / max(counts.sum(), 1)
    return pd.Series(dens, index=[f"size_{i}" for i in range(lo, hi + 1)])


def size_summary(frags: pd.DataFrame) -> pd.Series:
    """Scalar summaries of the size distribution."""
    L = frags["length"].to_numpy()
    if L.size == 0:
        return pd.Series(dtype=float)
    n_short = np.sum((L >= SHORT_RANGE[0]) & (L <= SHORT_RANGE[1]))
    n_long = np.sum((L >= LONG_RANGE[0]) & (L <= LONG_RANGE[1]))
    return pd.Series({
        "size_median": float(np.median(L)),
        "size_mean": float(L.mean()),
        "size_iqr": float(np.percentile(L, 75) - np.percentile(L, 25)),
        "frac_below_150": float(np.mean(L < 150)),
        "frac_below_120": float(np.mean(L < 120)),
        "frac_167_181": float(np.mean((L >= 167) & (L <= 181))),
        "frac_dinucleosomal": float(np.mean((L >= 275) & (L <= 400))),
        "short_long_ratio_global": float(n_short / n_long) if n_long else np.nan,
        "periodicity_10bp": _periodicity_score(L),
    })


def _periodicity_score(lengths: np.ndarray,
                       lo: int = 90, hi: int = 160) -> float:
    """Amplitude of the 10 bp oscillation in the size histogram.

    Computed as the magnitude of the DFT component nearest a 10 bp period,
    normalised by total spectral power. Loss of this periodicity indicates
    degraded or non-nucleosomal DNA and is a useful QC flag as much as a
    biological feature.
    """
    counts = np.bincount(np.clip(lengths, lo, hi) - lo, minlength=hi - lo + 1)
    if counts.sum() < 1000:
        return np.nan
    x = counts.astype(float)
    x = x - np.polyval(np.polyfit(np.arange(x.size), x, 3), np.arange(x.size))
    spec = np.abs(np.fft.rfft(x * np.hanning(x.size)))
    freqs = np.fft.rfftfreq(x.size, d=1.0)
    target = np.argmin(np.abs(freqs - 0.1))  # period 10 bp -> freq 0.1
    denom = spec[1:].sum()
    return float(spec[target] / denom) if denom > 0 else np.nan


# --------------------------------------------------------------------------
# Block 2: DELFI short/long profile
# --------------------------------------------------------------------------

def bin_short_long(frags: pd.DataFrame, bin_size: int = 100_000,
                   genome: str = "hg38") -> pd.DataFrame:
    """Count short, long and total fragments per 100 kb bin."""
    bins = make_bins(bin_size=bin_size, genome=genome)
    index = {bid: i for i, bid in enumerate(bins["bin_id"])}

    L = frags["length"].to_numpy()
    is_short = (L >= SHORT_RANGE[0]) & (L <= SHORT_RANGE[1])
    is_long = (L >= LONG_RANGE[0]) & (L <= LONG_RANGE[1])

    bin_idx = assign_bins(frags, bin_size)
    keys = frags["chrom"].astype(str).to_numpy() + ":" + \
        (bin_idx * bin_size).astype(str)
    pos = np.array([index.get(k, -1) for k in keys])
    valid = pos >= 0

    n_short = np.bincount(pos[valid & is_short], minlength=len(bins))
    n_long = np.bincount(pos[valid & is_long], minlength=len(bins))
    n_total = np.bincount(pos[valid], minlength=len(bins))

    out = bins.copy()
    out["n_short"] = n_short[: len(bins)]
    out["n_long"] = n_long[: len(bins)]
    out["n_total"] = n_total[: len(bins)]
    return out


def gc_correct(counts: np.ndarray, gc: np.ndarray,
               n_knots: int = 8, min_count: int = 1) -> np.ndarray:
    """GC-bias correction by piecewise-linear fit of count on GC content.

    A LOESS fit is standard here; a monotone-free spline over GC quantile knots
    behaves almost identically on binned counts and avoids a statsmodels
    dependency. Correction is multiplicative: corrected = observed * (median /
    predicted). Bins whose predicted count is degenerate are left untouched.
    """
    counts = counts.astype(float)
    ok = np.isfinite(gc) & (counts >= min_count)
    if ok.sum() < 100:
        return counts
    knots = np.quantile(gc[ok], np.linspace(0, 1, n_knots))
    knots = np.unique(knots)
    if knots.size < 3:
        return counts
    predicted = np.interp(
        gc, knots,
        [np.median(counts[ok & (np.abs(gc - k) < 0.02)]) if
         np.any(ok & (np.abs(gc - k) < 0.02)) else np.median(counts[ok])
         for k in knots],
    )
    med = np.median(counts[ok])
    with np.errstate(divide="ignore", invalid="ignore"):
        factor = np.where(predicted > 0, med / predicted, 1.0)
    return counts * factor


def delfi_profile(binned: pd.DataFrame, window_factor: int = 50,
                  gc: np.ndarray | None = None,
                  min_total: int = 200,
                  centre_scale: bool = True) -> pd.Series:
    """5 Mb short/long ratio profile, GC-corrected and z-scored within sample."""
    b = binned.copy()
    if gc is not None:
        b["n_short"] = gc_correct(b["n_short"].to_numpy(), gc)
        b["n_long"] = gc_correct(b["n_long"].to_numpy(), gc)
        b["n_total"] = gc_correct(b["n_total"].to_numpy(), gc)

    win = aggregate_to_windows(b, factor=window_factor,
                               value_cols=["n_short", "n_long", "n_total"])
    win = win[win["frac_centromeric"] < 0.5].copy()

    ratio = np.where(win["n_long"] >= min_total,
                     win["n_short"] / win["n_long"].replace(0, np.nan), np.nan)
    win["ratio"] = ratio

    vals = win["ratio"].to_numpy(dtype=float)
    if centre_scale:
        finite = np.isfinite(vals)
        if finite.sum() > 10:
            vals = (vals - np.nanmean(vals)) / (np.nanstd(vals) + 1e-9)
    return pd.Series(vals, index=["delfi_" + b for b in win["bin_id"]])


# --------------------------------------------------------------------------
# Block 3: end motifs
# --------------------------------------------------------------------------

def end_motifs(frags: pd.DataFrame, fasta, k: int = 4,
               max_fragments: int = 2_000_000,
               seed: int = 0) -> pd.Series:
    """Frequency of each k-mer at fragment 5' ends.

    For a fragment on the + strand the 5' end is the reference sequence at
    [start, start+k). For the - strand it is the reverse complement of
    [end-k, end). If the file has no strand column both ends are counted,
    which is the convention used when strand is unavailable.
    """
    rng = np.random.default_rng(seed)
    if len(frags) > max_fragments:
        frags = frags.iloc[rng.choice(len(frags), max_fragments, replace=False)]

    counts = Counter()
    has_strand = "strand" in frags.columns

    for chrom, sub in frags.groupby("chrom", observed=True):
        chrom = str(chrom)
        try:
            seq = fasta.fetch(chrom).upper()
        except (KeyError, ValueError):
            continue
        starts = sub["start"].to_numpy()
        ends = sub["end"].to_numpy()
        if has_strand:
            strands = sub["strand"].astype(str).to_numpy()
            plus = strands == "+"
        else:
            plus = np.ones(len(sub), dtype=bool)

        for s in starts[plus]:
            mer = seq[s: s + k]
            if len(mer) == k:
                counts[mer] += 1
        for e in ends[~plus]:
            mer = _revcomp(seq[e - k: e])
            if len(mer) == k:
                counts[mer] += 1
        if not has_strand:
            for e in ends:
                mer = _revcomp(seq[e - k: e])
                if len(mer) == k:
                    counts[mer] += 1

    total = sum(counts[m] for m in FOUR_MERS)
    if total == 0:
        return pd.Series(np.nan, index=[f"motif_{m}" for m in FOUR_MERS])
    freqs = np.array([counts[m] / total for m in FOUR_MERS])
    return pd.Series(freqs, index=[f"motif_{m}" for m in FOUR_MERS])


_COMP = str.maketrans("ACGTN", "TGCAN")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def motif_diversity_score(motif_freqs: pd.Series) -> float:
    """Normalised Shannon entropy of the 4-mer spectrum (MDS).

    Ranges 0-1; 1 means all 256 motifs equally likely. Reduced MDS has been
    reported in cancer plasma and reflects a narrower repertoire of cut sites.
    """
    p = motif_freqs.to_numpy(dtype=float)
    p = p[np.isfinite(p) & (p > 0)]
    if p.size == 0:
        return np.nan
    p = p / p.sum()
    return float(-(p * np.log(p)).sum() / np.log(len(FOUR_MERS)))
