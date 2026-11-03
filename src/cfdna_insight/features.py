"""Per-sample feature assembly.

One sample in, one row of features out. Everything here is computed from a
single sample in isolation — no cross-sample normalisation, no reference
panel. That constraint is what makes it safe to run this step once, cache the
result, and do all leak-sensitive normalisation later inside CV folds.

The only cross-sample step in the whole pipeline is the healthy coverage
reference, and it lives in `cnv.HealthyReference`, fit inside each training
fold.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

from . import cnv, fragmentomics as fx
from .bins import make_bins
from .io import read_fragments


def compute_bin_gc(fasta_path: str, bin_size: int = 100_000,
                   genome: str = "hg38", cache: str | None = None
                   ) -> pd.DataFrame:
    """GC fraction per bin, computed once and cached.

    N-heavy bins get GC = NaN and are excluded from the GC fit; they would
    otherwise anchor the low-GC end of the correction curve with garbage.
    """
    if cache and os.path.exists(cache):
        return pd.read_parquet(cache)

    import pysam

    bins = make_bins(bin_size=bin_size, genome=genome)
    fasta = pysam.FastaFile(fasta_path)
    gc_vals = np.full(len(bins), np.nan)

    for chrom, sub in bins.groupby("chrom", observed=True):
        name = chrom if chrom in fasta.references else chrom.replace("chr", "")
        if name not in fasta.references:
            continue
        seq = fasta.fetch(name).upper()
        for i, (s, e) in zip(sub.index, zip(sub["start"], sub["end"])):
            window = seq[s:e]
            n_acgt = len(window) - window.count("N")
            if n_acgt < 0.5 * (e - s):
                continue
            gc_vals[i] = (window.count("G") + window.count("C")) / n_acgt

    bins["gc"] = gc_vals
    if cache:
        os.makedirs(os.path.dirname(cache) or ".", exist_ok=True)
        bins.to_parquet(cache, index=False)
    return bins


def extract_sample_features(
    frag_path: str,
    sample_id: str,
    gc_track: pd.DataFrame | None = None,
    fasta_path: str | None = None,
    bin_size: int = 100_000,
    window_factor: int = 50,
    genome: str = "hg38",
    min_mapq: int = 30,
    include_size_histogram: bool = True,
) -> dict[str, pd.Series]:
    """Return feature blocks for one sample, keyed by block name.

    Blocks are kept separate rather than concatenated so that the ablation
    study in `modeling` can train on any subset without string-matching on
    column prefixes.
    """
    frags = read_fragments(frag_path, min_mapq=min_mapq)
    if len(frags) == 0:
        raise ValueError(f"{sample_id}: no fragments survived filtering")

    binned = fx.bin_short_long(frags, bin_size=bin_size, genome=genome)
    gc = gc_track["gc"].to_numpy() if gc_track is not None else None

    blocks: dict[str, pd.Series] = {}
    blocks["size_summary"] = fx.size_summary(frags)
    if include_size_histogram:
        blocks["size_hist"] = fx.size_distribution(frags)
    blocks["delfi"] = fx.delfi_profile(binned, window_factor=window_factor, gc=gc)

    cov = cnv.coverage_windows(binned, window_factor=window_factor, gc=gc)
    blocks["_coverage"] = pd.Series(cov["cov"].to_numpy(),
                                    index=cov["bin_id"].to_numpy())
    blocks["_window_arms"] = pd.Series(cov["arm"].to_numpy(),
                                       index=cov["bin_id"].to_numpy())
    blocks["technical"] = cnv.technical_covariates(frags, binned)
    blocks["aneuploidy"] = cnv.aneuploidy_covariates(frags, binned)

    if fasta_path:
        import pysam
        fasta = pysam.FastaFile(fasta_path)
        motifs = fx.end_motifs(frags, fasta)
        blocks["motif"] = motifs
        blocks["size_summary"]["mds"] = fx.motif_diversity_score(motifs)

    for name, series in blocks.items():
        series.name = sample_id
    return blocks


def stack_blocks(per_sample: dict[str, dict[str, pd.Series]]
                 ) -> dict[str, pd.DataFrame]:
    """Turn {sample: {block: series}} into {block: samples x features frame}."""
    block_names = sorted({b for s in per_sample.values() for b in s})
    out = {}
    for block in block_names:
        rows = {sid: blocks[block] for sid, blocks in per_sample.items()
                if block in blocks}
        if not rows:
            continue
        out[block] = pd.DataFrame(rows).T
    return out


def save_features(blocks: dict[str, pd.DataFrame], outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    for name, df in blocks.items():
        df.to_parquet(os.path.join(outdir, f"{name.lstrip('_')}.parquet"))


def load_features(outdir: str) -> dict[str, pd.DataFrame]:
    blocks = {}
    for fn in sorted(os.listdir(outdir)):
        if fn.endswith(".parquet"):
            blocks[fn[: -len(".parquet")]] = pd.read_parquet(
                os.path.join(outdir, fn))
    return blocks


# Feature directories written before the technical/aneuploidy split carry
# chromosome entropy and coverage dispersion inside the technical block, where
# they corrupt the negative control. This moves them without re-extracting.
LEGACY_ANEUPLOIDY = {"tech_chr_entropy": "an_chr_entropy",
                     "tech_coverage_cv": "an_coverage_cv"}


def split_legacy_technical(blocks: dict[str, pd.DataFrame]
                           ) -> dict[str, pd.DataFrame]:
    tech = blocks.get("technical")
    if tech is None:
        return blocks
    present = [c for c in LEGACY_ANEUPLOIDY if c in tech.columns]
    if not present:
        return blocks
    moved = tech[present].rename(columns=LEGACY_ANEUPLOIDY)
    blocks = dict(blocks)
    blocks["technical"] = tech.drop(columns=present)
    blocks["aneuploidy"] = (
        pd.concat([blocks["aneuploidy"], moved], axis=1)
        if "aneuploidy" in blocks else moved)
    print(f"[info] moved {present} from technical -> aneuploidy "
          f"(they are copy-number readouts, not library metrics)")
    return blocks
