"""Synthetic fragment files.

Purpose: let the test suite and a first end-to-end smoke run work before any
real data is downloaded. These are *not* a substitute for real cfDNA and no
result computed on them belongs in the write-up. They exist so that a broken
function fails in CI rather than three hours into a real run.

The generator models the two things the pipeline is supposed to detect:
a shift toward shorter fragments, and arm-level coverage imbalance.
"""

from __future__ import annotations

import gzip
import os

import numpy as np

from .bins import HG38_SIZES

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def _sample_lengths(n: int, short_shift: float,
                    rng: np.random.Generator) -> np.ndarray:
    """Mono- plus di-nucleosomal mixture with 10 bp periodicity."""
    n_di = rng.binomial(n, 0.06)
    n_mono = n - n_di

    mono = rng.normal(167 - 22 * short_shift, 28 + 6 * short_shift, n_mono)
    # superimpose the helical-periodicity ripple
    mono += 4.0 * np.sin(2 * np.pi * mono / 10.0)
    di = rng.normal(334, 45, n_di)

    lengths = np.concatenate([mono, di])
    return np.clip(lengths, 50, 1000).astype(int)


def simulate_sample(
    path: str,
    n_fragments: int = 200_000,
    is_cancer: bool = False,
    effect_size: float = 1.0,
    seed: int = 0,
    chroms: list[str] | None = None,
    genome_scale: float = 0.02,
) -> str:
    """Write a gzipped FinaleDB-style fragment file.

    `genome_scale` shrinks the coordinate space so test files stay small while
    still exercising the binning code; set to 1.0 for realistic coordinates.
    """
    rng = np.random.default_rng(seed)
    chroms = chroms or AUTOSOMES
    sizes = {c: int(HG38_SIZES[c] * genome_scale) for c in chroms}

    weights = np.array([sizes[c] for c in chroms], dtype=float)
    if is_cancer:
        # arm-level gains/losses on a few chromosomes
        for i in (0, 7, 12):
            if i < len(weights):
                weights[i] *= 1.0 + 0.25 * effect_size
        for i in (3, 16):
            if i < len(weights):
                weights[i] *= 1.0 - 0.20 * effect_size
    weights /= weights.sum()

    counts = rng.multinomial(n_fragments, weights)
    short_shift = effect_size if is_cancer else 0.0

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with gzip.open(path, "wt") as out:
        for chrom, k in zip(chroms, counts):
            if k == 0:
                continue
            lengths = _sample_lengths(k, short_shift, rng)
            starts = rng.integers(0, max(sizes[chrom] - 1000, 1), k)
            ends = starts + lengths
            mapqs = rng.integers(30, 61, k)
            strands = rng.choice(["+", "-"], k)
            order = np.argsort(starts)
            for i in order:
                out.write(f"{chrom}\t{starts[i]}\t{ends[i]}\t"
                          f"{mapqs[i]}\t{strands[i]}\n")
    return path


def simulate_cohort(outdir: str, n_per_group: int = 12,
                    n_fragments: int = 150_000, effect_size: float = 1.0,
                    study: str = "SIM01", seed: int = 0) -> str:
    """Write a cohort plus a matching sample sheet. Returns the sheet path."""
    os.makedirs(outdir, exist_ok=True)
    rows = ["sample_id\tpath\tlabel\tstudy"]
    for i in range(n_per_group):
        for label, is_cancer in (("healthy", False), ("cancer", True)):
            sid = f"{study}_{label}_{i:02d}"
            fn = os.path.join(outdir, f"{sid}.frag.gz")
            simulate_sample(fn, n_fragments=n_fragments, is_cancer=is_cancer,
                            effect_size=effect_size,
                            seed=seed + i * 2 + int(is_cancer))
            rows.append(f"{sid}\t{fn}\t{label}\t{study}")
    sheet = os.path.join(outdir, "samples.tsv")
    with open(sheet, "w") as fh:
        fh.write("\n".join(rows) + "\n")
    return sheet
