"""Genomic binning.

The DELFI design uses 100 kb sub-bins aggregated into 5 Mb windows. The 100 kb
level is where GC correction happens (GC content varies at that scale); the
5 Mb level is where features are reported, because short/long ratios in a
single 100 kb bin at 1-2x coverage are far too noisy to be a feature.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .io import AUTOSOMES

# hg38 primary autosome lengths (GRCh38.p14, no alts).
HG38_SIZES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468,
}

HG19_SIZES = {
    "chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276,
    "chr5": 180915260, "chr6": 171115067, "chr7": 159138663, "chr8": 146364022,
    "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753,
    "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
    "chr21": 48129895, "chr22": 51304566,
}

# Centromere / acrocentric-arm spans, used to assign p/q arms and to drop
# bins that straddle the centromere. Approximate, hg38.
HG38_CENTROMERE = {
    "chr1": (121700000, 125100000), "chr2": (91800000, 96000000),
    "chr3": (87800000, 94000000), "chr4": (48200000, 51800000),
    "chr5": (46100000, 50100000), "chr6": (58500000, 62600000),
    "chr7": (58100000, 62100000), "chr8": (43200000, 47200000),
    "chr9": (42200000, 45500000), "chr10": (38000000, 41600000),
    "chr11": (51000000, 55800000), "chr12": (33200000, 37800000),
    "chr13": (16500000, 18900000), "chr14": (16100000, 18200000),
    "chr15": (17500000, 20500000), "chr16": (35300000, 38400000),
    "chr17": (22700000, 27400000), "chr18": (15400000, 21500000),
    "chr19": (24200000, 28100000), "chr20": (25700000, 30400000),
    "chr21": (10900000, 13000000), "chr22": (13700000, 17400000),
}

# hg19 centromere spans (UCSC gap track). These differ from hg38 by several
# megabases on some chromosomes, which is enough to misassign whole 5 Mb
# windows to the wrong arm — so the build must be selected, not assumed.
HG19_CENTROMERE = {
    "chr1": (121535434, 124535434), "chr2": (92326171, 95326171),
    "chr3": (90504854, 93504854), "chr4": (49660117, 52660117),
    "chr5": (46405641, 49405641), "chr6": (58830166, 61830166),
    "chr7": (58054331, 61054331), "chr8": (43838887, 46838887),
    "chr9": (47367679, 50367679), "chr10": (39254935, 42254935),
    "chr11": (51644205, 54644205), "chr12": (34856694, 37856694),
    "chr13": (16000000, 19000000), "chr14": (16000000, 19000000),
    "chr15": (17000000, 20000000), "chr16": (35335801, 38335801),
    "chr17": (22263006, 25263006), "chr18": (15460898, 18460898),
    "chr19": (24681782, 27681782), "chr20": (26369569, 29369569),
    "chr21": (11288129, 14288129), "chr22": (13000000, 16000000),
}

SUPPORTED_GENOMES = ("hg38", "hg19")


def _check_genome(genome: str) -> str:
    if genome not in SUPPORTED_GENOMES:
        raise ValueError(
            f"unsupported genome {genome!r}; expected one of {SUPPORTED_GENOMES}")
    return genome


def get_sizes(genome: str = "hg38") -> dict[str, int]:
    return dict(HG38_SIZES if _check_genome(genome) == "hg38" else HG19_SIZES)


def get_centromeres(genome: str = "hg38") -> dict[str, tuple[int, int]]:
    return dict(HG38_CENTROMERE if _check_genome(genome) == "hg38"
                else HG19_CENTROMERE)


def make_bins(bin_size: int = 100_000, genome: str = "hg38",
              chroms: list[str] | None = None) -> pd.DataFrame:
    """Fixed-width bins tiling the autosomes, with arm annotation."""
    sizes = get_sizes(genome)
    centromeres = get_centromeres(genome)
    chroms = chroms or AUTOSOMES
    rows = []
    for chrom in chroms:
        n = sizes[chrom] // bin_size
        starts = np.arange(n, dtype=np.int64) * bin_size
        cen_start, cen_end = centromeres[chrom]
        mids = starts + bin_size // 2
        arm = np.where(mids < cen_start, "p", "q")
        in_cen = (starts + bin_size > cen_start) & (starts < cen_end)
        rows.append(pd.DataFrame({
            "chrom": chrom,
            "start": starts,
            "end": starts + bin_size,
            "arm": [f"{chrom}{a}" for a in arm],
            "in_centromere": in_cen,
        }))
    bins = pd.concat(rows, ignore_index=True)
    bins["bin_id"] = bins["chrom"] + ":" + bins["start"].astype(str)
    return bins


def assign_bins(frags: pd.DataFrame, bin_size: int) -> np.ndarray:
    """Bin index by fragment midpoint.

    Midpoint rather than start: using the start biases assignment toward the
    upstream bin by ~85 bp on average, which is negligible at 5 Mb but not at
    the sub-kb scales used for nucleosome work.
    """
    mid = (frags["start"].to_numpy() + frags["end"].to_numpy()) // 2
    return mid // bin_size


def aggregate_to_windows(binned: pd.DataFrame, factor: int = 50,
                         value_cols: list[str] | None = None) -> pd.DataFrame:
    """Roll 100 kb bins up into 5 Mb windows (factor=50) by summing counts.

    Aggregation is within-chromosome; a window is only emitted if it is
    complete, so the ragged tail of each chromosome is dropped rather than
    silently producing a short, low-count window.
    """
    value_cols = value_cols or [c for c in binned.columns
                                if c.startswith(("n_", "count"))]
    out = []
    for chrom, sub in binned.groupby("chrom", observed=True, sort=False):
        sub = sub.sort_values("start").reset_index(drop=True)
        n_full = len(sub) // factor
        if n_full == 0:
            continue
        sub = sub.iloc[: n_full * factor]
        grp = np.repeat(np.arange(n_full), factor)
        agg = sub.groupby(grp)[value_cols].sum()
        agg["chrom"] = chrom
        agg["start"] = sub.groupby(grp)["start"].first().to_numpy()
        agg["end"] = sub.groupby(grp)["end"].last().to_numpy()
        # a window inherits the arm of its majority sub-bin
        agg["arm"] = sub.groupby(grp)["arm"].agg(
            lambda s: s.value_counts().index[0]).to_numpy()
        agg["frac_centromeric"] = sub.groupby(grp)["in_centromere"].mean().to_numpy()
        out.append(agg)
    res = pd.concat(out, ignore_index=True)
    res["bin_id"] = res["chrom"] + ":" + res["start"].astype(str)
    return res


def load_blacklist(path: str | None) -> pd.DataFrame | None:
    """ENCODE blacklist BED. Optional but strongly recommended."""
    if path is None:
        return None
    bl = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end"], usecols=[0, 1, 2],
                     comment="#", dtype={"chrom": str})
    if not bl["chrom"].str.startswith("chr").all():
        bl["chrom"] = "chr" + bl["chrom"].str.replace("^chr", "", regex=True)
    return bl


def blacklist_overlap_fraction(bins: pd.DataFrame,
                               blacklist: pd.DataFrame | None) -> np.ndarray:
    """Fraction of each bin covered by blacklist intervals."""
    if blacklist is None:
        return np.zeros(len(bins))
    frac = np.zeros(len(bins))
    for chrom, sub in bins.groupby("chrom", observed=True):
        bl = blacklist[blacklist["chrom"] == chrom]
        if bl.empty:
            continue
        s = sub["start"].to_numpy()[:, None]
        e = sub["end"].to_numpy()[:, None]
        bs = bl["start"].to_numpy()[None, :]
        be = bl["end"].to_numpy()[None, :]
        ov = np.clip(np.minimum(e, be) - np.maximum(s, bs), 0, None).sum(axis=1)
        frac[sub.index.to_numpy()] = ov / (sub["end"] - sub["start"]).to_numpy()
    return np.clip(frac, 0, 1)
