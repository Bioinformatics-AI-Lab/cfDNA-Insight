"""Reading cfDNA fragment coordinate files.

FinaleDB / FinaleToolkit distribute cfDNA as *fragment* files rather than BAMs:
one row per sequenced fragment, with the genotype-bearing sequence stripped out.
That is what makes them openly downloadable. Columns are BED-like:

    chrom  start  end  [name]  mapq  strand

Some releases omit the `name` column, so column count is auto-detected.
Coordinates are 0-based half-open (BED convention); fragment length = end - start.
"""

from __future__ import annotations

import gzip
import io as _io
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd

FRAG_DTYPES = {"chrom": "category", "start": "int64", "end": "int64",
               "mapq": "int16", "strand": "category"}

# The 24 primary assemblies. Everything else (alts, decoys, chrM) is dropped:
# chrM in particular is wildly variable in cfDNA and is a classic source of
# spurious between-sample signal.
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def _open_maybe_gzip(path: str):
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "rt")


# Columns this pipeline actually reads. Everything else in the file is named
# but left unparsed, which also keeps the big string columns (CIGAR, name) out
# of memory on files with tens of millions of rows.
WANTED = ("chrom", "start", "end", "mapq", "strand")

# Header spellings seen in the wild, mapped to our names.
_ALIASES = {
    "chrom": "chrom", "chr": "chrom", "chromosome": "chrom", "#chrom": "chrom",
    "start": "start", "chromstart": "start",
    "end": "end", "chromend": "end", "stop": "end",
    "name": "name", "id": "name",
    "mapq": "mapq", "score": "mapq", "quality": "mapq",
    "strand": "strand",
}


def _sniff_columns(path: str) -> list[str]:
    """Name every column in the file, so pandas is never short a name.

    Two layouts occur. FinaleDB ships bare BED-like rows with no header. The
    CRAG fragment files carry a `#`-prefixed header and two extra CIGAR
    columns:

        #chrom  start  end  name  mapq  strand  cigar1  cigar2

    A fixed 6-name guess silently breaks on the second: pandas rejects a
    `names` list shorter than the row when `usecols` is in play. So the header
    is used when there is one, and otherwise names are assigned positionally
    with any surplus columns named `extra_N` rather than dropped.
    """
    header: list[str] | None = None
    fields: list[str] | None = None

    with _open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                if header is None:
                    header = line.lstrip("#").split("\t")
                continue
            if line.startswith(("track", "browser")):
                continue
            fields = line.split("\t")
            break

    if fields is None:
        raise ValueError(f"{path} contains no data rows")
    n = len(fields)
    if n < 3:
        raise ValueError(f"{path}: expected >=3 columns, found {n}")

    if header is not None and len(header) == n:
        names, seen = [], set()
        for i, h in enumerate(header):
            nm = _ALIASES.get(h.strip().lower().lstrip("#"), None)
            if nm is None or nm in seen:
                nm = f"extra_{i}"
            seen.add(nm)
            names.append(nm)
        if {"chrom", "start", "end"} <= set(names):
            return names
        # header present but unrecognisable — fall through to positional

    if n == 3:
        return ["chrom", "start", "end"]
    if n == 4:
        # 4th column is mapq if it parses as an integer, else a name field
        try:
            int(fields[3])
            return ["chrom", "start", "end", "mapq"]
        except ValueError:
            return ["chrom", "start", "end", "name"]
    if n == 5:
        try:
            int(fields[3])
            return ["chrom", "start", "end", "mapq", "strand"]
        except ValueError:
            return ["chrom", "start", "end", "name", "mapq"]

    base = ["chrom", "start", "end", "name", "mapq", "strand"]
    return base + [f"extra_{i}" for i in range(6, n)]


def read_fragments(
    path: str,
    min_mapq: int = 30,
    min_length: int = 50,
    max_length: int = 1000,
    chroms: list[str] | None = None,
    add_chr_prefix: bool = True,
) -> pd.DataFrame:
    """Load a fragment file into a DataFrame with a `length` column.

    Parameters mirror the filters used in the DELFI and FinaleDB papers.
    `max_length=1000` discards chimeric/discordant pairs; `min_length=50`
    discards sub-nucleosomal noise that is mostly adapter artefact.
    """
    names = _sniff_columns(path)
    usecols = [c for c in names if c in WANTED]

    df = pd.read_csv(
        path, sep="\t", header=None, names=names, usecols=usecols,
        comment="#", engine="c",
        dtype={"chrom": str, "start": "int64", "end": "int64",
               "mapq": "int16", "strand": "category"},
    )

    if add_chr_prefix:
        # FinaleDB ships both hg19 (no prefix) and hg38 (prefixed) builds
        needs_prefix = ~df["chrom"].str.startswith("chr")
        df.loc[needs_prefix, "chrom"] = "chr" + df.loc[needs_prefix, "chrom"]

    keep = chroms if chroms is not None else AUTOSOMES
    df = df[df["chrom"].isin(keep)]

    if "mapq" in df.columns and min_mapq is not None:
        df = df[df["mapq"] >= min_mapq]

    df["length"] = df["end"] - df["start"]
    df = df[(df["length"] >= min_length) & (df["length"] <= max_length)]

    df["chrom"] = pd.Categorical(df["chrom"], categories=keep, ordered=True)
    return df.reset_index(drop=True)


@dataclass
class SampleSheet:
    """Sample manifest. `study` is mandatory and load-bearing.

    Every downstream split is grouped by `study` so that a model is never
    evaluated on the same cohort/batch it was trained on without that being
    an explicit, labelled choice.
    """

    table: pd.DataFrame

    REQUIRED = ("sample_id", "path", "label", "study")

    @classmethod
    def load(cls, path: str) -> "SampleSheet":
        sep = "," if str(path).endswith(".csv") else "\t"
        df = pd.read_csv(path, sep=sep)
        missing = [c for c in cls.REQUIRED if c not in df.columns]
        if missing:
            raise ValueError(f"sample sheet missing required columns: {missing}")
        if df["sample_id"].duplicated().any():
            dupes = df.loc[df["sample_id"].duplicated(), "sample_id"].tolist()
            raise ValueError(f"duplicate sample_id: {dupes}")
        bad = set(df["label"].unique()) - {"cancer", "healthy"}
        if bad:
            raise ValueError(f"label must be cancer/healthy, found: {bad}")
        return cls(df)

    def resolve(self, root: str = ".") -> "SampleSheet":
        df = self.table.copy()
        df["path"] = df["path"].apply(
            lambda p: p if os.path.isabs(p) else os.path.join(root, p)
        )
        return SampleSheet(df)

    def existing(self) -> "SampleSheet":
        mask = self.table["path"].apply(os.path.exists)
        if not mask.all():
            missing = self.table.loc[~mask, "sample_id"].tolist()
            print(f"[warn] {len(missing)} sample files not found, skipping: "
                  f"{missing[:5]}{'...' if len(missing) > 5 else ''}")
        return SampleSheet(self.table[mask].reset_index(drop=True))

    @property
    def y(self) -> np.ndarray:
        return (self.table["label"] == "cancer").astype(int).to_numpy()

    @property
    def groups(self) -> np.ndarray:
        return self.table["study"].to_numpy()

    def __len__(self) -> int:
        return len(self.table)
