#!/usr/bin/env python3
"""Decide whether MAL#### and MAL####m are two people or one library twice.

    python scripts/00_check_pairing.py --raw data/raw/crag/liver

THE QUESTION

Every CRAG validation sample appears twice: `MAL0646.hg19.frag.bed.gz` and
`MAL0646m.hg19.frag.bed.gz`. Two readings fit the file counts equally well:

  (a) case + matched control  -> 33 cancer + 33 healthy. Usable.
  (b) one sample, two versions (merged runs, replicate libraries)
      -> 33 cancer and no controls at all. Not usable for this project.

THE TEST

Compare the exact genomic coordinates of fragments, not their lengths.

Two sequencing runs of the same library sample the same physical cfDNA
molecules, so a large fraction of (start, end) pairs recur exactly. Two
different people share a given (start, end) pair only when the same molecule
boundaries arose independently in both.

Fragment length is never used, so this says nothing about the cancer-vs-healthy
hypothesis and cannot contaminate it.

THE BASELINE MATTERS

cfDNA endpoints are not uniformly distributed — they cluster at nucleosome
boundaries, so even unrelated people share coordinates well above what a
uniform model would predict. The absolute overlap is therefore uninterpretable
on its own. Each pair is scored against an unrelated sample from the same
cohort, and only the *ratio* is read.

WHAT IT DOES NOT ANSWER

If the result is (a), this tells you the two files are different people. It
does NOT tell you which one is the cancer case. That still has to come from
Table S2 of the paper or from the authors.

Comparison is restricted to one chromosome (chr21 by default) to keep memory
bounded; that is a few hundred thousand fragments, far more than enough.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys

PATTERN = re.compile(r"^(MAL\d+)(m?)\.hg19\.frag\.bed\.gz$")


def load_coords(path: str, chrom: str, min_mapq: int = 30) -> set[int]:
    """Exact (start, end) pairs on one chromosome, packed into ints.

    Streams the gzip rather than loading a DataFrame: these files hold tens of
    millions of rows and only one chromosome is needed.
    """
    want = {chrom, chrom.replace("chr", "")}
    out: set[int] = set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 5)
            if f[0] not in want:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except (ValueError, IndexError):
                continue
            length = end - start
            if not (50 <= length <= 1000):
                continue
            if min_mapq and len(f) > 4:
                try:
                    if int(f[3]) < min_mapq:
                        continue
                except ValueError:
                    pass
            out.add((start << 11) | length
                    if length < 2048 else (start << 11))
    return out


def containment(a: set[int], b: set[int]) -> float:
    """|A and B| / min(|A|, |B|).

    Preferred over Jaccard here because sequencing depth differs several-fold
    between files in this cohort, which deflates Jaccard for reasons that have
    nothing to do with sample identity.
    """
    if not a or not b:
        return float("nan")
    return len(a & b) / min(len(a), len(b))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="data/raw/crag",
                    help="directory holding the extracted *.frag.bed.gz files")
    ap.add_argument("--chrom", default="chr21")
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--max-pairs", type=int, default=8,
                    help="stop after this many pairs (each pair reads 2 files)")
    args = ap.parse_args()

    pairs: dict[str, dict[str, str]] = {}
    for dirpath, _, filenames in os.walk(args.raw):
        for fn in sorted(filenames):
            m = PATTERN.match(fn)
            if m:
                pairs.setdefault(m.group(1), {})[m.group(2) or "base"] = \
                    os.path.join(dirpath, fn)

    complete = sorted(k for k, v in pairs.items() if {"base", "m"} <= set(v))
    if not complete:
        print(f"[error] no complete MAL####/MAL####m pairs under {args.raw}",
              file=sys.stderr)
        return 1
    complete = complete[: args.max_pairs]
    print(f"[info] {len(complete)} pairs on {args.chrom}, MAPQ>={args.min_mapq}\n")

    cache: dict[str, set[int]] = {}

    def coords(path: str) -> set[int]:
        if path not in cache:
            cache[path] = load_coords(path, args.chrom, args.min_mapq)
        return cache[path]

    print(f"{'pair':<10}{'n_base':>10}{'n_m':>10}{'self':>9}"
          f"{'unrelated':>11}{'ratio':>8}")
    print("-" * 58)

    ratios = []
    for i, pid in enumerate(complete):
        a = coords(pairs[pid]["base"])
        b = coords(pairs[pid]["m"])
        self_c = containment(a, b)

        other = complete[(i + 1) % len(complete)]
        if other == pid:
            print(f"{pid:<10} only one pair — no unrelated baseline available")
            continue
        unrel_c = containment(a, coords(pairs[other]["base"]))

        ratio = self_c / unrel_c if unrel_c > 0 else float("inf")
        ratios.append(ratio)
        print(f"{pid:<10}{len(a):>10,}{len(b):>10,}{self_c:>9.3f}"
              f"{unrel_c:>11.3f}{ratio:>8.1f}")

    if not ratios:
        return 1

    med = sorted(ratios)[len(ratios) // 2]
    print("\n" + "=" * 58)
    print(f"median self/unrelated containment ratio: {med:.1f}")
    print()
    if med > 3:
        print("VERDICT: reading (b) — X and Xm look like the SAME library.")
        print("  Self-overlap far exceeds the unrelated baseline, which means")
        print("  the two files sample the same physical molecules. If so this")
        print("  record contains 33 cancer samples and NO healthy controls,")
        print("  and cannot support a case-control design.")
    elif med < 1.5:
        print("VERDICT: reading (a) — X and Xm look like DIFFERENT people.")
        print("  Self-overlap is indistinguishable from the unrelated")
        print("  baseline, consistent with case + matched control.")
        print("  This does NOT say which one is the case — that still needs")
        print("  Table S2 or the authors.")
    else:
        print("VERDICT: inconclusive. Re-run on another chromosome")
        print("  (--chrom chr20) and with more pairs before concluding.")
    print("=" * 58)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
