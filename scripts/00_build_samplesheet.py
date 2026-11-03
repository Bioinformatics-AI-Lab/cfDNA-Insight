#!/usr/bin/env python3
"""Build samples.tsv from the extracted CRAG fragment files.

    python scripts/00_build_samplesheet.py --m-is control

THE LABEL MAPPING IS NOT DOCUMENTED ANYWHERE, SO THIS SCRIPT WILL NOT GUESS IT.

Every CRAG validation sample appears twice: `MAL0638.hg19.frag.bed.gz` and
`MAL0638m.hg19.frag.bed.gz`. The counts line up exactly with the paper's
cohort — 25 breast and 8 liver tumour plasma samples, each with a matched
control — so the pair is almost certainly case and matched control. But
neither the paper, the Zenodo record, nor the file names say which member of
the pair is which, and getting it backwards silently inverts every label in
the study.

Do not settle this by comparing fragment lengths between the two groups and
assuming the shorter one is cancer. That is the hypothesis under test; using
it to assign labels guarantees you "confirm" it. Settle it from Additional
file 1 Table S2 of the paper, or by emailing the corresponding author
(lyping1986@gmail.com), and record the answer in data/README.md.

Two sheets are written:

  config/samples.tsv
      study = crag_validation for every row. This is the truthful grouping:
      one cohort, one lab, one protocol. Leave-one-study-out will correctly
      report that cross-study generalisation is untested.

  config/samples_by_cancer_type.tsv
      study = crag_breast / crag_liver. Running the pipeline against this
      sheet performs leave-one-CANCER-TYPE-out instead: train on breast plus
      its controls, test on liver plus its controls. Weaker than a true
      cross-cohort test, since the protocol and lab are shared, but it does
      test whether the signal transfers across tumour biology. Report it
      under that name, never as cross-study validation.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

PATTERN = re.compile(r"^(MAL\d+)(m?)\.hg19\.frag\.bed\.gz$")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="data/raw/crag")
    ap.add_argument("--outdir", default="config")
    ap.add_argument(
        "--m-is", required=True, choices=["control", "cancer"],
        help="What the trailing 'm' means. REQUIRED — see the module "
             "docstring. Verify before trusting any result.")
    args = ap.parse_args()

    rows = []
    for dirpath, _, filenames in os.walk(args.raw):
        group = os.path.basename(dirpath)
        cancer_type = "liver" if group.startswith("liver") else "breast"
        for fn in sorted(filenames):
            m = PATTERN.match(fn)
            if not m:
                continue
            base, suffix = m.group(1), m.group(2)
            is_m = suffix == "m"
            if args.m_is == "control":
                label = "healthy" if is_m else "cancer"
            else:
                label = "cancer" if is_m else "healthy"
            rows.append({
                "sample_id": base + suffix,
                "path": os.path.join(dirpath, fn),
                "label": label,
                "pair_id": base,
                "cancer_type": cancer_type,
                "size_mb": round(os.path.getsize(os.path.join(dirpath, fn)) / 1e6, 1),
            })

    if not rows:
        print(f"[error] no *.hg19.frag.bed.gz under {args.raw}", file=sys.stderr)
        print("        run scripts/00_fetch_crag.sh first", file=sys.stderr)
        return 1

    rows.sort(key=lambda r: (r["cancer_type"], r["pair_id"], r["label"]))
    os.makedirs(args.outdir, exist_ok=True)

    cols = ["sample_id", "path", "label", "study", "cancer_type", "pair_id", "size_mb"]

    def write(path: str, study_of) -> None:
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                r = {**r, "study": study_of(r)}
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    main_sheet = os.path.join(args.outdir, "samples.tsv")
    alt_sheet = os.path.join(args.outdir, "samples_by_cancer_type.tsv")
    write(main_sheet, lambda r: "crag_validation")
    write(alt_sheet, lambda r: "crag_" + r["cancer_type"])

    n_cancer = sum(r["label"] == "cancer" for r in rows)
    unpaired = [p for p in {r["pair_id"] for r in rows}
                if sum(r["pair_id"] == p for r in rows) != 2]

    print(f"[ok] {len(rows)} samples: {n_cancer} cancer, "
          f"{len(rows) - n_cancer} healthy")
    for ct in sorted({r["cancer_type"] for r in rows}):
        sub = [r for r in rows if r["cancer_type"] == ct]
        print(f"     {ct}: {len(sub)} files, {len(sub) // 2} pairs")
    if unpaired:
        print(f"[warn] {len(unpaired)} ids without exactly 2 files: "
              f"{sorted(unpaired)[:5]}")
    print(f"[ok] wrote {main_sheet} and {alt_sheet}")
    print()
    print(f"[UNVERIFIED] labels assume trailing 'm' = {args.m_is}. "
          f"Confirm against Table S2 before reporting anything.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
