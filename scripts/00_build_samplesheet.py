#!/usr/bin/env python3
"""Build samples.tsv from the extracted CRAG fragment files.

    python scripts/00_extract_metadata.py --xlsx <Additional file 1>.xlsx
    python scripts/00_build_samplesheet.py

LABELS

Labels come from data/metadata/crag_table_s2.tsv, extracted from Supplementary
Table 2 of Zhou et al., Genome Medicine 2022. That table is laid out one row
per case-control pair, with `Subject_ID` the cancer patient and
`Control Subject_ID` the matched control — and the control ID is the case ID
with `m` appended, for all 33 pairs.

So `MAL0638` is the case and `MAL0638m` is the matched control. This was an
open question earlier in the project and is now settled from the source.

If the metadata file is absent, `--m-is` applies the same rule by convention
and every row is flagged unverified. Prefer the metadata file: it also carries
stage, grade, age, gender, smoking, alcohol and batch.

STUDY GROUPING

Two sheets are written, because the honest grouping and the useful one differ:

  config/samples.tsv
      study = crag_validation for every row. One cohort, one lab, one
      protocol. Leave-one-study-out will correctly report cross-study
      generalisation as untested.

  config/samples_by_cancer_type.tsv
      study = crag_breast / crag_liver. Running against this sheet performs
      leave-one-CANCER-TYPE-out: train on breast plus its controls, test on
      liver plus its controls. Weaker than a true cross-cohort test — protocol
      and lab are shared — but it does test whether the signal transfers
      across tumour biology. Report it under that name, never as cross-study.

A third grouping is available from the metadata: `seq_batch` (Cabat1-5). Note
that grouping by batch is NOT a confounding test here, because each pair's
case and control were processed in the same batch, so batch is orthogonal to
label by construction. That is a property worth stating in the write-up, and
it is why the technical-covariate negative control is a real test.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys

PATTERN = re.compile(r"^(MAL\d+)(m?)\.hg19\.frag\.bed\.gz$")
META_DEFAULT = "data/metadata/crag_table_s2.tsv"

EXTRA = ["stage", "grade", "age", "gender", "smoking", "alcohol",
         "lib_batch", "seq_batch", "seq_id"]


def load_metadata(path: str) -> dict[str, dict[str, str]]:
    with open(path) as fh:
        return {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="data/raw/crag")
    ap.add_argument("--outdir", default="config")
    ap.add_argument("--metadata", default=META_DEFAULT)
    ap.add_argument("--m-is", choices=["control", "cancer"], default=None,
                    help="Fallback only, used when --metadata is missing.")
    args = ap.parse_args()

    meta = {}
    if os.path.exists(args.metadata):
        meta = load_metadata(args.metadata)
        print(f"[info] labels from {args.metadata} ({len(meta)} samples)")
    elif args.m_is:
        print(f"[warn] {args.metadata} not found; falling back to "
              f"--m-is {args.m_is} (UNVERIFIED)")
    else:
        print(f"[error] {args.metadata} not found and no --m-is given.\n"
              f"        run scripts/00_extract_metadata.py first, or pass "
              f"--m-is control to apply the convention unverified.",
              file=sys.stderr)
        return 1

    rows, unknown = [], []
    for dirpath, _, filenames in os.walk(args.raw):
        for fn in sorted(filenames):
            m = PATTERN.match(fn)
            if not m:
                continue
            sid = m.group(1) + m.group(2)
            is_m = m.group(2) == "m"
            path = os.path.join(dirpath, fn)
            row = {
                "sample_id": sid, "path": path, "pair_id": m.group(1),
                "size_mb": round(os.path.getsize(path) / 1e6, 1),
            }
            if meta:
                if sid not in meta:
                    unknown.append(sid)
                    continue
                md = meta[sid]
                row["label"] = md["label"]
                row["cancer_type"] = md["cancer_type"]
                row.update({k: md.get(k, "NA") for k in EXTRA})
            else:
                row["label"] = ("healthy" if is_m else "cancer") \
                    if args.m_is == "control" else \
                    ("cancer" if is_m else "healthy")
                row["cancer_type"] = ("liver"
                                      if os.path.basename(dirpath).startswith("liver")
                                      else "breast")
                row.update({k: "NA" for k in EXTRA})
            rows.append(row)

    if not rows:
        print(f"[error] no *.hg19.frag.bed.gz under {args.raw}; run "
              f"scripts/00_fetch_crag.sh first", file=sys.stderr)
        return 1
    if unknown:
        print(f"[warn] {len(unknown)} files absent from the metadata, "
              f"skipped: {sorted(unknown)[:5]}")

    rows.sort(key=lambda r: (r["cancer_type"], r["pair_id"], r["label"]))
    os.makedirs(args.outdir, exist_ok=True)
    cols = (["sample_id", "path", "label", "study", "cancer_type", "pair_id"]
            + EXTRA + ["size_mb"])

    def write(path: str, study_of) -> None:
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                r = {**r, "study": study_of(r)}
                fh.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")

    main_sheet = os.path.join(args.outdir, "samples.tsv")
    alt_sheet = os.path.join(args.outdir, "samples_by_cancer_type.tsv")
    write(main_sheet, lambda r: "crag_validation")
    write(alt_sheet, lambda r: "crag_" + r["cancer_type"])

    n_cancer = sum(r["label"] == "cancer" for r in rows)
    print(f"[ok] {len(rows)} samples: {n_cancer} cancer, "
          f"{len(rows) - n_cancer} healthy")
    for ct in sorted({r["cancer_type"] for r in rows}):
        sub = [r for r in rows if r["cancer_type"] == ct]
        print(f"     {ct}: {len(sub)} files "
              f"({sum(s['label'] == 'cancer' for s in sub)} cancer)")

    broken = [p for p in {r["pair_id"] for r in rows}
              if sorted(r["label"] for r in rows if r["pair_id"] == p)
              != ["cancer", "healthy"]]
    if broken:
        print(f"[warn] {len(broken)} pair_ids lack exactly one cancer and one "
              f"healthy file: {sorted(broken)[:5]}")
        print("       (expected while only part of the cohort is downloaded)")

    if meta:
        batches = {r["seq_batch"] for r in rows if r["seq_batch"] != "NA"}
        mixed = [p for p in {r["pair_id"] for r in rows}
                 if len({r["seq_batch"] for r in rows if r["pair_id"] == p}) > 1]
        print(f"[ok] {len(batches)} sequencing batches; "
              f"{len(mixed)} pairs split across batches "
              f"({'none — batch is orthogonal to label' if not mixed else 'CHECK THIS'})")
    else:
        print(f"\n[UNVERIFIED] labels assume trailing 'm' = {args.m_is}.")

    print(f"[ok] wrote {main_sheet} and {alt_sheet}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
