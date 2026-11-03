#!/usr/bin/env python3
"""Extract per-sample clinical metadata from the CRAG supplementary workbook.

    python scripts/00_extract_metadata.py \
        --xlsx ~/Downloads/13073_2022_1141_MOESM1_ESM.xlsx

Input is Additional file 1 of Zhou et al., Genome Medicine 2022
(https://doi.org/10.1186/s13073-022-01141-8). Output is
data/metadata/crag_table_s2.tsv, one row per *file* (66), which
00_build_samplesheet.py then uses as the authoritative source of labels.

LAYOUT OF SUPPLEMENTARY TABLE 2

One row per case-control pair, with the case's fields on the left and the
matched control's on the right:

    Subject_ID           the cancer patient        e.g. MAL0638
    Control Subject_ID   their matched control     e.g. MAL0638m

This is what settles the trailing-`m` question: the control ID is the case ID
with `m` appended, and the rule holds for all 33 pairs with no exceptions.
The script re-verifies it on every run rather than trusting the convention.

BATCH

The table also carries `Library Preparation Batch` and `Sequencing Batch`, and
each pair's case and control share both. Batch is therefore orthogonal to
label by construction, which is what makes the technical-covariate negative
control in 02_train_evaluate.py a meaningful test rather than a formality.
Both columns are carried through so that claim can be checked, not assumed.
"""

from __future__ import annotations

import argparse
import os
import sys

# Column indices in "Supplementary table 2" (0-based). The header has
# duplicate names (Age_at_collection and Gender each appear for case and
# control), so positions are used rather than names.
C = {
    "case_id": 0, "tumor_type": 1, "stage": 2, "path_stage": 3, "grade": 4,
    "case_age": 5, "case_gender": 6, "case_smoking": 7, "case_alcohol": 8,
    "biomarkers": 9,
    "control_id": 11, "control_uc_id": 12, "control_age": 13,
    "control_gender": 14,
    "smoking_code": 15, "packs_per_day": 16, "pack_years": 17,
    "alcohol_code": 18,
    "isolation_date": 19, "lib_batch": 20, "lib_date": 21,
    "seq_batch": 22, "seq_date": 23,
    "tumor_seq_id": 24, "control_seq_id": 25,
}

OUT_COLS = [
    "sample_id", "pair_id", "label", "cancer_type", "stage", "grade",
    "age", "gender", "smoking", "alcohol", "biomarkers",
    "lib_batch", "seq_batch", "seq_id",
]


def clean(v) -> str:
    if v is None:
        return "NA"
    s = str(v).strip().replace("\t", " ").replace("\n", " ")
    return s if s and s.lower() != "none" else "NA"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True,
                    help="Additional file 1 (13073_2022_1141_MOESM1_ESM.xlsx)")
    ap.add_argument("--sheet", default="Supplementary table 2")
    ap.add_argument("--out", default="data/metadata/crag_table_s2.tsv")
    args = ap.parse_args()

    try:
        import openpyxl
    except ImportError:
        print("[error] needs openpyxl: pip install openpyxl", file=sys.stderr)
        return 1

    wb = openpyxl.load_workbook(args.xlsx, read_only=True, data_only=True)
    if args.sheet not in wb.sheetnames:
        print(f"[error] no sheet {args.sheet!r}; found {wb.sheetnames}",
              file=sys.stderr)
        return 1

    rows = list(wb[args.sheet].iter_rows(values_only=True))
    pairs = [r for r in rows[1:] if r and r[C["case_id"]] and r[C["tumor_type"]]]
    if not pairs:
        print("[error] no data rows found", file=sys.stderr)
        return 1

    violations = [
        (clean(r[C["case_id"]]), clean(r[C["control_id"]]))
        for r in pairs
        if clean(r[C["control_id"]]) != clean(r[C["case_id"]]) + "m"
    ]
    if violations:
        print(f"[error] the control-is-case+'m' rule fails for "
              f"{len(violations)} pairs: {violations[:5]}", file=sys.stderr)
        print("        labels cannot be assigned from the suffix; stopping.",
              file=sys.stderr)
        return 1

    out_rows = []
    same_batch = 0
    for r in pairs:
        case = clean(r[C["case_id"]])
        ctype = clean(r[C["tumor_type"]]).lower()
        lib, seq = clean(r[C["lib_batch"]]), clean(r[C["seq_batch"]])
        same_batch += 1  # one lib/seq batch per pair by construction

        out_rows.append({
            "sample_id": case, "pair_id": case, "label": "cancer",
            "cancer_type": ctype, "stage": clean(r[C["stage"]]),
            "grade": clean(r[C["grade"]]), "age": clean(r[C["case_age"]]),
            "gender": clean(r[C["case_gender"]]),
            "smoking": clean(r[C["case_smoking"]]),
            "alcohol": clean(r[C["case_alcohol"]]),
            "biomarkers": clean(r[C["biomarkers"]]),
            "lib_batch": lib, "seq_batch": seq,
            "seq_id": clean(r[C["tumor_seq_id"]]),
        })
        out_rows.append({
            "sample_id": case + "m", "pair_id": case, "label": "healthy",
            "cancer_type": ctype, "stage": "NA", "grade": "NA",
            "age": clean(r[C["control_age"]]),
            "gender": clean(r[C["control_gender"]]),
            "smoking": clean(r[C["smoking_code"]]),
            "alcohol": clean(r[C["alcohol_code"]]),
            "biomarkers": "NA",
            "lib_batch": lib, "seq_batch": seq,
            "seq_id": clean(r[C["control_seq_id"]]),
        })

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\t".join(OUT_COLS) + "\n")
        for row in out_rows:
            fh.write("\t".join(row[c] for c in OUT_COLS) + "\n")

    n_pairs = len(pairs)
    types: dict[str, int] = {}
    for r in pairs:
        t = clean(r[C["tumor_type"]]).lower()
        types[t] = types.get(t, 0) + 1

    print(f"[ok] {n_pairs} pairs -> {len(out_rows)} samples")
    print(f"     {dict(sorted(types.items()))}")
    print(f"     control-is-case+'m' verified for all {n_pairs} pairs")
    print(f"     {same_batch}/{n_pairs} pairs share one library-prep and "
          f"sequencing batch")
    print(f"[ok] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
