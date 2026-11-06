#!/usr/bin/env python3
"""In silico dilution series and limit-of-detection estimate.

    python scripts/03_dilution_lod.py \
        --features results/features --cancer MAL0646 --healthy MAL0646m \
        --fasta data/ref/hg19.fa --genome hg19

WHAT THIS ANSWERS

Case-control AUC says a classifier can separate two groups. It does not say
how much tumour signal has to be present before that works, which is the
question that decides whether an assay is useful for early detection. Mixing
cancer-patient plasma into healthy plasma at known ratios and rescoring gives
a direct answer without new sequencing.

WHICH FEATURES

By default the series uses `aneuploidy` (chromosome entropy and coverage
dispersion), because on this cohort those are the only features that separate
cases from controls at all, and because neither needs a healthy reference
panel. The second point is a hard constraint, not a convenience: arm-level
CNV z-scores are defined against a panel, and a panel built from the same
donors that produced the mixture would leak.

`--blocks` can select others. `size_summary,size_hist,delfi` reproduces the
fragmentomic version, which on this cohort produces a flat line, since those
features were at chance to begin with.

DONOR EXCLUSION

The model is trained on every donor except the two being mixed. A model that
has seen the cancer donor will recognise it at any dilution for the wrong
reason.

WHAT THE X-AXIS IS NOT

The mixing fraction is a fraction of cancer-patient *plasma*, not of tumour
DNA. The donor's plasma is itself mostly haematopoietic cfDNA, so a 10%
mixture of a sample carrying 5% ctDNA is roughly 0.5% ctDNA. Every LOD here
is an upper bound on the true ctDNA LOD, by a factor equal to the donor's own
tumour fraction, which this cohort does not report.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from cfdna_insight import modeling, plots  # noqa: E402
from cfdna_insight.dilution import (  # noqa: E402
    check_direction, dilution_series, dose_response, estimate_lod,
    threshold_at_specificity,
)
from cfdna_insight.features import (  # noqa: E402
    compute_bin_gc, load_features, split_legacy_technical,
)

REFERENCE_FREE = ("aneuploidy", "size_summary", "size_hist", "delfi")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default="results/features")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--cancer", required=True, help="sample_id of cancer donor")
    ap.add_argument("--healthy", required=True,
                    help="sample_id of healthy donor")
    ap.add_argument("--blocks", default="aneuploidy",
                    help="comma-separated feature blocks. Only blocks "
                         "recomputable without a reference panel are allowed: "
                         + ", ".join(REFERENCE_FREE))
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--genome", default="hg19", choices=["hg19", "hg38"])
    ap.add_argument("--gc-cache", default="results/cache/bin_gc.parquet")
    ap.add_argument("--total-fragments", type=int, default=5_000_000)
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--fractions", default=None,
                    help="comma-separated, e.g. 0,0.01,0.05,0.1,0.25,0.5,1")
    args = ap.parse_args()

    use = [b.strip() for b in args.blocks.split(",") if b.strip()]
    bad = [b for b in use if b not in REFERENCE_FREE]
    if bad:
        print(f"[error] {bad} need a healthy reference panel and cannot be "
              f"recomputed on a mixture without leaking. Allowed: "
              f"{list(REFERENCE_FREE)}", file=sys.stderr)
        return 1

    blocks = split_legacy_technical(load_features(args.features))
    samples = pd.read_csv(os.path.join(args.features, "samples_used.tsv"),
                          sep="\t")
    paths = dict(zip(samples["sample_id"], samples["path"]))
    for sid in (args.cancer, args.healthy):
        if sid not in paths:
            print(f"[error] {sid} not in the sample sheet", file=sys.stderr)
            return 1
    if not os.path.exists(paths[args.cancer]):
        print(f"[error] fragment file missing: {paths[args.cancer]}",
              file=sys.stderr)
        return 1

    # hold both donors out of training
    train = samples[~samples["sample_id"].isin([args.cancer, args.healthy])]
    order = [s for s in train["sample_id"] if s in blocks["size_summary"].index]
    y = (train.set_index("sample_id").loc[order, "label"] == "cancer"
         ).astype(int).to_numpy()

    X = modeling.assemble_matrix(blocks, use, order)
    print(f"[info] locked model: {len(order)} donors, {X.shape[1]} features "
          f"from {use}")
    print(f"[info] mixing {args.cancer} (cancer) into {args.healthy} (healthy)")

    model = modeling.fit_final(X, y)
    baseline = modeling.nested_cv(X, y, n_repeats=5)
    thresh = threshold_at_specificity(baseline.oof_scores[y == 0], 0.95)
    print(f"[info] held-out AUC on training donors {baseline.auc_mean:.3f}; "
          f"threshold at 95% specificity = {thresh:.3f}")
    if baseline.auc_mean < 0.65:
        print(f"[warn] the locked model is near chance on its own training "
              f"donors. A dilution curve from it describes noise, not a "
              f"limit of detection.")

    gc_track = None
    if args.fasta:
        gc_track = compute_bin_gc(args.fasta, genome=args.genome,
                                  cache=args.gc_cache)

    fractions = ([float(f) for f in args.fractions.split(",")]
                 if args.fractions else
                 (0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 1.0))

    series = dilution_series(
        paths[args.cancer], paths[args.healthy], model, list(X.columns),
        fractions=fractions, total_fragments=args.total_fragments,
        n_replicates=args.replicates, gc_track=gc_track, genome=args.genome,
    )

    tabdir = os.path.join(args.outdir, "tables")
    figdir = os.path.join(args.outdir, "figures")
    os.makedirs(tabdir, exist_ok=True)
    os.makedirs(figdir, exist_ok=True)
    series.to_csv(os.path.join(tabdir, "dilution_series.tsv"),
                  sep="\t", index=False)
    plots.plot_dilution(series, thresh,
                        os.path.join(figdir, "dilution_lod.png"))

    print()
    print(f"{'fraction':>10}{'mean score':>12}{'detected':>10}")
    for frac, grp in series.groupby("fraction"):
        rate = float((grp["score"] >= thresh).mean())
        print(f"{frac:>10.3f}{grp['score'].mean():>12.3f}"
              f"{rate:>9.0%}")

    lod = estimate_lod(series, thresh, min_detection_rate=0.95)

    # A population-threshold LOD on its own cannot distinguish "the feature
    # does not track tumour content" from "it tracks it precisely but the
    # between-person spread is wider than this donor's effect". The
    # dose-response fit separates them.
    dose = {}
    for feat in ("an_chr_entropy", "an_coverage_cv"):
        if feat in series.columns:
            try:
                dose[feat] = dose_response(series, feat)
            except (KeyError, ValueError):
                continue
    # Which way does each feature go across the whole cohort? A single donor
    # pair can sit on the wrong side of a real cohort-level effect.
    cohort_dir = {}
    if "aneuploidy" in blocks:
        an_all = blocks["aneuploidy"].reindex(order)
        lab = (train.set_index("sample_id").loc[order, "label"] == "cancer")
        for feat in an_all.columns:
            v = an_all[feat]
            cohort_dir[feat] = bool(v[lab.to_numpy()].median()
                                    > v[~lab.to_numpy()].median())

    if dose:
        print()
        print("dose-response")
        for feat, d in dose.items():
            print(f"  {feat}")
            print(f"    linear in fraction : R2 = {d['r2_linear']:.4f}, "
                  f"slope = {d['slope_per_unit_fraction']:+.3e}")
            print(f"    replicate SD       = {d['replicate_sd_at_zero']:.2e} "
                  f"({d['replicate_sd_as_pct_of_span']:.0f}% of the span)")
            print(f"    LOD vs own baseline= {d['lod_own_baseline_fraction']:.3f} "
                  f"plasma fraction at {d['n_sd']:.0f} SD")
        warnings = [w for w in
                    (check_direction(d, cohort_dir.get(f))
                     for f, d in dose.items()) if w]
        for w in warnings:
            print()
            print(f"  [warn] {w}")
        report_dir = {f: cohort_dir.get(f) for f in dose}

        concordant = {f: d for f, d in dose.items()
                      if cohort_dir.get(f) is None
                      or (d["slope_per_unit_fraction"] > 0) == cohort_dir[f]}
        if not concordant:
            print()
            print("  [warn] every feature in this series runs against the "
                  "cohort direction.")
            print("         Pick a different donor pair before quoting an LOD.")
            concordant = dose
        best = max(concordant.values(), key=lambda d: d["r2_linear"])
        if best["r2_linear"] > 0.95 and (np.isnan(lod) or lod >= 0.5):
            print()
            print(f"  {best['feature']} is linear in mixing fraction "
                  f"(R2 = {best['r2_linear']:.3f}) but was not detected "
                  f"against the")
            print(f"  population threshold until "
                  f"{'1.0' if lod == 1.0 else lod}. The feature is "
                  f"quantitative; what limits it is that")
            print(f"  between-person variation is larger than this donor's "
                  f"effect. Against the same")
            print(f"  donor's own baseline the limit is "
                  f"{best['lod_own_baseline_fraction']:.2f}, which is the "
                  f"monitoring case rather than screening.")
    out = {
        "cancer_donor": args.cancer, "healthy_donor": args.healthy,
        "blocks": use, "n_training_donors": len(order),
        "training_auc": baseline.auc_mean,
        "total_fragments_per_mixture": args.total_fragments,
        "replicates_per_level": args.replicates,
        "fractions_tested": list(fractions),
        "threshold_at_95_specificity": thresh,
        "lod95_plasma_fraction": None if np.isnan(lod) else lod,
        "dose_response": dose,
        "cohort_direction_case_higher": report_dir if dose else {},
        "note": ("LOD is a fraction of cancer-patient plasma, not tumour "
                 "fraction. The donor's plasma is itself mostly "
                 "haematopoietic cfDNA, so the true ctDNA LOD is lower by a "
                 "factor equal to that donor's tumour fraction, which this "
                 "cohort does not report. Treat this as an upper bound, and "
                 "as bounded below by the coarsest gap in the tested grid. "
                 "lod95_plasma_fraction is against a population threshold; "
                 "dose_response.lod_own_baseline_fraction is against the same "
                 "donor's own zero point, which is the monitoring case."),
    }
    with open(os.path.join(args.outdir, "dilution_summary.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print()
    print(json.dumps({k: v for k, v in out.items() if k != "note"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
