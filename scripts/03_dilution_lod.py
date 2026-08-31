#!/usr/bin/env python3
"""In silico dilution series and limit-of-detection estimate.

    python scripts/03_dilution_lod.py \
        --features results/features --cancer SAMPLE_A --healthy SAMPLE_B

Trains a locked model on all samples *except* the two donors used for mixing,
then scores synthetic mixtures. Excluding the donors matters: a model that has
seen the cancer donor will recognise it at any dilution for the wrong reason.
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
    dilution_series, estimate_lod, threshold_at_specificity,
)
from cfdna_insight.features import compute_bin_gc, load_features  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default="results/features")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--cancer", required=True, help="sample_id of cancer donor")
    ap.add_argument("--healthy", required=True, help="sample_id of healthy donor")
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--gc-cache", default="results/cache/bin_gc.parquet")
    ap.add_argument("--total-fragments", type=int, default=5_000_000)
    ap.add_argument("--replicates", type=int, default=5)
    args = ap.parse_args()

    blocks = load_features(args.features)
    samples = pd.read_csv(os.path.join(args.features, "samples_used.tsv"),
                          sep="\t")
    paths = dict(zip(samples["sample_id"], samples["path"]))
    for sid in (args.cancer, args.healthy):
        if sid not in paths:
            raise SystemExit(f"{sid} not in sample sheet")

    # hold out both donors from training
    train = samples[~samples["sample_id"].isin([args.cancer, args.healthy])]
    order = [s for s in train["sample_id"] if s in blocks["size_summary"].index]
    y = (train.set_index("sample_id").loc[order, "label"] == "cancer"
         ).astype(int).to_numpy()

    # only blocks recomputable from a mixture: no CNV (needs a healthy panel)
    use = [b for b in ("size_summary", "size_hist", "delfi") if b in blocks]
    X = modeling.assemble_matrix(blocks, use, order)
    print(f"[info] locked model on {len(order)} donors, {X.shape[1]} features")

    model = modeling.fit_final(X, y)
    baseline = modeling.nested_cv(X, y, n_repeats=5)
    thresh = threshold_at_specificity(baseline.oof_scores[y == 0], 0.95)
    print(f"[info] training AUC {baseline.auc_mean:.3f}, "
          f"threshold at 95% specificity = {thresh:.3f}")

    gc_track = compute_bin_gc(args.fasta, cache=args.gc_cache) \
        if args.fasta else None

    print(f"[run] mixing {args.cancer} into {args.healthy}")
    series = dilution_series(
        paths[args.cancer], paths[args.healthy], model, list(X.columns),
        total_fragments=args.total_fragments, n_replicates=args.replicates,
        gc_track=gc_track,
    )

    tabdir = os.path.join(args.outdir, "tables")
    figdir = os.path.join(args.outdir, "figures")
    os.makedirs(tabdir, exist_ok=True)
    os.makedirs(figdir, exist_ok=True)
    series.to_csv(os.path.join(tabdir, "dilution_series.tsv"),
                  sep="\t", index=False)
    plots.plot_dilution(series, thresh,
                        os.path.join(figdir, "dilution_lod.png"))

    lod = estimate_lod(series, thresh, min_detection_rate=0.95)
    out = {
        "cancer_donor": args.cancer, "healthy_donor": args.healthy,
        "total_fragments_per_mixture": args.total_fragments,
        "replicates_per_level": args.replicates,
        "threshold_at_95_specificity": thresh,
        "lod95_plasma_fraction": None if np.isnan(lod) else lod,
        "note": ("LOD is expressed as fraction of cancer-patient plasma, not "
                 "tumour fraction; true ctDNA LOD is lower by the donor's own "
                 "tumour fraction and this figure is an upper bound."),
    }
    with open(os.path.join(args.outdir, "dilution_summary.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
