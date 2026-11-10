#!/usr/bin/env python3
"""Run the dilution series over several donor pairs and pool the result.

    python scripts/05_dilution_panel.py --features results/features \
        --cancer-type liver --genome hg19 --outdir results/dilution_panel

WHY A PANEL

A limit of detection from one donor pair is close to meaningless. In this
cohort the case-control difference in chromosome entropy varies about sixfold
between liver pairs, and the own-baseline LOD is inversely proportional to
that difference, so a single pair can report anything from optimistic to
hopeless depending on which one happened to be picked. Two donors at the same
stage can sit at opposite ends of that range, so it is not a stage effect
either; it is how much that particular tumour shifts the feature.

What this reports is the distribution: median and range across pairs, which
is the honest form of the claim.

Pairs whose slope runs against the cohort direction are run and reported but
excluded from the pooled figure, with a note. They are not failures; they are
donors for whom this feature does not point the right way, and their share of
the cohort is itself worth stating.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from cfdna_insight import modeling  # noqa: E402
from cfdna_insight.dilution import (  # noqa: E402
    check_direction, dilution_series, dose_response, estimate_lod,
    threshold_at_specificity,
)
from cfdna_insight.features import (  # noqa: E402
    compute_bin_gc, load_features, split_legacy_technical,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default="results/features")
    ap.add_argument("--outdir", default="results/dilution_panel")
    ap.add_argument("--cancer-type", default="liver",
                    help="run every matched pair of this cancer type")
    ap.add_argument("--pairs", default=None,
                    help="comma-separated pair_ids instead of --cancer-type")
    ap.add_argument("--feature", default="an_chr_entropy")
    ap.add_argument("--blocks", default="aneuploidy")
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--genome", default="hg19", choices=["hg19", "hg38"])
    ap.add_argument("--gc-cache", default="results/cache/bin_gc.parquet")
    ap.add_argument("--total-fragments", type=int, default=5_000_000)
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--fractions", default="0,0.01,0.05,0.1,0.25,0.5,1.0",
                    help="a coarser grid than the single-pair script, since "
                         "the panel is about the spread across donors")
    args = ap.parse_args()

    use = [b.strip() for b in args.blocks.split(",") if b.strip()]
    fractions = [float(f) for f in args.fractions.split(",")]

    blocks = split_legacy_technical(load_features(args.features))
    samples = pd.read_csv(os.path.join(args.features, "samples_used.tsv"),
                          sep="\t")
    paths = dict(zip(samples["sample_id"], samples["path"]))

    if args.pairs:
        pair_ids = [p.strip() for p in args.pairs.split(",")]
    else:
        pair_ids = sorted(samples.loc[samples["cancer_type"] ==
                                      args.cancer_type, "pair_id"].unique())
    if not pair_ids:
        print(f"[error] no pairs found", file=sys.stderr)
        return 1

    # cohort direction for each feature, so a discordant donor is recognised
    cohort_dir = {}
    if "aneuploidy" in blocks:
        an_all = blocks["aneuploidy"]
        lab = samples.set_index("sample_id").reindex(an_all.index)["label"]
        is_c = (lab == "cancer").to_numpy()
        for feat in an_all.columns:
            v = an_all[feat]
            cohort_dir[feat] = bool(v[is_c].median() > v[~is_c].median())

    gc_track = (compute_bin_gc(args.fasta, genome=args.genome,
                               cache=args.gc_cache) if args.fasta else None)

    print(f"[info] {len(pair_ids)} pairs: {', '.join(pair_ids)}")
    print(f"[info] {len(fractions)} fractions x {args.replicates} replicates "
          f"each\n")

    # An eight-pair panel is about an hour. Results are written as each pair
    # finishes and completed pairs are skipped on a re-run, so a crash or an
    # interrupt costs one pair rather than the whole run.
    tabdir = os.path.join(args.outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)
    summary_path = os.path.join(tabdir, "panel_summary.tsv")
    series_path = os.path.join(tabdir, "panel_series.tsv")

    rows, all_series = [], []
    if os.path.exists(summary_path):
        prev = pd.read_csv(summary_path, sep="\t")
        rows = prev.to_dict("records")
        done = set(prev["pair_id"])
        if os.path.exists(series_path):
            all_series.append(pd.read_csv(series_path, sep="\t"))
        pair_ids = [p for p in pair_ids if p not in done]
        print(f"[info] resuming: {len(done)} pairs already done, "
              f"{len(pair_ids)} to go\n")
    for i, pid in enumerate(pair_ids, 1):
        case, ctrl = pid, pid + "m"
        if case not in paths or ctrl not in paths:
            print(f"  [{i}/{len(pair_ids)}] {pid}: missing a member, skipped")
            continue
        if not (os.path.exists(paths[case]) and os.path.exists(paths[ctrl])):
            print(f"  [{i}/{len(pair_ids)}] {pid}: fragment file missing, "
                  f"skipped")
            continue

        train = samples[~samples["sample_id"].isin([case, ctrl])]
        order = [s for s in train["sample_id"]
                 if s in blocks["size_summary"].index]
        y = (train.set_index("sample_id").loc[order, "label"] == "cancer"
             ).astype(int).to_numpy()
        X = modeling.assemble_matrix(blocks, use, order)
        model = modeling.fit_final(X, y)
        base = modeling.nested_cv(X, y, n_repeats=3)
        thresh = threshold_at_specificity(base.oof_scores[y == 0], 0.95)

        series = dilution_series(
            paths[case], paths[ctrl], model, list(X.columns),
            fractions=fractions, total_fragments=args.total_fragments,
            n_replicates=args.replicates, gc_track=gc_track,
            genome=args.genome)
        series["pair_id"] = pid
        all_series.append(series)

        d = dose_response(series, args.feature)
        warn = check_direction(d, cohort_dir.get(args.feature))
        lod_pop = estimate_lod(series, thresh, 0.95)

        rows.append({
            "pair_id": pid,
            "span": d["span_over_series"],
            "r2": d["r2_linear"],
            "replicate_sd": d["replicate_sd_at_zero"],
            "lod_own_baseline": d["lod_own_baseline_fraction"],
            "lod_population": None if np.isnan(lod_pop) else lod_pop,
            "concordant": warn is None,
        })
        flag = "" if warn is None else "   [against cohort direction]"
        print(f"  [{i}/{len(pair_ids)}] {pid}: span {d['span_over_series']:+.2e}, "
              f"R2 {d['r2_linear']:.3f}, own-baseline LOD "
              f"{d['lod_own_baseline_fraction']:.3f}{flag}", flush=True)

        pd.DataFrame(rows).to_csv(summary_path, sep="\t", index=False)
        pd.concat(all_series, ignore_index=True).to_csv(
            series_path, sep="\t", index=False)

    if not rows:
        print("[error] no pairs ran", file=sys.stderr)
        return 1

    res = pd.DataFrame(rows)

    ok = res[res["concordant"]]
    print()
    print("=" * 62)
    print(f"{args.feature} across {len(res)} pairs "
          f"({len(ok)} concordant with the cohort direction)")
    print("=" * 62)
    if len(ok):
        lods = ok["lod_own_baseline"]
        print(f"  own-baseline LOD : median {lods.median():.3f}, "
              f"range {lods.min():.3f} to {lods.max():.3f}")
        print(f"  linearity        : R2 median {ok['r2'].median():.4f}, "
              f"min {ok['r2'].min():.4f}")
        pop = ok["lod_population"].dropna()
        print(f"  population LOD   : "
              f"{'all 1.0 (never detected below undiluted)' if (pop == 1.0).all() else f'median {pop.median():.3f}'}")
        print()
        print(f"  The spread across donors is the point. Effect size varies "
              f"{ok['span'].abs().max() / ok['span'].abs().min():.0f}x between "
              f"pairs,")
        print(f"  and the own-baseline LOD tracks it inversely. One pair "
              f"cannot stand for the cohort.")
    if (~res["concordant"]).any():
        bad = res.loc[~res["concordant"], "pair_id"].tolist()
        print()
        print(f"  excluded, slope against cohort direction: {', '.join(bad)} "
              f"({len(bad)}/{len(res)})")

    out = {
        "feature": args.feature,
        "n_pairs": int(len(res)),
        "n_concordant": int(len(ok)),
        "lod_own_baseline_median": float(ok["lod_own_baseline"].median())
        if len(ok) else None,
        "lod_own_baseline_min": float(ok["lod_own_baseline"].min())
        if len(ok) else None,
        "lod_own_baseline_max": float(ok["lod_own_baseline"].max())
        if len(ok) else None,
        "r2_median": float(ok["r2"].median()) if len(ok) else None,
        "per_pair": rows,
    }
    with open(os.path.join(args.outdir, "panel_summary.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\n[ok] wrote {args.outdir}/panel_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
