#!/usr/bin/env python3
"""Train, evaluate, ablate, permute.

    python scripts/02_train_evaluate.py --features results/features

Produces every table and figure the write-up needs, including the ones that
make the result look worse. That is the point.
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
from cfdna_insight.cnv import HealthyReference, arm_summary  # noqa: E402
from cfdna_insight.features import load_features, split_legacy_technical  # noqa: E402


def build_cnv_block(coverage: pd.DataFrame, window_arms: pd.DataFrame,
                    y: np.ndarray) -> pd.DataFrame:
    """Arm-level CNV z-scores.

    NOTE ON LEAKAGE: the healthy reference is built here from all healthy
    samples, including ones that will land in test folds. For arm-level medians
    across dozens of samples the effect is small, but it is not zero, and the
    honest version recomputes the reference inside each training fold. The
    fold-safe path is `--cnv-fold-safe`, which is slower; both numbers are
    reported so the size of the difference is visible rather than assumed.
    """
    healthy_cov = coverage[y == 0]
    ref = HealthyReference.fit(healthy_cov)
    arms = window_arms.iloc[0]
    rows = {sid: arm_summary(ref.zscore(coverage.loc[sid]), arms)
            for sid in coverage.index}
    return pd.DataFrame(rows).T


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default="results/features")
    ap.add_argument("--grouping", default="study",
                    help="What the `study` column holds, for labelling only. "
                         "Pass 'cancer-type' with a by-cancer-type sheet so "
                         "the output and summary.json do not claim a "
                         "cross-study result that was never run.")
    ap.add_argument("--samples", default=None,
                    help="Sample sheet to use instead of the features "
                         "directory's own samples_used.tsv. The point of the "
                         "override is the `study` column: pass "
                         "config/samples_by_cancer_type.tsv to run "
                         "leave-one-cancer-type-out instead of "
                         "leave-one-study-out.")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--estimator", default="logistic",
                    choices=["logistic", "gbm"])
    ap.add_argument("--n-repeats", type=int, default=10)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--skip-permutation", action="store_true")
    ap.add_argument("--cnv-fold-safe", action="store_true",
                    help="Also rebuild the CNV healthy reference inside each "
                         "CV fold and report it next to the precomputed "
                         "number. Slower; the gap between them is the leak.")
    args = ap.parse_args()

    tabdir = os.path.join(args.outdir, "tables")
    figdir = os.path.join(args.outdir, "figures")
    os.makedirs(tabdir, exist_ok=True)
    os.makedirs(figdir, exist_ok=True)

    blocks = split_legacy_technical(load_features(args.features))
    sheet_path = args.samples or os.path.join(args.features,
                                              "samples_used.tsv")
    samples = pd.read_csv(sheet_path, sep="\t")
    if args.samples:
        print(f"[info] sample sheet: {sheet_path}")
        if args.grouping == "study":
            print("[warn] --samples given but --grouping left as 'study'. "
                  "If the sheet groups by something else (cancer type, "
                  "batch), pass --grouping so the output does not report a "
                  "cross-study result that was not run.")
    order = [s for s in samples["sample_id"]
             if s in blocks["size_summary"].index]
    samples = samples.set_index("sample_id").loc[order].reset_index()
    y = (samples["label"] == "cancer").astype(int).to_numpy()
    groups = samples["study"].to_numpy()

    print(f"[info] {len(order)} samples, {y.sum()} cancer, "
          f"{len(y) - y.sum()} non-cancer, "
          f"{len(np.unique(groups))} {args.grouping} group(s)")

    if "coverage" in blocks and "window_arms" in blocks:
        blocks["cnv"] = build_cnv_block(
            blocks["coverage"].reindex(order),
            blocks["window_arms"].reindex(order), y)
        print(f"[info] CNV block: {blocks['cnv'].shape}")

    # ---------------------------------------------------------------- ablation
    print("[run] feature-block ablation")
    abl = modeling.ablation(blocks, y, order, estimator=args.estimator,
                            n_repeats=args.n_repeats)
    abl.to_csv(os.path.join(tabdir, "ablation.tsv"), sep="\t", index=False)
    print(abl[["label", "n_features", "auc_mean", "auc_lo", "auc_hi"]]
          .to_string(index=False))
    plots.plot_ablation(abl, os.path.join(figdir, "ablation_auc.png"))

    # The technical control is judged against its own permutation null, not a
    # fixed AUC. At 8 per class the null AUC has sd near 0.15, so 0.75 means
    # nothing; at n=200 it would be damning. A constant threshold gets this
    # backwards at both ends.
    tech = abl[abl["label"].str.contains("TECHNICAL")]
    tech_p = None
    if len(tech) and "technical" in blocks:
        obs = float(tech["auc_mean"].iloc[0])
        Xt = modeling.assemble_matrix(blocks, ["technical"], order)
        n_perm_tech = max(30, args.n_perm // 4)
        print(f"[run] technical-control null ({n_perm_tech} permutations)")
        tech_p, null_mu, null_sd = modeling.null_pvalue(
            Xt, y, obs, n_perm=n_perm_tech, estimator=args.estimator)
        print(f"       technical AUC {obs:.3f}; null {null_mu:.3f} "
              f"+/- {null_sd:.3f}; p = {tech_p:.3g}")
        if tech_p < 0.05:
            print(f"\n[WARNING] the technical-covariate model beats its own "
                  f"null (p={tech_p:.3g}). Library characteristics alone "
                  f"separate the groups; the biological AUC is confounded "
                  f"and must be reported with that caveat.\n")
        else:
            print(f"       not distinguishable from chance at this n -- no "
                  f"evidence of depth/coverage confounding\n")

    # ------------------------------------------------------------- main model
    bio = [b for b in ("size_summary", "size_hist", "delfi", "cnv", "motif",
                       "aneuploidy") if b in blocks]
    X = modeling.assemble_matrix(blocks, bio, order)
    print(f"[run] full model: {X.shape[1]} features")
    main_res = modeling.nested_cv(X, y, estimator=args.estimator,
                                  n_repeats=args.n_repeats, label="full model")
    lo, hi = main_res.auc_ci
    print(f"       AUC {main_res.auc_mean:.3f} [{lo:.3f}, {hi:.3f}]")

    # -------------------------------------------------------- learning curve
    print("[run] learning curve")
    lc = modeling.learning_curve(X, y, estimator=args.estimator)
    lc.to_csv(os.path.join(tabdir, "learning_curve.tsv"), sep="\t", index=False)
    print(lc.to_string(index=False))

    # --------------------------------------------------------- cross-study
    cross = None
    if len(np.unique(groups)) >= 2:
        print(f"[run] leave-one-{args.grouping}-out")
        try:
            cross = modeling.cross_study_cv(
                X, y, groups, estimator=args.estimator,
                label=f"leave-one-{args.grouping}-out")
            print(f"       held-out-{args.grouping} AUC "
                  f"{cross.auc_mean:.3f}")
            print(f"       per {args.grouping}: "
                  f"{cross.extra['per_study_auc']}")
        except ValueError as exc:
            print(f"       skipped: {exc}")
    else:
        print("[info] single study; cross-study generalisation UNTESTED. "
              "Any AUC below is within-cohort only.")

    # ------------------------------------------------------ fold-safe CNV
    foldsafe = {}
    if args.cnv_fold_safe and "coverage" in blocks and "cnv" in blocks:
        print("[run] fold-safe CNV (reference rebuilt per training fold)")
        cov = blocks["coverage"].reindex(order)
        arms = blocks["window_arms"].reindex(order).iloc[0]
        try:
            fs_cnv = modeling.nested_cv_foldsafe_cnv(
                cov, arms, y, estimator=args.estimator,
                n_repeats=max(3, args.n_repeats // 2),
                label="CNV arms only (fold-safe)")
            leaky = float(abl.loc[abl["label"] == "CNV arms only",
                                  "auc_mean"].iloc[0])
            foldsafe["cnv_only_leaky"] = leaky
            foldsafe["cnv_only_foldsafe"] = fs_cnv.auc_mean
            foldsafe["cnv_leak_delta"] = leaky - fs_cnv.auc_mean
            print(f"       CNV precomputed reference: {leaky:.3f}")
            print(f"       CNV fold-safe reference:   {fs_cnv.auc_mean:.3f}")
            print(f"       leak: {leaky - fs_cnv.auc_mean:+.3f}")
            if leaky - fs_cnv.auc_mean > 0.1:
                print(f"\n[WARNING] the precomputed CNV reference inflates "
                      f"AUC by {leaky - fs_cnv.auc_mean:.3f}. Report the "
                      f"fold-safe number.\n")
        except ValueError as exc:
            print(f"       skipped: {exc}")

    # ---------------------------------------------------------- permutation
    perm_p = None
    if not args.skip_permutation:
        print(f"[run] permutation null ({args.n_perm} permutations)")
        null, perm_p = modeling.permutation_null(
            X, y, n_perm=args.n_perm, estimator=args.estimator)
        print(f"       null mean {null.mean():.3f}, p = {perm_p:.4g}")
        np.savetxt(os.path.join(tabdir, "permutation_null.txt"), null)
        plots.plot_permutation_null(
            null, main_res.auc_mean, perm_p,
            os.path.join(figdir, "permutation_null.png"))

    # -------------------------------------------------------------- figures
    if "size_hist" in blocks:
        plots.plot_size_distributions(
            blocks["size_hist"].reindex(order), y,
            os.path.join(figdir, "size_distributions.png"))
    if "delfi" in blocks:
        plots.plot_delfi_heatmap(blocks["delfi"].reindex(order), y,
                                 os.path.join(figdir, "delfi_heatmap.png"))

    # -------------------------------------------------------- final + coefs
    final = modeling.fit_final(X, y, estimator=args.estimator)
    if args.estimator == "logistic":
        coefs = modeling.coefficient_report(final, list(X.columns))
        coefs.to_csv(os.path.join(tabdir, "top_coefficients.tsv"),
                     sep="\t", index=False)

    scores_healthy = main_res.oof_scores[y == 0]
    thresh = float(np.quantile(scores_healthy, 0.95))
    sens = float(np.mean(main_res.oof_scores[y == 1] >= thresh))

    summary = {
        "n_samples": len(order), "n_cancer": int(y.sum()),
        "n_groups": int(len(np.unique(groups))),
        "n_features": int(X.shape[1]),
        "estimator": args.estimator,
        "auc_mean": main_res.auc_mean,
        "auc_interval": list(main_res.auc_ci),
        "permutation_p": perm_p,
        "grouping": args.grouping,
        f"cross_{args.grouping.replace('-', '_')}_auc":
            cross.auc_mean if cross else None,
        f"per_{args.grouping.replace('-', '_')}_auc":
            cross.extra["per_study_auc"] if cross else None,
        "threshold_at_95_specificity": thresh,
        "sensitivity_at_95_specificity": sens,
        "technical_control_auc": float(tech["auc_mean"].iloc[0])
        if len(tech) else None,
        "technical_control_p": tech_p,
        **foldsafe,
    }
    with open(os.path.join(args.outdir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    pd.DataFrame({"sample_id": order, "label": samples["label"],
                  "study": groups, "oof_score": main_res.oof_scores}
                 ).to_csv(os.path.join(tabdir, "oof_scores.tsv"),
                          sep="\t", index=False)

    print("\n" + json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
