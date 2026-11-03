#!/usr/bin/env python3
"""Work out what the technical-covariate control is actually detecting.

    python scripts/04_confounder_audit.py --features results/features

WHY THIS EXISTS

On the full CRAG cohort the technical-covariate model reached AUC 0.886
(p=0.020) — higher than the 1191-feature biological model at 0.850. Something
about the libraries separates cases from controls, and until it is named, no
biological number from this cohort means anything.

Batch is not the explanation. Table S2 shows every case and its matched
control went through the same library-prep and sequencing batch, so batch is
orthogonal to label by construction. Whatever this is, it operates *within*
pair.

WHAT IT CHECKS

1. Which covariate carries the signal, by univariate AUC.

2. Whether it holds within matched pairs. Each case has one control matched on
   age, sex, smoking and alcohol and run in the same batch, so a paired test
   removes everything the matching controlled for. A difference that survives
   is a property of the sample, not of cohort assembly.

3. Whether the biological blocks are just reading something else. Each
   block's cross-validated AUC is recomputed after regressing out, in turn,
   log sequencing depth and chromosome entropy. A block that collapses was
   measuring that covariate; one that holds was not.

   The second column matters as much as the first. End motifs are read off
   the reference genome at fragment ends, so a sample whose coverage has been
   redistributed by arm-level gains and losses samples a different mix of
   sequence context and shifts its 4-mer spectrum with no change in nuclease
   activity at all. Motif and aneuploidy scoring highly in the same cohort is
   therefore not automatically two independent findings.

HOW TO READ IT

A depth difference between cases and controls is not automatically an
artefact. Cancer patients have higher plasma cfDNA concentration, and at fixed
library input that changes library complexity and therefore unique fragment
count. That is a real biological difference showing up through a technical
channel — which is still a confound for *this* feature set, because nothing
here measures it deliberately, but it is a different claim from "the batches
differed" and should be written up as such.
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
from cfdna_insight.features import load_features, split_legacy_technical  # noqa: E402

BIO_BLOCKS = ["size_summary", "size_hist", "delfi", "motif", "cnv",
              "aneuploidy"]


def residualise(X: pd.DataFrame, covariate: np.ndarray) -> pd.DataFrame:
    """Remove the linear component of `covariate` from every column.

    Ordinary least squares, column by column: each feature is replaced by its
    residual about its own regression on the covariate. Only the linear part
    is removed, so a feature coupled to depth non-linearly survives this and
    the result is a lower bound on depth's contribution, not an exact one.
    """
    c = np.asarray(covariate, dtype=float)
    c = (c - c.mean()) / (c.std() + 1e-12)
    design = np.column_stack([np.ones_like(c), c])
    out = X.copy()
    for col in X.columns:
        v = X[col].to_numpy(dtype=float)
        ok = np.isfinite(v)
        if ok.sum() < 5 or len(np.unique(v[ok])) < 2:
            continue
        beta, *_ = np.linalg.lstsq(design[ok], v[ok], rcond=None)
        resid = v.copy()
        resid[ok] = v[ok] - design[ok] @ beta
        out[col] = resid
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default="results/features")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--n-repeats", type=int, default=5)
    args = ap.parse_args()

    from scipy.stats import mannwhitneyu, wilcoxon
    from sklearn.metrics import roc_auc_score

    blocks = split_legacy_technical(load_features(args.features))
    samples = pd.read_csv(os.path.join(args.features, "samples_used.tsv"),
                          sep="\t")
    order = [s for s in samples["sample_id"]
             if s in blocks["size_summary"].index]
    samples = samples.set_index("sample_id").loc[order].reset_index()
    y = (samples["label"] == "cancer").astype(int).to_numpy()
    tech = blocks["technical"].reindex(order)

    report: dict = {"n_samples": len(order), "n_cancer": int(y.sum())}

    # -------------------------------------------------- 1. which covariate
    print("=" * 70)
    print("1. UNIVARIATE SEPARATION BY TECHNICAL COVARIATE")
    print("=" * 70)
    uni = []
    for c in tech.columns:
        v = tech[c].to_numpy(dtype=float)
        finite = v[np.isfinite(v)]
        if len(np.unique(finite)) < 2:
            print(f"  {c:30s} constant — dropped before modelling")
            continue
        auc = roc_auc_score(y, v)
        _, p = mannwhitneyu(v[y == 1], v[y == 0])
        uni.append({"covariate": c, "auc": auc, "p": p,
                    "cancer_median": float(np.median(v[y == 1])),
                    "healthy_median": float(np.median(v[y == 0]))})
        print(f"  {c:30s} AUC={auc:.3f}  p={p:<10.3g} "
              f"cancer={np.median(v[y == 1]):<12.4g} "
              f"healthy={np.median(v[y == 0]):.4g}")
    report["univariate"] = uni

    # ------------------------------------------------------ 2. within pair
    print()
    print("=" * 70)
    print("2. WITHIN MATCHED PAIRS (same batch, matched donor)")
    print("=" * 70)
    paired = []
    if "pair_id" in samples.columns:
        d = tech.join(samples.set_index("sample_id")[["pair_id", "label"]])
        for c in tech.columns:
            piv = d.pivot_table(index="pair_id", columns="label",
                                values=c).dropna()
            if piv.empty or "cancer" not in piv or "healthy" not in piv:
                continue
            diff = piv["cancer"] - piv["healthy"]
            if np.allclose(diff, 0):
                continue
            try:
                _, p = wilcoxon(piv["cancer"], piv["healthy"])
            except ValueError:
                continue
            frac = float((diff > 0).mean())
            paired.append({"covariate": c, "n_pairs": int(len(piv)),
                           "frac_cancer_higher": frac, "p": float(p)})
            flag = "  <-- holds within pair" if p < 0.05 else ""
            print(f"  {c:30s} cancer>control in {frac:5.0%} of "
                  f"{len(piv):2d} pairs  p={p:.4g}{flag}")
    report["paired"] = paired

    # ------------------------------------- 3. biology after removing depth
    print()
    print("=" * 70)
    print("3. BIOLOGICAL BLOCKS, RAW AND WITH EACH COVARIATE REMOVED")
    print("=" * 70)
    covariates = {"depth": tech["tech_log10_n_fragments"].to_numpy(dtype=float)}
    if "aneuploidy" in blocks and "an_chr_entropy" in blocks["aneuploidy"]:
        covariates["aneuploidy"] = (blocks["aneuploidy"]
                                    .reindex(order)["an_chr_entropy"]
                                    .to_numpy(dtype=float))

    rows = []
    header = "  ".join(f"-{k}" for k in covariates)
    print(f"  {'block':16s} {'feat':>5s}   {'raw':>6s}   {header}")
    for name in BIO_BLOCKS + ["technical"]:
        if name not in blocks:
            continue
        try:
            X = modeling.assemble_matrix(blocks, [name], order)
        except ValueError:
            continue
        before = modeling.nested_cv(X, y, n_repeats=args.n_repeats).auc_mean
        row = {"block": name, "n_features": int(X.shape[1]), "auc": before}
        cells = []
        for cname, cvals in covariates.items():
            if name == "aneuploidy" and cname == "aneuploidy":
                cells.append("    -")        # don't residualise on itself
                continue
            after = modeling.nested_cv(residualise(X, cvals), y,
                                       n_repeats=args.n_repeats).auc_mean
            row[f"auc_minus_{cname}"] = after
            cells.append(f"{after:.3f}")
        rows.append(row)
        print(f"  {name:16s} {X.shape[1]:>5d}   {before:.3f}   "
              + "   ".join(f"{c:>6s}" for c in cells))
    report["adjusted"] = rows

    # Is the motif signal independent of aneuploidy, or downstream of it?
    # End motifs are read off the reference at fragment ends, so a sample
    # whose coverage is redistributed by arm-level gains and losses samples a
    # different mix of sequence context and shifts its 4-mer spectrum without
    # any change in nuclease activity. If motif AUC collapses once chromosome
    # entropy is partialled out, that is what happened.
    mrow = next((r for r in rows if r["block"] == "motif"), None)
    if mrow and "auc_minus_aneuploidy" in mrow:
        drop = mrow["auc"] - mrow["auc_minus_aneuploidy"]
        report["motif_independent_of_aneuploidy"] = bool(drop < 0.1)
        print()
        if drop > 0.1:
            print(f"  End-motif AUC falls {drop:.3f} once chromosome entropy "
                  f"is removed:\n  the motif spectrum is largely a "
                  f"consequence of where coverage sits,\n  not independent "
                  f"nuclease biology. Report them together, not as two "
                  f"findings.")
        else:
            print(f"  End-motif AUC holds ({mrow['auc']:.3f} -> "
                  f"{mrow['auc_minus_aneuploidy']:.3f}) with chromosome "
                  f"entropy removed:\n  the motif signal is not merely "
                  f"downstream of aneuploidy.")

    # Partialling out a covariate that is itself a consequence of the label
    # is conditioning on a mediator, and that can manufacture association
    # rather than remove it. Blocks whose AUC RISES under an adjustment are
    # the symptom. The negative control is the diagnostic: nothing about
    # removing aneuploidy can reveal genuine library signal, so if the
    # technical block climbs, every rise in that column is suspect and only
    # the falls should be read.
    for cname in covariates:
        key = f"auc_minus_{cname}"
        trow = next((r for r in rows if r["block"] == "technical"), None)
        if not trow or key not in trow:
            continue
        induced = trow[key] - trow["auc"]
        risers = [r["block"] for r in rows
                  if key in r and r[key] - r["auc"] > 0.1
                  and r["block"] != "technical"]
        if induced > 0.05 and risers:
            report[f"induced_association_{cname}"] = True
            print()
            print(f"  CAUTION: removing {cname} lifts the NEGATIVE CONTROL "
                  f"by {induced:+.3f}.")
            print(f"  {cname} is itself label-associated, so partialling it "
                  f"out conditions on a")
            print(f"  mediator and can manufacture association. Treat the "
                  f"rises in that column")
            print(f"  ({', '.join(risers)}) as artefacts, not as hidden "
                  f"signal. The falls are still")
            print(f"  informative: an adjustment cannot destroy signal that "
                  f"was never there.")

    # ------------------------------------------------------------ verdict
    print()
    print("=" * 70)
    rows_depth = {r["block"]: r for r in rows}
    depth_auc = next((r["auc"] for r in uni
                      if r["covariate"] == "tech_log10_n_fragments"), None)
    depth_paired = next((r for r in paired
                         if r["covariate"] == "tech_log10_n_fragments"), None)
    if depth_auc is not None and depth_auc > 0.65:
        print(f"Sequencing depth alone separates the groups (AUC "
              f"{depth_auc:.3f}).")
        if depth_paired and depth_paired["p"] < 0.05:
            print(f"It holds within matched pairs (p={depth_paired['p']:.3g}, "
                  f"cancer higher in "
                  f"{depth_paired['frac_cancer_higher']:.0%}), so it is not an "
                  f"artefact of how the cohort was assembled — the case and "
                  f"control libraries themselves differ.")
        else:
            print("It does NOT hold within matched pairs, which points at "
                  "cohort assembly rather than the samples.")
    survivors = [r for r in rows
                 if r["block"] != "technical"
                 and r.get("auc_minus_depth", r["auc"]) > 0.65]
    if survivors:
        print("\nBlocks still above 0.65 with depth removed: "
              + ", ".join(f"{r['block']} ({r.get('auc_minus_depth', r['auc']):.3f})"
                          for r in survivors))
    else:
        print("\nNo biological block stays above 0.65 once depth is removed. "
              "On this cohort the separation is a depth effect.")
    print("=" * 70)

    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, "confounder_audit.json")
    with open(out, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"\n[ok] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
