"""Figures. Every plot here is one that has to appear in the write-up."""

from __future__ import annotations

import numpy as np
import pandas as pd


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_size_distributions(size_hist: pd.DataFrame, y: np.ndarray,
                            out: str, lo: int = 50) -> str:
    """Mean size density by group, with the 167 bp mononucleosome peak marked."""
    plt = _mpl()
    x = np.arange(lo, lo + size_hist.shape[1])
    fig, ax = plt.subplots(figsize=(7, 4))
    for val, name, colour in ((0, "non-cancer", "#3b7dd8"),
                              (1, "cancer", "#d1495b")):
        sub = size_hist.to_numpy()[y == val]
        if sub.size == 0:
            continue
        m, sd = sub.mean(axis=0), sub.std(axis=0)
        ax.plot(x, m, color=colour, label=f"{name} (n={sub.shape[0]})")
        ax.fill_between(x, m - sd, m + sd, color=colour, alpha=0.18, lw=0)
    ax.axvline(167, ls="--", lw=0.8, color="0.4")
    ax.annotate("167 bp", (167, ax.get_ylim()[1] * 0.95), fontsize=8,
                ha="left", color="0.3")
    ax.set_xlabel("fragment length (bp)")
    ax.set_ylabel("density")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def plot_delfi_heatmap(delfi: pd.DataFrame, y: np.ndarray, out: str) -> str:
    """Samples x 5 Mb windows, cancer and non-cancer blocked separately."""
    plt = _mpl()
    order = np.argsort(y, kind="stable")
    mat = delfi.to_numpy()[order]
    fig, ax = plt.subplots(figsize=(10, 4))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3,
                   interpolation="nearest")
    split = int(np.sum(y == 0))
    ax.axhline(split - 0.5, color="k", lw=1)
    ax.set_xlabel("5 Mb window (genome order)")
    ax.set_ylabel("sample")
    ax.set_yticks([split / 2, split + (len(y) - split) / 2])
    ax.set_yticklabels(["non-cancer", "cancer"])
    fig.colorbar(im, ax=ax, label="z-scored short/long ratio", shrink=0.8)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def plot_permutation_null(null: np.ndarray, observed: float, p: float,
                          out: str) -> str:
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.hist(null, bins=30, color="0.75", edgecolor="white")
    ax.axvline(observed, color="#d1495b", lw=2)
    ax.annotate(f"observed = {observed:.3f}\npermutation p = {p:.3g}",
                (observed, ax.get_ylim()[1] * 0.8), fontsize=9,
                ha="right" if observed > np.median(null) else "left")
    ax.set_xlabel("AUC under shuffled labels")
    ax.set_ylabel("permutations")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def plot_ablation(table: pd.DataFrame, out: str) -> str:
    """Horizontal AUC bars, negative control highlighted."""
    plt = _mpl()
    t = table.sort_values("auc_mean")
    colours = ["#c47f2a" if "NEGATIVE" in str(l).upper() else "#3b7dd8"
               for l in t["label"]]
    fig, ax = plt.subplots(figsize=(7, 0.45 * len(t) + 1.5))
    err = np.vstack([t["auc_mean"] - t["auc_lo"], t["auc_hi"] - t["auc_mean"]])
    ax.barh(t["label"], t["auc_mean"], xerr=err, color=colours,
            error_kw={"lw": 1, "ecolor": "0.3"})
    ax.axvline(0.5, ls="--", lw=0.9, color="0.4")
    ax.set_xlim(0.3, 1.0)
    ax.set_xlabel("nested CV AUC")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def plot_dilution(series: pd.DataFrame, threshold: float, out: str) -> str:
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(6.5, 4))
    grouped = series.groupby("fraction")["score"]
    x = np.array(sorted(series["fraction"].unique()))
    ax.plot(x, grouped.mean().reindex(x), "o-", color="#3b7dd8")
    ax.fill_between(x, grouped.min().reindex(x), grouped.max().reindex(x),
                    color="#3b7dd8", alpha=0.18, lw=0)
    ax.axhline(threshold, ls="--", color="#d1495b", lw=1,
               label="threshold at 95% specificity")
    ax.set_xscale("symlog", linthresh=0.005)
    ax.set_xlabel("fraction of cancer-patient plasma in mixture")
    ax.set_ylabel("model score")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def plot_paired_feature(values: pd.Series, pair_id: pd.Series,
                        label: pd.Series, out: str,
                        feature_name: str = "", ylabel: str = "") -> str:
    """One line per matched pair, control on the left and case on the right.

    This is the right picture for a matched design and a scatter or boxplot is
    the wrong one: the comparison that the study supports is within pair, and
    between-pair spread is exactly the variance the matching was meant to
    remove. A feature can separate cleanly pair by pair while the two group
    distributions overlap almost completely, and only this plot shows that.
    """
    plt = _mpl()
    df = pd.DataFrame({"v": values, "pair": pair_id, "label": label})
    piv = df.pivot_table(index="pair", columns="label", values="v").dropna()
    if piv.empty or "cancer" not in piv or "healthy" not in piv:
        raise ValueError("need both cancer and healthy per pair")

    up = piv["cancer"] > piv["healthy"]
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    for _, r in piv.iterrows():
        rising = r["cancer"] > r["healthy"]
        ax.plot([0, 1], [r["healthy"], r["cancer"]],
                color="#d1495b" if rising else "#3b7dd8",
                alpha=0.65, lw=1.2, marker="o", ms=3.5,
                markerfacecolor="white", markeredgewidth=1.0)

    ax.set_xlim(-0.25, 1.25)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["matched\ncontrol", "case"])
    ax.set_ylabel(ylabel or feature_name or "value")
    ax.set_title(f"{feature_name}\n{int(up.sum())} of {len(piv)} pairs higher "
                 f"in the case", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out
