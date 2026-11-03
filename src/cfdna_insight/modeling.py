"""Evaluation machinery.

The claim this project has to defend is not "we got AUC 0.9". It is "we got
AUC X, and here is what X survives". So the evaluation is built around four
things that most cfDNA classifier write-ups omit:

1. Nested, repeated, group-aware CV. Hyperparameters are chosen on an inner
   loop so the reported AUC is not optimistically biased by model selection.

2. A permutation null. With n ~ 50 and p ~ 500, an AUC of 0.75 is not
   self-evidently meaningful. Shuffling labels and re-running the entire
   pipeline gives the null distribution it should be judged against.

3. Feature-block ablation. If size summaries alone match the full model, the
   500-dimensional DELFI profile is decoration.

4. A technical-covariate negative control. A model trained only on library
   depth, mapping quality and coverage dispersion should be at chance. If it
   isn't, the headline number is partly batch effect and must be reported that
   way.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

RANDOM_STATE = 20260831


# --------------------------------------------------------------------------
# Estimators
# --------------------------------------------------------------------------

def make_estimator(kind: str = "logistic") -> tuple[Pipeline, dict]:
    """Return an (estimator, param_grid) pair.

    Ridge logistic regression is the default rather than a gradient-boosted
    tree. With n in the tens and p in the hundreds, a strongly regularised
    linear model is both better calibrated and far less prone to fitting a
    single confounded window. GBM is available for comparison, not as default.
    """
    if kind == "logistic":
        pipe = Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("clf", LogisticRegression(solver="liblinear",
                                       max_iter=5000,
                                       random_state=RANDOM_STATE)),
        ])
        grid = {"clf__C": [0.001, 0.01, 0.1, 1.0]}
        return pipe, grid

    if kind == "gbm":
        pipe = Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("clf", GradientBoostingClassifier(random_state=RANDOM_STATE)),
        ])
        grid = {"clf__n_estimators": [100, 300],
                "clf__max_depth": [1, 2],
                "clf__learning_rate": [0.05, 0.1]}
        return pipe, grid

    raise ValueError(f"unknown estimator: {kind}")


# --------------------------------------------------------------------------
# Nested CV
# --------------------------------------------------------------------------

@dataclass
class CVResult:
    auc_per_repeat: np.ndarray
    oof_scores: np.ndarray
    y: np.ndarray
    label: str = ""
    extra: dict = field(default_factory=dict)

    @property
    def auc_mean(self) -> float:
        return float(np.mean(self.auc_per_repeat))

    @property
    def auc_ci(self) -> tuple[float, float]:
        """Percentile interval across repeats.

        This is a stability interval over resampling of the CV partition, not
        a confidence interval for the population AUC. With n=50 the true
        sampling interval is considerably wider; both are reported in the
        write-up so the distinction is not quietly lost.
        """
        return (float(np.percentile(self.auc_per_repeat, 2.5)),
                float(np.percentile(self.auc_per_repeat, 97.5)))

    def summary(self) -> dict:
        lo, hi = self.auc_ci
        return {"label": self.label, "auc_mean": self.auc_mean,
                "auc_lo": lo, "auc_hi": hi,
                "n": int(len(self.y)), "n_pos": int(self.y.sum()),
                **self.extra}


def nested_cv(
    X: pd.DataFrame,
    y: np.ndarray,
    estimator: str = "logistic",
    outer_folds: int = 5,
    inner_folds: int = 4,
    n_repeats: int = 10,
    label: str = "",
    random_state: int = RANDOM_STATE,
) -> CVResult:
    """Repeated stratified nested CV, returning per-repeat AUC and OOF scores."""
    X = X.copy()
    aucs, oof_accum = [], np.zeros((n_repeats, len(y)))

    for rep in range(n_repeats):
        outer = StratifiedKFold(n_splits=outer_folds, shuffle=True,
                                random_state=random_state + rep)
        oof = np.full(len(y), np.nan)

        for train_idx, test_idx in outer.split(X, y):
            pipe, grid = make_estimator(estimator)
            inner = StratifiedKFold(n_splits=inner_folds, shuffle=True,
                                    random_state=random_state + rep)
            search = GridSearchCV(pipe, grid, scoring="roc_auc", cv=inner,
                                  n_jobs=1, refit=True)
            search.fit(X.iloc[train_idx], y[train_idx])
            oof[test_idx] = search.predict_proba(X.iloc[test_idx])[:, 1]

        aucs.append(roc_auc_score(y, oof))
        oof_accum[rep] = oof

    return CVResult(np.array(aucs), oof_accum.mean(axis=0), y, label=label)


def cross_study_cv(X: pd.DataFrame, y: np.ndarray, groups: np.ndarray,
                   estimator: str = "logistic",
                   label: str = "leave-one-study-out") -> CVResult:
    """Leave-one-study-out evaluation.

    This is the number that actually matters for a screening assay. Within-
    cohort CV mostly measures whether a model can learn one lab's batch
    signature. Held-out-study AUC measures whether the signal transfers.
    Expect it to be substantially lower; report it anyway.
    """
    unique = np.unique(groups)
    if len(unique) < 2:
        raise ValueError("cross-study CV needs >= 2 studies in the sample sheet")

    oof = np.full(len(y), np.nan)
    per_study = {}
    for study in unique:
        test = groups == study
        train = ~test
        if len(np.unique(y[train])) < 2 or len(np.unique(y[test])) < 2:
            continue
        pipe, grid = make_estimator(estimator)
        inner = StratifiedKFold(n_splits=4, shuffle=True,
                                random_state=RANDOM_STATE)
        search = GridSearchCV(pipe, grid, scoring="roc_auc", cv=inner, n_jobs=1)
        search.fit(X[train], y[train])
        oof[test] = search.predict_proba(X[test])[:, 1]
        per_study[str(study)] = float(roc_auc_score(y[test], oof[test]))

    valid = np.isfinite(oof)
    auc = roc_auc_score(y[valid], oof[valid]) if valid.sum() else np.nan
    return CVResult(np.array([auc]), oof, y, label=label,
                    extra={"per_study_auc": per_study})


# --------------------------------------------------------------------------
# Permutation null
# --------------------------------------------------------------------------

def permutation_null(X: pd.DataFrame, y: np.ndarray, n_perm: int = 200,
                     estimator: str = "logistic", n_repeats: int = 2,
                     random_state: int = RANDOM_STATE
                     ) -> tuple[np.ndarray, float]:
    """Null AUC distribution under label shuffling, plus a one-sided p-value.

    The full CV is re-run inside each permutation. Permuting labels and reusing
    a model fit on the true labels would test nothing.
    """
    rng = np.random.default_rng(random_state)
    observed = nested_cv(X, y, estimator=estimator,
                         n_repeats=n_repeats).auc_mean

    null = np.empty(n_perm)
    for i in range(n_perm):
        y_perm = rng.permutation(y)
        null[i] = nested_cv(X, y_perm, estimator=estimator, n_repeats=1,
                            random_state=random_state + i).auc_mean

    p = (np.sum(null >= observed) + 1) / (n_perm + 1)
    return null, float(p)


# --------------------------------------------------------------------------
# Feature blocks and ablation
# --------------------------------------------------------------------------

def assemble_matrix(blocks: dict[str, pd.DataFrame],
                    use: list[str],
                    sample_order: list[str],
                    min_frac_present: float = 0.8) -> pd.DataFrame:
    """Concatenate selected blocks, aligned on sample id.

    Two kinds of column are dropped:

      - constant ones, which carry no information;
      - sparsely observed ones, present in fewer than `min_frac_present` of
        samples.

    The second matters more than it looks. `build_cnv_block` takes the union
    of chromosome arms across samples, so an arm summarised for only a
    handful of them becomes a column that is entirely missing inside any CV
    fold where none of those samples land in training. The imputer then has
    nothing to impute from, warns, and silently drops the column in some
    folds but not others — which makes the feature set quietly
    fold-dependent. Dropping them up front keeps every fold working from the
    same columns.
    """
    parts = []
    for name in use:
        if name not in blocks:
            continue
        parts.append(blocks[name].reindex(sample_order))
    if not parts:
        raise ValueError(f"none of {use} present in blocks {list(blocks)}")
    X = pd.concat(parts, axis=1)

    informative = X.nunique(dropna=True) > 1
    present = X.notna().mean(axis=0) >= min_frac_present
    keep = X.columns[informative & present]
    if len(keep) == 0:
        raise ValueError(
            f"no usable features from {use}: {int((~informative).sum())} "
            f"constant, {int((~present).sum())} too sparsely observed")
    return X[keep]


def ablation(blocks: dict[str, pd.DataFrame], y: np.ndarray,
             sample_order: list[str],
             block_sets: dict[str, list[str]] | None = None,
             estimator: str = "logistic",
             n_repeats: int = 10) -> pd.DataFrame:
    """Run nested CV on each feature-block combination."""
    block_sets = block_sets or {
        "size summary only": ["size_summary"],
        "size histogram only": ["size_hist"],
        "DELFI profile only": ["delfi"],
        "CNV arms only": ["cnv"],
        "aneuploidy proxies only": ["aneuploidy"],
        "end motifs only": ["motif"],
        "all biological": ["size_summary", "size_hist", "delfi", "cnv",
                           "motif", "aneuploidy"],
        "TECHNICAL ONLY (negative control)": ["technical"],
    }
    rows = []
    for name, use in block_sets.items():
        available = [b for b in use if b in blocks]
        if not available:
            continue
        X = assemble_matrix(blocks, available, sample_order)
        res = nested_cv(X, y, estimator=estimator, n_repeats=n_repeats,
                        label=name)
        rows.append({**res.summary(), "n_features": X.shape[1],
                     "blocks": ",".join(available)})
    return pd.DataFrame(rows).sort_values("auc_mean", ascending=False)


def learning_curve(X: pd.DataFrame, y: np.ndarray,
                   fractions=(0.4, 0.6, 0.8, 1.0),
                   n_repeats: int = 5,
                   estimator: str = "logistic",
                   random_state: int = RANDOM_STATE) -> pd.DataFrame:
    """AUC as a function of training-set size.

    A curve that is still climbing steeply at n=all means the estimate is
    sample-size limited and the headline AUC should not be quoted as if it
    were converged.
    """
    rng = np.random.default_rng(random_state)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    rows = []
    for frac in fractions:
        n_pos = int(round(len(pos) * frac))
        n_neg = int(round(len(neg) * frac))
        # need enough of each class for the inner 4-fold CV to be defined
        if min(n_pos, n_neg) < 4:
            continue
        aucs = []
        for r in range(n_repeats):
            idx = np.concatenate([rng.choice(pos, n_pos, replace=False),
                                  rng.choice(neg, n_neg, replace=False)])
            res = nested_cv(X.iloc[idx], y[idx], estimator=estimator,
                            n_repeats=1, random_state=random_state + r)
            aucs.append(res.auc_mean)
        if aucs:
            rows.append({"fraction": frac, "n": n_pos + n_neg,
                         "n_pos": n_pos, "n_neg": n_neg,
                         "auc_mean": float(np.mean(aucs)),
                         "auc_sd": float(np.std(aucs))})
    return pd.DataFrame(rows)


def fit_final(X: pd.DataFrame, y: np.ndarray,
              estimator: str = "logistic") -> Pipeline:
    """Fit on everything, for use as a locked scorer in the dilution series."""
    pipe, grid = make_estimator(estimator)
    search = GridSearchCV(pipe, grid, scoring="roc_auc",
                          cv=StratifiedKFold(5, shuffle=True,
                                             random_state=RANDOM_STATE))
    search.fit(X, y)
    return search.best_estimator_


def coefficient_report(model: Pipeline, feature_names: list[str],
                       top_n: int = 30) -> pd.DataFrame:
    """Standardised coefficients, largest magnitude first."""
    clf = model.named_steps["clf"]
    if not hasattr(clf, "coef_"):
        raise TypeError("coefficient report requires a linear model")
    coefs = clf.coef_.ravel()
    df = pd.DataFrame({"feature": feature_names, "coef": coefs})
    df["abs_coef"] = df["coef"].abs()
    return df.sort_values("abs_coef", ascending=False).head(top_n)


# --------------------------------------------------------------------------
# Fold-safe CNV evaluation
# --------------------------------------------------------------------------

def make_foldsafe_cnv_pipeline(window_arms, n_coverage_cols: int,
                               estimator: str = "logistic"):
    """Pipeline that rebuilds the healthy CNV reference inside every fold.

    The first `n_coverage_cols` columns of X are raw per-window coverage and
    go through `HealthyReferenceTransformer`; anything after them is already a
    finished feature and passes through untouched. That split lets the
    combined model stay fold-safe without precomputing CNV features from the
    whole cohort.
    """
    from sklearn.compose import ColumnTransformer

    from .cnv import HealthyReferenceTransformer

    cov_idx = list(range(n_coverage_cols))
    rest_idx = None  # filled per-call below

    def build(total_cols: int):
        rest = list(range(n_coverage_cols, total_cols))
        steps = [("cnv", HealthyReferenceTransformer(window_arms), cov_idx)]
        if rest:
            steps.append(("passthrough", "passthrough", rest))
        pre = ColumnTransformer(steps, remainder="drop")
        pipe, grid = make_estimator(estimator)
        full = Pipeline([("pre", pre)] + list(pipe.steps))
        return full, grid

    return build


def nested_cv_foldsafe_cnv(
    coverage: pd.DataFrame,
    window_arms,
    y: np.ndarray,
    other: pd.DataFrame | None = None,
    estimator: str = "logistic",
    outer_folds: int = 5,
    inner_folds: int = 4,
    n_repeats: int = 10,
    label: str = "",
    random_state: int = RANDOM_STATE,
) -> CVResult:
    """Nested CV with the CNV reference rebuilt per training fold.

    Slower than `nested_cv` on precomputed CNV features, because the reference
    and every arm summary are recomputed for each fold. That cost is the whole
    point: the difference between this number and the precomputed one is the
    size of the leak.
    """
    X = coverage if other is None else pd.concat([coverage, other], axis=1)
    X = X.to_numpy(dtype=float)
    n_cov = coverage.shape[1]
    build = make_foldsafe_cnv_pipeline(window_arms, n_cov, estimator)

    aucs, oof_accum = [], np.zeros((n_repeats, len(y)))
    for rep in range(n_repeats):
        outer = StratifiedKFold(n_splits=outer_folds, shuffle=True,
                                random_state=random_state + rep)
        oof = np.full(len(y), np.nan)
        for train_idx, test_idx in outer.split(X, y):
            pipe, grid = build(X.shape[1])
            inner = StratifiedKFold(n_splits=inner_folds, shuffle=True,
                                    random_state=random_state + rep)
            search = GridSearchCV(pipe, grid, scoring="roc_auc", cv=inner,
                                  n_jobs=1, refit=True)
            search.fit(X[train_idx], y[train_idx])
            oof[test_idx] = search.predict_proba(X[test_idx])[:, 1]
        aucs.append(roc_auc_score(y, oof))
        oof_accum[rep] = oof

    return CVResult(np.array(aucs), oof_accum.mean(axis=0), y, label=label)


def null_pvalue(X: pd.DataFrame, y: np.ndarray, observed: float,
                n_perm: int = 100, estimator: str = "logistic",
                random_state: int = RANDOM_STATE) -> tuple[float, float, float]:
    """One-sided p for an already-computed AUC, plus the null's mean and sd.

    Used to judge the technical-covariate control. A fixed AUC threshold is
    the wrong test at small n: with 8 per class the null AUC has a standard
    deviation near 0.15, so 0.75 is unremarkable, while at n=200 the same
    value would be damning.
    """
    rng = np.random.default_rng(random_state)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = nested_cv(X, rng.permutation(y), estimator=estimator,
                            n_repeats=1, random_state=random_state + i).auc_mean
    p = (np.sum(null >= observed) + 1) / (n_perm + 1)
    return float(p), float(null.mean()), float(null.std())
