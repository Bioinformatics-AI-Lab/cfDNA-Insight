"""Tests run entirely on simulated fragment files — no downloads required."""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from cfdna_insight import bins, cnv, fragmentomics as fx, modeling, simulate
from cfdna_insight.dilution import estimate_lod, mix_fragments
from cfdna_insight.features import extract_sample_features, stack_blocks
from cfdna_insight.io import SampleSheet, read_fragments


@pytest.fixture(scope="module")
def healthy_file(tmp_path_factory):
    p = tmp_path_factory.mktemp("frags") / "healthy.frag.gz"
    return simulate.simulate_sample(str(p), n_fragments=60_000,
                                    is_cancer=False, seed=1)


@pytest.fixture(scope="module")
def cancer_file(tmp_path_factory):
    p = tmp_path_factory.mktemp("frags") / "cancer.frag.gz"
    return simulate.simulate_sample(str(p), n_fragments=60_000,
                                    is_cancer=True, seed=2)


# ------------------------------------------------------------------------ io

def test_read_fragments_filters(healthy_file):
    df = read_fragments(healthy_file)
    assert len(df) > 0
    assert (df["length"] >= 50).all() and (df["length"] <= 1000).all()
    assert df["chrom"].astype(str).str.startswith("chr").all()
    assert (df["mapq"] >= 30).all()


def test_read_fragments_length_matches_coords(healthy_file):
    df = read_fragments(healthy_file)
    assert (df["length"] == df["end"] - df["start"]).all()


def test_sample_sheet_rejects_bad_label(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("sample_id\tpath\tlabel\tstudy\nA\tx\tmaybe\tS1\n")
    with pytest.raises(ValueError, match="label must be"):
        SampleSheet.load(str(p))


def test_sample_sheet_rejects_duplicates(tmp_path):
    p = tmp_path / "dup.tsv"
    p.write_text("sample_id\tpath\tlabel\tstudy\n"
                 "A\tx\tcancer\tS1\nA\ty\thealthy\tS1\n")
    with pytest.raises(ValueError, match="duplicate"):
        SampleSheet.load(str(p))


# ---------------------------------------------------------------------- bins

def test_bins_tile_without_overlap():
    b = bins.make_bins(bin_size=1_000_000)
    for _, sub in b.groupby("chrom", observed=True):
        sub = sub.sort_values("start")
        assert (sub["start"].to_numpy()[1:] == sub["end"].to_numpy()[:-1]).all()


def test_arm_assignment_splits_at_centromere():
    b = bins.make_bins(bin_size=1_000_000)
    c1 = b[b["chrom"] == "chr1"]
    cen_start = bins.HG38_CENTROMERE["chr1"][0]
    p_arm = c1[c1["arm"] == "chr1p"]
    assert p_arm["start"].max() < cen_start


def test_aggregate_drops_incomplete_windows():
    b = bins.make_bins(bin_size=100_000)
    b["n_total"] = 1
    win = bins.aggregate_to_windows(b, factor=50, value_cols=["n_total"])
    assert (win["n_total"] == 50).all()


# -------------------------------------------------------------- fragmentomics

def test_size_distribution_sums_to_one(healthy_file):
    d = fx.size_distribution(read_fragments(healthy_file))
    assert np.isclose(d.sum(), 1.0)


def test_cancer_fragments_are_shorter(healthy_file, cancer_file):
    h = fx.size_summary(read_fragments(healthy_file))
    c = fx.size_summary(read_fragments(cancer_file))
    assert c["size_median"] < h["size_median"]
    assert c["frac_below_150"] > h["frac_below_150"]


def test_periodicity_detected_in_simulated_data(healthy_file):
    score = fx._periodicity_score(read_fragments(healthy_file)["length"].to_numpy())
    assert np.isfinite(score) and score > 0


def test_delfi_profile_is_centred(healthy_file):
    binned = fx.bin_short_long(read_fragments(healthy_file))
    prof = fx.delfi_profile(binned, min_total=1)
    finite = prof[np.isfinite(prof)]
    assert len(finite) > 20
    assert abs(finite.mean()) < 0.2
    assert 0.7 < finite.std() < 1.4


def test_gc_correction_preserves_scale():
    rng = np.random.default_rng(0)
    gc = rng.uniform(0.3, 0.6, 2000)
    counts = 100 * (1 + 2 * (gc - 0.45)) + rng.normal(0, 5, 2000)
    corrected = fx.gc_correct(counts, gc)
    # correlation with GC should drop substantially
    before = abs(np.corrcoef(counts, gc)[0, 1])
    after = abs(np.corrcoef(corrected, gc)[0, 1])
    assert after < before


def test_mds_bounds():
    uniform = pd.Series(np.ones(256) / 256,
                        index=[f"motif_{m}" for m in fx.FOUR_MERS])
    assert np.isclose(fx.motif_diversity_score(uniform), 1.0)
    spike = pd.Series(np.zeros(256), index=uniform.index)
    spike.iloc[0] = 1.0
    assert fx.motif_diversity_score(spike) < 0.01


# ----------------------------------------------------------------------- cnv

def test_healthy_reference_zscore_centred():
    rng = np.random.default_rng(0)
    mat = pd.DataFrame(rng.normal(1.0, 0.05, (30, 100)))
    ref = cnv.HealthyReference.fit(mat)
    z = ref.zscore(mat.iloc[0])
    assert abs(np.median(z)) < 1.0


def test_arm_summary_flags_alteration():
    idx = [f"chr1:{i}" for i in range(20)]
    z = pd.Series(np.concatenate([np.full(10, 4.0), np.zeros(10)]), index=idx)
    arms = pd.Series(["chr1p"] * 10 + ["chr1q"] * 10, index=idx)
    out = cnv.arm_summary(z, arms)
    assert out["cnv_n_arms_altered"] >= 1


def test_technical_covariates_present(healthy_file):
    frags = read_fragments(healthy_file)
    binned = fx.bin_short_long(frags)
    tech = cnv.technical_covariates(frags, binned)
    assert tech["tech_n_fragments"] == len(frags)
    assert 0 <= tech["tech_frac_bins_zero"] <= 1


# ------------------------------------------------------------------ features

def test_extract_sample_features_blocks(healthy_file):
    blocks = extract_sample_features(healthy_file, "S1")
    for name in ("size_summary", "size_hist", "delfi", "technical",
                 "_coverage", "_window_arms"):
        assert name in blocks
    assert blocks["size_hist"].size == 351


def test_stack_blocks_shape(healthy_file, cancer_file):
    per = {"h": extract_sample_features(healthy_file, "h"),
           "c": extract_sample_features(cancer_file, "c")}
    stacked = stack_blocks(per)
    assert stacked["size_hist"].shape[0] == 2
    assert set(stacked["size_hist"].index) == {"h", "c"}


# ------------------------------------------------------------------ modeling

def test_nested_cv_on_separable_data():
    rng = np.random.default_rng(0)
    y = np.array([0] * 20 + [1] * 20)
    X = pd.DataFrame(rng.normal(0, 1, (40, 8)))
    X.iloc[20:, 0] += 3.0
    res = modeling.nested_cv(X, y, n_repeats=2)
    assert res.auc_mean > 0.9


def test_nested_cv_on_pure_noise_is_near_chance():
    rng = np.random.default_rng(1)
    y = np.array([0] * 20 + [1] * 20)
    X = pd.DataFrame(rng.normal(0, 1, (40, 50)))
    res = modeling.nested_cv(X, y, n_repeats=3)
    # the whole point of nested CV: high-dimensional noise must not look good
    assert res.auc_mean < 0.75


def test_permutation_null_centres_on_chance():
    rng = np.random.default_rng(2)
    y = np.array([0] * 15 + [1] * 15)
    X = pd.DataFrame(rng.normal(0, 1, (30, 10)))
    null, p = modeling.permutation_null(X, y, n_perm=25, n_repeats=1)
    assert 0.3 < null.mean() < 0.7
    assert 0 < p <= 1


def test_assemble_matrix_drops_constant_features():
    blocks = {"a": pd.DataFrame({"x": [1, 2, 3], "const": [5, 5, 5]},
                                index=["s1", "s2", "s3"])}
    X = modeling.assemble_matrix(blocks, ["a"], ["s1", "s2", "s3"])
    assert "const" not in X.columns and "x" in X.columns


# ------------------------------------------------------------------ dilution

def test_mix_fragments_respects_fraction(healthy_file, cancer_file):
    rng = np.random.default_rng(0)
    h = read_fragments(healthy_file)
    c = read_fragments(cancer_file)
    mixed = mix_fragments(c, h, 0.5, 10_000, rng)
    assert len(mixed) == 10_000


def test_mixture_size_is_monotone(healthy_file, cancer_file):
    rng = np.random.default_rng(0)
    h = read_fragments(healthy_file)
    c = read_fragments(cancer_file)
    medians = [np.median(mix_fragments(c, h, f, 20_000, rng)["length"])
               for f in (0.0, 0.5, 1.0)]
    assert medians[0] >= medians[1] >= medians[2]


def test_estimate_lod_picks_lowest_passing_level():
    df = pd.DataFrame({
        "fraction": [0.01] * 4 + [0.1] * 4 + [0.5] * 4,
        "score": [0.1, 0.1, 0.1, 0.9] + [0.9] * 4 + [0.9] * 4,
    })
    assert estimate_lod(df, threshold=0.5) == 0.1
