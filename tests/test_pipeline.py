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
from cfdna_insight.io import SampleSheet, _sniff_columns, read_fragments


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


# --------------------------------------------------- genome build selection

def test_hg19_and_hg38_centromeres_differ():
    """Regression: make_bins once used hg38 centromeres for every build."""
    b38 = bins.make_bins(bin_size=1_000_000, genome="hg38")
    b19 = bins.make_bins(bin_size=1_000_000, genome="hg19")
    arms38 = b38.set_index("bin_id")["arm"]
    arms19 = b19.set_index("bin_id")["arm"]
    common = arms38.index.intersection(arms19.index)
    assert (arms38[common] != arms19[common]).any(), \
        "hg19 and hg38 arm assignments are identical - build is being ignored"


def test_hg19_arm_split_matches_hg19_centromere():
    b = bins.make_bins(bin_size=1_000_000, genome="hg19")
    for chrom in ("chr1", "chr12", "chr22"):
        cen_start = bins.HG19_CENTROMERE[chrom][0]
        p_arm = b[(b["chrom"] == chrom) & (b["arm"] == f"{chrom}p")]
        assert p_arm["start"].max() < cen_start


def test_hg19_sizes_used_for_hg19():
    b19 = bins.make_bins(bin_size=1_000_000, genome="hg19")
    b38 = bins.make_bins(bin_size=1_000_000, genome="hg38")
    assert (b19["chrom"] == "chr1").sum() == bins.HG19_SIZES["chr1"] // 1_000_000
    assert (b38["chrom"] == "chr1").sum() == bins.HG38_SIZES["chr1"] // 1_000_000


def test_unsupported_genome_rejected():
    with pytest.raises(ValueError, match="unsupported genome"):
        bins.make_bins(genome="hg17")


# ------------------------------------------- real-world fragment file layouts

def _write_crag_style(path):
    """Replica of the CRAG layout: '#' header, 8 columns, no 'chr' prefix."""
    import gzip
    rng = np.random.default_rng(7)
    with gzip.open(path, "wt") as fh:
        fh.write("#chrom\tstart\tend\tname\tmapq\tstrand\tcigar1\tcigar2\n")
        for _ in range(3000):
            c = int(rng.integers(1, 23))
            s = int(rng.integers(10000, 40_000_000))
            ln = int(rng.integers(80, 250))
            q = int(rng.integers(0, 61))
            st = "+" if rng.random() < 0.5 else "-"
            fh.write(f"{c}\t{s}\t{s+ln}\t.\t{q}\t{st}\t100M\t100M\n")
    return str(path)


def test_reads_crag_eight_column_format(tmp_path):
    p = _write_crag_style(tmp_path / "MAL9999.hg19.frag.bed.gz")
    df = read_fragments(p, min_mapq=30)
    assert len(df) > 0
    assert list(df.columns) == ["chrom", "start", "end", "mapq", "strand", "length"]
    assert df["chrom"].astype(str).str.startswith("chr").all()
    assert (df["mapq"] >= 30).all()
    assert (df["length"] == df["end"] - df["start"]).all()


def test_sniff_names_every_column(tmp_path):
    p = _write_crag_style(tmp_path / "s.frag.bed.gz")
    names = _sniff_columns(p)
    assert len(names) == 8
    assert names[:6] == ["chrom", "start", "end", "name", "mapq", "strand"]


def test_cigar_columns_not_loaded(tmp_path):
    """The big string columns must stay out of memory."""
    p = _write_crag_style(tmp_path / "s.frag.bed.gz")
    df = read_fragments(p)
    assert not any(c.startswith("extra") for c in df.columns)


def test_headerless_bed_still_works(tmp_path):
    import gzip
    p = tmp_path / "plain.frag.gz"
    with gzip.open(p, "wt") as fh:
        for i in range(500):
            fh.write(f"chr1\t{10000+i*300}\t{10000+i*300+160}\t40\t+\n")
    df = read_fragments(str(p))
    assert len(df) == 500
    assert set(["chrom", "start", "end", "mapq", "strand"]) <= set(df.columns)


# ------------------------------------------------- fold-safe CNV transformer

def _toy_coverage(n_healthy=12, n_cancer=12, n_win=60, seed=0):
    rng = np.random.default_rng(seed)
    cov = rng.normal(1.0, 0.03, (n_healthy + n_cancer, n_win))
    cov[n_healthy:, :15] *= 1.15          # arm-level gain in the cancer group
    idx = [f"s{i}" for i in range(len(cov))]
    cols = [f"chr1:{i}" for i in range(n_win)]
    arms = pd.Series(["chr1p"] * 30 + ["chr1q"] * 30, index=cols)
    y = np.array([0] * n_healthy + [1] * n_cancer)
    return pd.DataFrame(cov, index=idx, columns=cols), arms, y


def test_reference_transformer_uses_training_healthy_only():
    cov, arms, y = _toy_coverage()
    tr = cnv.HealthyReferenceTransformer(arms)
    train = np.arange(0, 18)
    tr.fit(cov.iloc[train], y[train])
    # the reference must be the median of TRAINING healthy rows only
    expected = cov.iloc[[i for i in train if y[i] == 0]].median(axis=0)
    assert np.allclose(tr.reference_.median.to_numpy(), expected.to_numpy())


def test_reference_transformer_needs_y():
    cov, arms, y = _toy_coverage()
    with pytest.raises(ValueError, match="needs y"):
        cnv.HealthyReferenceTransformer(arms).fit(cov)


def test_reference_transformer_rejects_too_few_healthy():
    cov, arms, y = _toy_coverage()
    with pytest.raises(ValueError, match="need >=3 healthy"):
        cnv.HealthyReferenceTransformer(arms).fit(cov.iloc[-4:], y[-4:])


def test_transform_output_is_stable_shape():
    cov, arms, y = _toy_coverage()
    tr = cnv.HealthyReferenceTransformer(arms).fit(cov, y)
    out = tr.transform(cov)
    assert len(out) == len(cov)
    assert list(out.columns) == tr.feature_names_


def test_foldsafe_cnv_runs_and_detects_real_signal():
    cov, arms, y = _toy_coverage()
    res = modeling.nested_cv_foldsafe_cnv(cov, arms, y, n_repeats=2)
    assert 0.0 <= res.auc_mean <= 1.0
    assert res.auc_mean > 0.7      # the planted gain is real, not a leak


def test_foldsafe_cnv_near_chance_on_label_free_coverage():
    """No planted signal -> fold-safe CNV must not manufacture one."""
    rng = np.random.default_rng(3)
    cols = [f"chr1:{i}" for i in range(60)]
    cov = pd.DataFrame(rng.normal(1.0, 0.03, (24, 60)),
                       index=[f"s{i}" for i in range(24)], columns=cols)
    arms = pd.Series(["chr1p"] * 30 + ["chr1q"] * 30, index=cols)
    y = np.array([0] * 12 + [1] * 12)
    res = modeling.nested_cv_foldsafe_cnv(cov, arms, y, n_repeats=2)
    assert res.auc_mean < 0.75


def test_null_pvalue_shape():
    rng = np.random.default_rng(4)
    y = np.array([0] * 10 + [1] * 10)
    X = pd.DataFrame(rng.normal(0, 1, (20, 5)))
    p, mu, sd = modeling.null_pvalue(X, y, observed=0.75, n_perm=20)
    assert 0 < p <= 1 and 0.3 < mu < 0.7 and sd > 0


def test_transformer_feature_names_span_training_fold():
    """Names come from the whole training fold, not one probe row."""
    cov, arms, y = _toy_coverage()
    tr = cnv.HealthyReferenceTransformer(arms).fit(cov, y)
    out = tr.transform(cov)
    assert list(out.columns) == tr.feature_names_
    assert out.notna().any(axis=0).all()


def test_assemble_matrix_drops_sparse_columns():
    """An arm seen in only a few samples must not survive into CV."""
    X = pd.DataFrame({
        "dense": np.arange(10.0),
        "sparse": [1.0, 2.0] + [np.nan] * 8,
        "const": np.ones(10),
    }, index=[f"s{i}" for i in range(10)])
    out = modeling.assemble_matrix({"a": X}, ["a"], list(X.index))
    assert list(out.columns) == ["dense"]


def test_assemble_matrix_keeps_mostly_present_columns():
    X = pd.DataFrame({
        "mostly": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, np.nan],
        "dense": np.arange(10.0),
    }, index=[f"s{i}" for i in range(10)])
    out = modeling.assemble_matrix({"a": X}, ["a"], list(X.index))
    assert set(out.columns) == {"mostly", "dense"}


def test_assemble_matrix_errors_when_nothing_usable():
    X = pd.DataFrame({"const": np.ones(6)}, index=[f"s{i}" for i in range(6)])
    with pytest.raises(ValueError, match="no usable features"):
        modeling.assemble_matrix({"a": X}, ["a"], list(X.index))


# ------------------------------------------------------ confounder residuals

def _load_audit():
    import importlib.util, os
    p = os.path.join(os.path.dirname(__file__), "..", "scripts",
                     "04_confounder_audit.py")
    spec = importlib.util.spec_from_file_location("audit", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_residualise_removes_linear_dependence():
    rng = np.random.default_rng(0)
    depth = rng.normal(7, 0.3, 60)
    X = pd.DataFrame({
        "driven": 3.0 * depth + rng.normal(0, 0.05, 60),
        "free": rng.normal(0, 1, 60),
    })
    out = _load_audit().residualise(X, depth)
    assert abs(np.corrcoef(out["driven"], depth)[0, 1]) < 0.05
    assert abs(np.corrcoef(out["free"], depth)[0, 1]) < 0.35


def test_residualise_preserves_independent_signal():
    rng = np.random.default_rng(1)
    depth = rng.normal(7, 0.3, 60)
    signal = np.r_[np.zeros(30), np.ones(30)]
    X = pd.DataFrame({"f": signal + rng.normal(0, 0.1, 60)})
    out = _load_audit().residualise(X, depth)
    assert abs(np.corrcoef(out["f"], signal)[0, 1]) > 0.8


def test_residualise_tolerates_constant_and_nan():
    rng = np.random.default_rng(2)
    depth = rng.normal(7, 0.3, 40)
    X = pd.DataFrame({
        "const": np.ones(40),
        "sparse": [np.nan] * 36 + [1.0, 2.0, 3.0, 4.0],
        "ok": rng.normal(0, 1, 40),
    })
    out = _load_audit().residualise(X, depth)
    assert out.shape == X.shape
    assert (out["const"] == 1.0).all()


# ---------------------------------------- technical / aneuploidy separation

def test_technical_block_holds_no_biological_features(healthy_file):
    """Regression: chromosome entropy once sat in the negative control."""
    frags = read_fragments(healthy_file)
    binned = fx.bin_short_long(frags)
    tech = cnv.technical_covariates(frags, binned)
    for banned in ("entropy", "coverage_cv"):
        assert not any(banned in k for k in tech.index), \
            f"{banned} is a copy-number readout and must not be in the " \
            f"negative control"
    assert set(tech.index) == {
        "tech_n_fragments", "tech_log10_n_fragments",
        "tech_median_mapq", "tech_frac_bins_zero"}


def test_aneuploidy_block_present(healthy_file):
    frags = read_fragments(healthy_file)
    binned = fx.bin_short_long(frags)
    an = cnv.aneuploidy_covariates(frags, binned)
    assert set(an.index) == {"an_chr_entropy", "an_coverage_cv"}
    assert np.isfinite(an["an_chr_entropy"])


def test_split_legacy_technical_moves_columns():
    from cfdna_insight.features import split_legacy_technical
    idx = ["s1", "s2"]
    blocks = {"technical": pd.DataFrame(
        {"tech_n_fragments": [1.0, 2.0], "tech_chr_entropy": [0.9, 0.8],
         "tech_coverage_cv": [0.3, 0.4]}, index=idx)}
    out = split_legacy_technical(blocks)
    assert list(out["technical"].columns) == ["tech_n_fragments"]
    assert set(out["aneuploidy"].columns) == {"an_chr_entropy", "an_coverage_cv"}


def test_split_legacy_technical_is_noop_when_already_split():
    from cfdna_insight.features import split_legacy_technical
    blocks = {"technical": pd.DataFrame({"tech_n_fragments": [1.0]},
                                        index=["s1"])}
    out = split_legacy_technical(blocks)
    assert "aneuploidy" not in out
