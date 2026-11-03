# cfDNA-Insight

**Can cancer-associated cfDNA be distinguished from non-cancer cfDNA using genomic and fragmentomic features — and how much of the measured performance survives proper validation?**

Short answer, on 33 matched case–control pairs: **large-scale chromosomal imbalance separates them; fragment length does not.** The best-known feature in cfDNA fragmentomics returns a null here, a two-feature aneuploidy proxy beats a 1,193-feature model, and the apparent end-motif signal turns out to be the aneuploidy signal measured a second way.

The second half of the title is the point of the project. Reproducing a fragmentomics classifier is a weekend's work; the reason cfDNA classifiers fail in the clinic is almost never that nobody could fit one. So this repo reports the classifier *and* the analyses that try to break it — including the two occasions when those analyses caught errors in the work itself.

---

## Results

33 cancer patients (25 breast, 8 liver), each with a control matched on age, sex, smoking and alcohol use and processed in the same library-prep and sequencing batch. hg19 fragment coordinates, 11.7–45.2M fragments per sample, ~1–3× WGS.

Nested cross-validation (5 outer × 4 inner, 10 repeats, ridge logistic regression):

| Feature block | Features | AUC | 95% band |
|---|---:|---:|---|
| End motifs (4-mer at fragment 5′ ends) | 256 | 0.899 | 0.881–0.919 |
| **Aneuploidy proxies** | **2** | **0.866** | 0.822–0.898 |
| All biological | 1193 | 0.849 | 0.814–0.895 |
| CNV arm z-scores — precomputed reference | 44 | 0.835 | 0.786–0.864 |
| CNV arm z-scores — fold-safe reference | 44 | 0.680 | — |
| DELFI short/long profile | 530 | 0.625 | 0.524–0.684 |
| **Technical covariates (negative control)** | **3** | **0.533** | 0.474–0.596 |
| Fragment size summary | 10 | 0.500 | 0.436–0.569 |
| Fragment size histogram | 351 | 0.462 | 0.360–0.523 |

Full model: AUC 0.849, permutation p = 0.005 (200 permutations — the floor).
Sensitivity 0.70 at 95% specificity.
Leave-one-cancer-type-out: 0.824 overall; trained on breast, tested on liver, 0.969.

### Fragment size does not separate these pairs

Size summary 0.500, size histogram 0.462, DELFI short/long profile 0.625 — and DELFI is unchanged by depth adjustment (0.629), so it isn't a coverage artefact.

This is the substantive negative result. The DELFI feature (Cristiano et al. 2019) is the most cited construct in cfDNA fragmentomics, and this cohort is an unusually fair test of it: age, sex, smoking, alcohol and batch are all held by matching rather than adjusted away. It still returns nothing.

Three explanations are consistent with the data and these data cannot separate them: coverage too low for stable per-window short/long ratios; disease too early to shed detectably (9 Stage IA, 9 Stage IIA); or the within-sample z-scoring discarding the global size shift that carries the real signal. The third is a design choice this repo made deliberately and documents in `docs/methodology.md`.

### The signal that is there is simple

`an_chr_entropy` — the normalised Shannon entropy of per-chromosome fragment proportions — separates cases from **their own matched controls in 30 of 33 pairs** (Wilcoxon p = 3×10⁻⁷). With `an_coverage_cv` it reaches AUC 0.866 on two features, beating the full 1,193-feature model.

Being a paired result is what makes it strong. Each control was matched on age, sex, smoking and alcohol and run in the same batch, so those explanations are excluded by construction rather than by argument. And neither feature needs a reference panel, so neither can leak information across CV folds.

### End motifs are the same finding, not a second one

End motifs score highest in the table (0.899) and that number is misleading. Residualising on chromosome entropy drops them to 0.593.

The mechanism: 4-mers are read off the reference genome at fragment termini. A sample whose coverage has been redistributed by arm-level gains and losses is sampling a different mix of sequence context, and its 4-mer spectrum shifts with no change in nuclease activity at all. Motif and aneuploidy are one mechanism observed through two instruments, and the write-up says so rather than counting them twice.

---

## What the stress tests caught

These stand independent of whether the classifier works, and two of them are corrections to this repo's own earlier output.

**The negative control was contaminated, and that is how we found out.** A model trained only on technical covariates initially reached AUC 0.886 — higher than the biological model — which would have invalidated everything. Investigation showed two of its six covariates were copy-number readouts: chromosome entropy and bin-coverage dispersion are exactly what aneuploidy perturbs. Moved to their own block, the control fell to **0.533 (p = 0.37)**, where a negative control belongs. A test in `tests/test_pipeline.py` now fails if anything matching `entropy` or `coverage_cv` reappears in the technical block.

**A standard shortcut inflates AUC by 0.155.** Fitting the healthy coverage reference across all healthy samples — rather than rebuilding it inside each training fold — gives 0.835 instead of 0.680 at n=66, and 0.963 instead of 0.786 at n=16. Each test-fold healthy sample contributes roughly 1/n of the reference it is later scored against, pulling it toward z = 0. Both numbers are reported; `--cnv-fold-safe` computes them side by side.

**Partialling out aneuploidy manufactures signal.** Removing chromosome entropy raised fragment-size AUC from 0.489 to 0.713 and size-histogram from 0.442 to 0.765 — which looks like hidden signal emerging. It isn't. The negative control rose in the same column (0.507 → 0.606), and nothing about removing aneuploidy can reveal genuine library signal. Aneuploidy is a *consequence* of the label, so conditioning on it is conditioning on a mediator, which induces association. In that column, falls are interpretable and rises are not. `04_confounder_audit.py` now flags this automatically using the negative control as the diagnostic.

**The case/control pairing was verified, not assumed.** Table S2 states each `MAL####m` is the matched control for `MAL####`. `00_check_pairing.py` confirms it from the data: exact fragment-coordinate overlap between paired files is indistinguishable from the unrelated-sample baseline (median ratio 1.0 across 8 liver pairs), so they are different people rather than two versions of one library. Coordinates only — fragment length never enters, so the check cannot contaminate the hypothesis under test.

**Batch is orthogonal to label by construction.** Table S2 records 5 library-prep batches mapping 1:1 to 5 sequencing batches, with every pair's case and control in the same one. This is why the negative control is a real test here rather than a formality: in the usual public cfDNA cohort — patients collected prospectively in clinic, controls pulled from a biobank — it fires by construction and tells you nothing.

---

## Scope limits, stated up front

- **Starts from fragment coordinates, not FASTQ.** No alignment is performed. This is a deliberate choice with a real cost: alignment parameters are fixed upstream and cannot be varied here.
- **No somatic variant calling.** At 1–3× coverage, per-base variant calling in plasma is not meaningful, and clonal haematopoiesis would confound it even at depth. There is no `variants/` module and no mutation-based feature.
- **No tumour-fraction estimation.** The CNV module produces GC-corrected, healthy-referenced arm-level z-scores. It is not ichorCNA and does not solve for purity or ploidy.
- **Retrospective case–control, n=66, one cohort.** Cross-cohort generalisation is untested and the pipeline refuses to report a number for it. Leave-one-cancer-type-out is not a substitute: same lab, same protocol, same five batches.
- **The learning curve has not converged** (0.861 at n=52, 0.835 at n=66), so the point estimate should be quoted with that noted.
- **The reported band is stability across CV repartitions**, not a sampling confidence interval. At n=66 the true interval is wider.
- **Nothing here speaks to screening performance.** At population prevalence, the 95%-specificity operating point (sensitivity 0.70) would return far more false positives than true ones.

`docs/limitations.md` has the full list and what would fix each.

---

## Data

[Zenodo 6914806](https://zenodo.org/records/6914806), CC-BY-4.0 — the independent validation cohort from Zhou, Zheng, Fu et al., *Genome Medicine* 2022 (CRAG). Openly downloadable because the sequence is stripped and only fragment coordinates remain, which is all fragmentation analysis needs.

See `data/README.md` for the download procedure, the Table S2 label mapping, and why FinaleDB — the natural source for a genuine second cohort — was unavailable (its API returned 500s throughout).

## Quick start

```bash
conda env create -f environment.yml && conda activate cfdna-insight
pip install -e ".[dev,full]"

python scripts/00_smoke_test.py                      # ~2 min, no downloads
bash   scripts/00_fetch_data.sh                      # hg19 + ENCODE blacklist
bash   scripts/00_fetch_crag.sh all                  # 14.3 GB
python scripts/00_extract_metadata.py --xlsx <Additional file 1>.xlsx
python scripts/00_build_samplesheet.py
python scripts/00_check_pairing.py --raw data/raw/crag/liver

python scripts/01_extract_features.py --samples config/samples.tsv \
    --outdir results/features --fasta data/ref/hg19.fa --genome hg19
python scripts/02_train_evaluate.py --features results/features \
    --outdir results/full --cnv-fold-safe
python scripts/04_confounder_audit.py --features results/features
```

`00_smoke_test.py` runs the whole pipeline on synthetic fragment files first, so a broken environment fails in two minutes rather than three hours into a real run. Its numbers are meaningless — the simulator plants the effect it then finds.

Feature extraction is ~13 s/sample; all 66 run serially in about 15 minutes.

## Repository layout

```
src/cfdna_insight/     io, bins, fragmentomics, cnv, features, modeling,
                       dilution, plots, simulate
scripts/               00 smoke/fetch/metadata/samplesheet/pairing,
                       01 features, 02 evaluate, 03 dilution,
                       04 confounder audit
config/                config.yaml, samples.tsv, samples_by_cancer_type.tsv
tests/                 51 tests, all on simulated data — no downloads
docs/                  methodology, biological background, limitations
results/               summary.json, confounder_audit.json, figures, tables
```

Reproducibility note: every AUC in this README was reproduced to three decimal places independently on macOS/arm64 and Linux/x86-64.

## License

MIT for the code. The cfDNA data is CC-BY-4.0 and must be cited as Zhou, Zheng & Liu (2022) alongside the CRAG paper. `data/raw/` is gitignored.

## References

Zhou, Zheng, Fu et al., *Genome Medicine* 2022 (CRAG) · Cristiano et al., *Nature* 2019 (DELFI) · Snyder et al., *Cell* 2016 · Jiang et al., *Cancer Discovery* 2020 (end motifs) · Adalsteinsson et al., *Nature Communications* 2017 (ichorCNA)
