# Methodology

## 1. Input and filtering

Fragment files (FinaleDB `.frag.gz` / `.frag.bed.gz`): BED-like, one row per
sequenced fragment, columns `chrom start end [name] mapq strand`, 0-based
half-open.

Filters, and why each one:

| Filter | Value | Reason |
|---|---|---|
| MAPQ | ≥ 30 | Multi-mapping fragments cluster in repeats, where fragment length is unreliable. |
| Length | 50–1000 bp | Below 50 is largely adapter artefact; above 1000 is chimeric or discordant pairs. |
| Chromosomes | autosomes | chrX/chrY create sex-driven separation if the case and control groups are not sex-balanced. chrM copy number varies enormously between individuals and adds variance without signal. |
| Centromeric windows | dropped | Mappability collapses; coverage there is an artefact of the aligner, not the sample. |

## 2. Binning

100 kb bins for GC correction, aggregated 50× into 5 Mb windows for reporting.
The two scales serve different purposes. GC bias operates at the scale over
which GC content actually varies, which is tens of kilobases; correcting at
5 Mb averages the bias away before it can be removed. But a short/long ratio in
a single 100 kb bin at 1–2× coverage is computed from a few hundred fragments
and is mostly Poisson noise, so features are only reported at 5 Mb.

Fragments are assigned by midpoint rather than start. Using the start
systematically shifts assignment upstream by half a fragment length — negligible
at 5 Mb, but the habit matters if the same code is ever reused at nucleosome
scale.

Incomplete windows at chromosome ends are dropped rather than emitted short. A
3 Mb window carries 60% of the counts of a 5 Mb one and would look like a
consistent deletion in every sample.

## 3. GC correction

Multiplicative, against a piecewise-linear fit of count on GC content across
quantile knots: `corrected = observed × median(counts) / predicted(GC)`.

This is a LOESS fit in the DELFI paper. The quantile-knot interpolation used
here behaves near-identically on binned count data and avoids a statsmodels
dependency; the difference is well inside the noise at these coverages. Bins
with fewer than 50% non-N bases get GC = NaN and are excluded from the fit,
since they would otherwise anchor the low-GC end of the curve with garbage.

GC bias is not a nuisance to be waved at here. It varies with library prep and
PCR cycle number, both of which differ between studies and often between case
and control batches within a study. Uncorrected GC bias is one of the most
common routes by which a "biological" cfDNA classifier is actually a batch
classifier.

## 4. Feature blocks

### Size distribution
Density-normalised 1 bp histogram over 50–400 bp, plus nine scalars: median,
mean, IQR, fraction below 150 and 120, fraction in the 167–181 mononucleosome
peak, dinucleosomal fraction, global short/long ratio, and a 10 bp periodicity
score.

The periodicity score is the magnitude of the DFT component at period 10 bp
after cubic detrending, normalised by total spectral power. Physically it
measures how cleanly the DNA was wrapped on histone octamers when it was cut —
the helical pitch means nuclease access is periodic. It's a biological feature,
but its main practical use is QC: a sample that has lost its 10 bp ripple has
been degraded, and its size distribution should not be trusted.

### DELFI short/long profile
Short (100–150 bp) to long (151–220 bp) count ratio per 5 Mb window,
GC-corrected at 100 kb, then **z-scored within sample** across windows.

The within-sample centring is the single most consequential choice in the
feature set. It throws away the mean short/long ratio — a real biological signal
— and keeps only the profile shape. That trade is worth making because the mean
also tracks time-to-plasma-separation, storage temperature, freeze-thaw cycles
and extraction kit, all of which routinely differ between a prospectively
collected cancer cohort and a convenience-sampled healthy cohort. The global
mean is still available separately in the size-summary block, so the ablation
can show what it contributes on its own.

### Copy number
Depth-normalised, GC-corrected 5 Mb coverage, converted to z-scores against a
healthy reference panel (per-window median and MAD), then collapsed to per-arm
medians. Arm medians are rescaled by √k, where k is the number of windows in the
arm, so that a |z| > 3 threshold means the same thing on chr21p as on chr1q.

Two scalars summarise the profile: mean |z| across arms (aneuploidy burden) and
the count of arms exceeding |z| > 3.

This is not ichorCNA. There is no HMM, no segmentation, no purity/ploidy
solution, no tumour fraction. It is a coverage-imbalance feature, and it is
labelled as one.

### End motifs
Frequency of each of 256 4-mers at fragment 5' ends, taken from the reference at
`[start, start+4)` for + strand fragments and the reverse complement of
`[end-4, end)` for − strand. Motif diversity score is the Shannon entropy of the
spectrum normalised by log(256), so 1 means a maximally diverse cut repertoire.

The biology: DNASE1L3 and DFFB cut with different sequence preferences, and
DNASE1L3 activity is reduced in several cancers. A narrowed motif spectrum is a
readout of which nuclease did the cutting.

### Technical covariates
Fragment count, log fragment count, median MAPQ, fraction of empty bins,
coverage coefficient of variation, chromosome-distribution entropy.

These are never used in the biological model. Their only job is the negative
control in §5.

## 5. Evaluation

### Nested cross-validation
Outer 5-fold stratified, inner 4-fold for hyperparameter selection, repeated 10
times with different partitions. Reporting the best inner-loop score as if it
were an out-of-sample estimate is the most common inflation mechanism in small-n
omics papers; nesting is what prevents it.

Ridge logistic regression is the default, not gradient boosting. With n in the
tens and p in the hundreds, a strongly regularised linear model is better
calibrated, more interpretable, and far less able to latch onto one confounded
window. GBM is available for comparison.

The interval reported around the AUC is a percentile interval **across CV
repeats**. It measures stability of the estimate under repartitioning, not
sampling variability of the population AUC. At n = 50 the true confidence
interval is considerably wider. Both are stated in the results so the
distinction is not quietly lost.

### Permutation null
Labels are shuffled and the **entire nested CV is re-run**, 200 times. Permuting
labels while reusing a model fit on the true labels tests nothing.

The resulting null is the reference distribution the observed AUC should be read
against. At n = 50 with 500 features, a null centred at 0.5 with a standard
deviation near 0.08 means an AUC of 0.65 is unremarkable.

### Feature-block ablation
Each block is evaluated alone and in combination. The specific question: does
the 500-dimensional DELFI profile beat the 9 scalar size summaries? If it
doesn't, the profile is decoration and the write-up says so.

### Technical-covariate negative control
A model trained on library depth, MAPQ and coverage dispersion alone. It should
sit at chance. If it reaches AUC > 0.65 the pipeline emits an explicit warning,
because at that point the biological model's AUC is partly measuring which batch
a sample came from and any headline number must carry that caveat.

This is the most informative single analysis in the repo and the one most
comparable published work omits.

### Leave-one-study-out
Train on all studies but one, test on the held out. This is the number that
matters for a screening assay: within-cohort CV mostly measures whether a model
can learn one lab's signature. Expect a substantial drop. Report it anyway.

With a single study in the sample sheet the pipeline refuses to produce a
cross-study number and states that generalisation is untested, rather than
letting a within-cohort AUC be read as if it generalised.

### Learning curve
AUC against training-set size. A curve still climbing steeply at n = all means
the estimate has not converged and the headline AUC should not be quoted as if
it had.

## 6. In silico dilution

Synthetic samples are built by drawing a fixed total number of fragments,
a fraction *f* from a cancer donor and 1−*f* from a healthy donor, then
recomputing features and scoring with a model trained on donors excluding both.
Excluding the donors is essential: a model that has seen the cancer donor will
recognise it at any dilution for the wrong reason.

The detection threshold is set on healthy samples at 95% specificity first, and
sensitivity is read off afterwards. Screening assays are specificity-anchored —
at population prevalence a 95%-specific test generates far more false than true
positives — so choosing a threshold to maximise sensitivity and reporting the
resulting specificity inverts the actual design constraint.

Limits of this simulation, all of which belong in the results text:

- *f* is a **plasma fraction, not a tumour fraction**. The cancer donor's plasma
  is itself mostly haematopoietic cfDNA, so mixing 10% of a sample carrying 5%
  ctDNA gives roughly 0.5% ctDNA. The reported LOD is an upper bound.
- Fragments are drawn with replacement, so at high depth the same fragment can
  appear twice. Real libraries at that depth have more molecular diversity.
- Duplication rate, error profile and library complexity do not change with
  dilution here. In a real dilution they all would.
- LOD is reported as "LOD95 ≤ x" because with five replicates per level the
  estimate is coarse and bounded by the grid of fractions tested.
