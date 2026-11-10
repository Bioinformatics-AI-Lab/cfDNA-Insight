# Limitations

Stated here rather than at the end of a results section, because they determine
what this project can and cannot be claimed to show.

## Limitations specific to the result obtained

**Stage dependence.** The cohort is heavily early-stage (9 Stage IA, 9 Stage
IIA of 25 breast). Per-stage breakdown at this n would not be interpretable,
so none is reported, and whether detection degrades at Stage I specifically is
not answerable from these data.

**The fragment-size null.** DELFI gives 0.625 and fragment-size summaries
0.500. Three explanations remain open: coverage too low for stable per-window
ratios at 1-3x, disease too early to shed detectably, or the within-sample
z-scoring in `delfi_profile` discarding the global size shift. The last is a
choice made here (see methodology) and could be tested by reporting the
uncentred profile, at the cost of reintroducing the pre-analytical
confounding the centring removes.

**End motifs are not a separate modality.** 0.894 raw, 0.593 with chromosome
entropy partialled out. Any downstream use should treat the block as a second
view of copy number.

**Mediator-adjustment columns.** Partialling out aneuploidy inflates unrelated
blocks including the negative control. Only decreases in those columns carry
information.

**Two features on 66 samples.** `an_chr_entropy` and `an_coverage_cv` at 0.866
is consistent with a real effect and also with a lucky pair of summary
statistics. The paired test (30 of 33, p = 3e-07) is the stronger evidence.
Replication in an independent cohort would settle it; none was available.

**The LOD is donor-specific and the panel is small.** Across the eight liver
pairs the own-baseline LOD ranges 0.067 to 0.289, a 4.3x spread that tracks
each donor's case-control difference rather than stage: the best and worst
pairs are both stage 1a. The median, 0.108, comes from eight donors of one
cancer type; the breast pairs were not run. Any single-donor figure is close
to meaningless, which an earlier version of the README learned the hard way by
quoting 0.08 from a better-than-typical pair.

**Population-threshold detection essentially fails.** Four of eight pairs were
not flagged at any tested fraction, undiluted included. That is the honest
read of this feature as a tumour-naive screening marker on this cohort.

**Leave-one-cancer-type-out is not external validation.** 0.824 overall,
0.969 breast to liver, from one lab and five shared batches. It tests transfer
across tumour type only. The breast-held-out direction (0.773) trains on 16
samples and should not be quoted alone.

## Design

**Retrospective, case-control, small n.** Cohorts assembled from published
studies compare diagnosed patients against healthy volunteers. A screening assay
has to work on an undiagnosed, unselected population where prevalence is well
under 1% and the differential diagnosis includes benign disease. Case-control
AUC systematically overstates screening performance, sometimes by a lot.

**Cases and controls may not be matched.** Age, sex, collection site, collection
protocol and processing batch are recorded where the source study published them
(`data/metadata/`) but cannot be balanced retrospectively. Any of them can
produce separation that has nothing to do with disease.

**Batch and label are partially confounded by construction.** In most source
studies, cancer and healthy samples were collected through different pathways.
No amount of downstream modelling fixes this; the technical-covariate control
measures how bad it is, but it cannot remove it.

## Data

**Fragment files, not raw reads.** Alignment parameters are fixed by FinaleDB's
pipeline. Their uniform reprocessing is a genuine advantage for cross-study
comparison, but it means this project cannot investigate whether the choice of
aligner or trimming affects the result ,  the very thing it would want to check
when a cross-study comparison fails.

**Reads trimmed to 50 bp.** Fragment length inference is preserved (trimming
removes inner bases, not the outer 5' ends), but anything requiring read
sequence beyond the ends is unavailable.

**Low coverage.** Most public cfDNA WGS is 1-3×. That is sufficient for
megabase-scale fragmentation profiles and arm-level copy number and insufficient
for anything at finer resolution.

**No matched buffy coat.** Clonal haematopoiesis cannot be excluded. This is one
reason no mutation-based feature exists here.

## Methods

**GC correction is approximate.** A quantile-knot interpolation rather than
LOESS. The difference is small at these coverages but it is a difference.

**The CNV module is not a copy-number caller.** No segmentation, no HMM, no
purity/ploidy solution, no tumour-fraction estimate. Arm-level z-scores against a
healthy panel, nothing more.

**The healthy reference panel is fit outside CV folds by default.** For arm-level
medians across dozens of samples the leak is small, but it is not zero. Both the
default and the fold-safe variant are reported so the size of the difference is
visible rather than assumed.

**Chromosomes X and Y are excluded.** This removes sex as a trivial confounder
but also discards real signal, particularly for sex-specific cancers.

**The reported AUC interval is a stability interval across CV repeats**, not a
confidence interval for the population AUC. At n ≈ 50 the true interval is
considerably wider. Both are given.

## The dilution series

**The x-axis is plasma fraction, not tumour fraction.** The reported LOD is an
upper bound on the true ctDNA LOD, by a factor equal to the donor's own tumour
fraction.

**Mixing is with replacement**, so fragment diversity at high simulated depth is
lower than a real library's would be.

**Duplication rate, error profile and library complexity do not vary with
dilution.** In a real dilution experiment all three would.

**One donor pair per series.** Estimates from a single cancer donor and a single
healthy donor do not generalise; running several pairs would give a range.

## What would fix each of these

| Limitation | Fix |
|---|---|
| Case-control design | Prospective screening cohort with clinical follow-up |
| Batch confounding | Cases and controls collected and processed in the same batches |
| Small n | Hundreds to thousands of samples per class |
| No tumour fraction | ichorCNA, or a tumour-informed panel with matched tissue |
| CHIP confounding | Matched buffy coat sequencing |
| Simulated LOD | Real contrived dilution series with orthogonally measured tumour fraction |
| Single-study evaluation | Three or more independent cohorts, held out cleanly |
| Aneuploidy vs fragmentation unresolved | Higher-coverage WGS on the same donors |
| Two-feature model unreplicated | An independent matched cohort |
| LOD from 8 liver donors only | Run the 25 breast pairs; add a second cohort |

None of these are available for a public-data portfolio project. That is the
reason the evaluation is built the way it is: if the design cannot be fixed, the
next best thing is to measure precisely how much the design is doing.
