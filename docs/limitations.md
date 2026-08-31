# Limitations

Stated here rather than at the end of a results section, because they determine
what this project can and cannot be claimed to show.

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
aligner or trimming affects the result — the very thing it would want to check
when a cross-study comparison fails.

**Reads trimmed to 50 bp.** Fragment length inference is preserved (trimming
removes inner bases, not the outer 5' ends), but anything requiring read
sequence beyond the ends is unavailable.

**Low coverage.** Most public cfDNA WGS is 1–3×. That is sufficient for
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

None of these are available for a public-data portfolio project. That is the
reason the evaluation is built the way it is: if the design cannot be fixed, the
next best thing is to measure precisely how much the design is doing.
