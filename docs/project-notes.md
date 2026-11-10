# Project notes

Written for my own reference: why this project exists, how it was planned,
what actually happened including the wrong turns, and what the numbers mean.

## Why this project

I was moving toward liquid-biopsy and cfDNA roles and had no cfDNA work to
point at. My production experience is bulk and single-cell RNA-seq, WGS/WES
and methylation arrays at cohort scale, plus a lot of real-world clinical data
harmonisation. All adjacent, none of it cfDNA.

The obvious move was to reproduce a published fragmentomics classifier. I
decided against leading with that, because reproducing DELFI is a weekend of
work and would demonstrate almost nothing. Plenty of people can fit a
classifier. The reason cfDNA classifiers fail when they reach the clinic is
almost never that nobody could fit one; it is that the AUC was measuring
cohort assembly, library prep or batch.

So the project became: build the classifier, then spend most of the effort
trying to break it, and report whatever survives. That reframing is what made
it worth doing, and it turned out to be the right call, because the
stress tests caught two real errors in my own work.

## Getting the data

Almost all cfDNA WGS sits in dbGaP or EGA under controlled access, because raw
reads carry recoverable genotype. Months of application, usually an
institutional DTA.

The way around it: fragmentation analysis does not need sequence, only where
each fragment started and ended. Several groups publish fragment coordinates
with the sequence stripped, which makes the data openly downloadable.

FinaleDB is the obvious source and it was down. Its React frontend loaded but
both live API routes returned 500, which is why the page rendered blank. Not
fixable from outside.

The fallback turned out to be better. Zenodo record 6914806 holds the
independent validation cohort from the CRAG paper (Zhou et al., *Genome
Medicine* 2022) under CC-BY-4.0. Before downloading 14 GB I read the tar
headers over HTTP range requests to get the full file manifest without
transferring anything.

**The property that made this cohort worth using**, found by reading
Supplementary Table 2: 33 cancer patients, each with a control matched on age,
sex, smoking and alcohol, and *every pair processed in the same library-prep
and sequencing batch*. Five batches, case and control together in all 33.

That matters more than the sample count. The usual public cfDNA cohort has
patients collected prospectively in a clinic and controls pulled from a
biobank, which makes batch a proxy for label and quietly guarantees that any
negative control fires. Here it does not, which is what makes the negative
control an actual test rather than a formality.

### The sample-ID problem

Every sample appeared twice, `MAL0638` and `MAL0638m`, with nothing documented
about what the suffix meant. Two readings fit the file counts equally well:
case plus matched control (usable), or one sample sequenced twice (33 cancer
samples, no controls, unusable).

I wrote a script to settle it from the data before finding the answer on
paper: compare *exact fragment coordinates* between the paired files against
an unrelated-sample baseline. Two runs of one library share molecule
boundaries; two people share them at chance. Median ratio came out 1.0 across
8 pairs, so different people. Coordinates only, never fragment length, so the
check could not contaminate the hypothesis under test.

Table S2 then confirmed it directly: `Subject_ID` is the case,
`Control Subject_ID` is the matched control, and the control ID is always the
case ID plus `m`, for all 33 pairs with no exceptions.

Lesson worth keeping: I nearly used the wrong file. The Zenodo record contains
`Supplemantary_Table 2_Hotspots.zip`, which is per-cohort hotspot intervals
and has no sample IDs. The table I needed was Additional file 1 on the journal
page. Almost identical names, completely different content.

## What was built

Four feature blocks, each capturing a different physical property:

- **Fragment size.** Normal cfDNA peaks near 167 bp (147 bp wrapped on the
  nucleosome plus linker) with a 10 bp ripple from the helical pitch of DNA on
  the histone octamer. Tumour cfDNA is reported to be shorter.
- **DELFI short/long profile.** Ratio of 100-150 bp to 151-220 bp fragments in
  5 Mb windows, GC-corrected at 100 kb, z-scored within sample.
- **End motifs.** The 4-mer at each fragment's 5' end, which reflects which
  nuclease cut it.
- **Copy number.** Arm-level coverage z-scores against a healthy panel, plus
  two reference-free aneuploidy scalars.

Plus a **technical block** that is never used for biology: fragment count,
median MAPQ, fraction of empty bins. Its only job is to be the negative
control.

Evaluation: nested repeated CV, permutation null, feature-block ablation,
leave-one-cancer-type-out, learning curve, and an in silico dilution series.

## What happened, in order

**The negative control fired at 0.886**, higher than the 1,193-feature
biological model. If true, the whole result was an artefact.

It was my error. Two of the six technical covariates were chromosome entropy
and bin-coverage dispersion, and both are copy-number readouts: aneuploidy is
exactly what perturbs per-chromosome fragment proportions and coverage
variance. I had put a cancer signal inside the negative control, so it fired
by construction. Moved to their own block, the control fell to **0.533,
p = 0.37**. There is now a test that fails if anything matching `entropy` or
`coverage_cv` reappears in the technical block.

**The CNV reference panel leaked 0.155 AUC.** Fitting the healthy coverage
reference across all healthy samples, including ones that land in test folds,
gives 0.835; rebuilding it inside each training fold gives 0.680. Each
test-fold healthy sample is roughly 1/n of the reference it is later scored
against, which pulls it toward z = 0. Same comparison at n=16 gave 0.963 vs
0.786. I had documented a `--cnv-fold-safe` flag in a docstring and never
implemented it; implementing it is how the number got measured.

**End motifs turned out not to be a second finding.** They scored highest at
0.899, and dropped to 0.593 once chromosome entropy was partialled out.
4-mers are read off the reference at fragment ends, so a sample whose coverage
has been redistributed by arm-level gains and losses is sampling a different
mix of sequence context; its spectrum shifts with no change in nuclease
activity at all. One mechanism, two instruments.

**Partialling out aneuploidy manufactured signal.** Removing chromosome
entropy raised fragment-size AUC from 0.489 to 0.713, which looks like hidden
signal emerging. The negative control rose in the same column, 0.507 to 0.606,
and nothing about removing aneuploidy can reveal library signal. Aneuploidy is
a *consequence* of the label, so conditioning on it is conditioning on a
mediator, which induces association. In that column, falls are interpretable
and rises are not. The audit script now uses the negative control to detect
this automatically.

## What the results mean

**Fragment length does not separate these pairs.** Size summary 0.500, size
histogram 0.462, DELFI 0.625. The per-feature paired test is the cleaner
statement: all ten size features split near 50/50 within matched pairs. This
is not "the model failed to find a size signal"; no size feature differs
between a case and its own matched control.

That is a real negative result on the field's best-known fragmentomic
construct, in a cohort where age, sex, smoking, alcohol and batch are held by
matching rather than adjusted away. Three explanations remain open and these
data cannot separate them: 1-3x coverage may be too low for stable per-window
ratios; the cohort is early-stage (9 Stage IA, 9 Stage IIA of 25 breast); or
the within-sample z-scoring I chose discards the global size shift that
carries the signal.

**What works is almost embarrassingly simple.** Chromosome entropy, the
evenness of the per-chromosome fragment distribution, is higher in the case in
**30 of 33 matched pairs** (Wilcoxon p = 3e-07). With coverage dispersion it
gives AUC 0.866 on two features, beating the 1,193-feature model.

The paired test is stronger evidence than the AUC, and this is the part worth
remembering: a signed test over 33 matched pairs makes no use of the
between-sample variance that a small-n AUC is hostage to. The paired figure
shows it, the two group distributions overlap almost completely and nearly
every line still rises.

**The dilution series found the distinction that matters.** Mixing one case's
plasma into its own matched control at known ratios: chromosome entropy is
linear in mixing fraction, R² = 0.998, replicate SD 4.6e-05 or about 3% of the
span. The feature measures tumour content precisely. But against a population
threshold it is not detected at any dilution, because this donor's whole
case-control difference spans only 28% of the between-person range. Against
the same donor's own baseline the limit is 0.08 plasma fraction.

That maps onto the real split in the field between tumour-naive screening and
tumour-informed longitudinal monitoring.

Running all eight liver pairs afterwards changed the number and taught the
more useful lesson. The own-baseline LOD ranges 0.067 to 0.289 across donors,
median 0.108, a 4.3x spread. MAL1246, the pair I happened to pick first, came
in at 0.080, better than typical, and I had picked it for sequencing depth
rather than effect size. The spread tracks how much each tumour shifts the
feature, not stage: the best and worst pairs are both stage 1a. Against a
population threshold, four of the eight were never detected at all, undiluted
included.

I had predicted median 0.08 and a range of 0.03 to 0.20 before running it.
Optimistic at both ends. Worth remembering that the first donor I pick is not
a cohort.

**One caveat I should not forget**, because it is the sort of thing an
interviewer would press on: R² = 0.998 is weaker evidence than it looks.
Mixing two samples that differ in any feature produces a linear response in
that feature by construction. Linearity says the series resolves the mixture,
not that the feature tracks cancer. The second aneuploidy feature proves it:
equally linear, better LOD, and sloping the wrong way for this donor, because
this case is one of the 10 of 33 pairs discordant on coverage dispersion. The
script now checks each slope against the cohort direction.

## What I would not claim

This does not make me a cfDNA assay developer. It is public retrospective
data, n=66, one cohort, one protocol. No alignment from FASTQ, no somatic
variant calling, no tumour-fraction estimation. Cross-cohort generalisation is
untested and the pipeline refuses to produce a number for it.
Leave-one-cancer-type-out is not a substitute: same lab, same five batches.

What it does show is that I can build cfDNA feature extraction and modelling
from primary data, and that I will find and report the ways my own result
might be wrong before someone else does.

## The shortest version

Built a cfDNA cancer classifier on 33 matched case-control pairs. The famous
fragmentomic feature returned nothing. A two-feature aneuploidy proxy beat a
1,193-feature model. The negative control caught a biological feature I had
mistakenly put inside it. A standard shortcut was measured to inflate AUC by
0.155. And the limit-of-detection analysis separated screening from monitoring:
undetectable against a population threshold, median 0.11 plasma fraction
against the patient's own baseline, across eight donors.
