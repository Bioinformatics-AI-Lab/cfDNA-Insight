# Biological background

## Where cfDNA comes from

Circulating cell-free DNA is released into plasma mainly by apoptotic cells.
In a healthy adult the overwhelming majority is haematopoietic — neutrophils,
lymphocytes, monocytes — with a smaller contribution from vascular endothelium
and hepatocytes. Total concentration is typically a few nanograms per millilitre
of plasma, and the half-life is short, on the order of tens of minutes to a
couple of hours. That short half-life is what makes cfDNA a near-real-time
readout rather than a historical one.

In a cancer patient some fraction of that pool is tumour-derived (ctDNA). The
fraction is the central difficulty of the entire field: it ranges from tens of
percent in advanced disease down to well below 0.1% in early-stage,
low-shedding tumours. **The clinically valuable regime — early detection — is
exactly the regime where the signal is smallest.** Any method that reports
excellent performance on late-stage cohorts has demonstrated very little about
the screening use case.

## Why fragmentation is informative at all

cfDNA is not sheared at random. Before an apoptotic cell's DNA is released, it
is digested by nucleases, and the nucleosome protects the ~147 bp it wraps.
What survives is a fragment length distribution dominated by a peak near 167 bp
(147 bp of core particle plus roughly 20 bp of linker), with a smaller
dinucleosomal peak near 334 bp.

Superimposed on the main peak is a ~10 bp periodicity. DNA on the histone
octamer has a helical pitch of about 10.4 bp, so the minor groove faces outward
periodically and nuclease access oscillates with the same period. Its presence
is a signature of genuine nucleosome-protected DNA; its absence means the sample
is degraded.

Because nucleosome positioning tracks chromatin state, and chromatin state is
cell-type specific, the *positions* at which cfDNA fragments start and end carry
information about which tissue released them. This is the conceptual bridge from
chromatin accessibility work to cfDNA: both are readouts of nucleosome
occupancy, one by transposase insertion preference and the other by nuclease
protection. The physical quantity being measured is closely related; the assay
that measures it is inverted.

## What changes in cancer

**Fragments get shorter.** Tumour-derived cfDNA is consistently shifted toward
shorter lengths, with enrichment below 150 bp. The mechanism is not fully
settled — proposed contributors include altered chromatin compaction in tumour
cells, different nuclease exposure, and shorter effective protection. The
observation is robust across studies even where the mechanism is disputed.

**The fragmentation profile becomes uneven across the genome.** This is the
DELFI observation: in healthy plasma the short/long ratio is relatively uniform
across megabase windows, because the pool is dominated by one fairly homogeneous
cell population. In cancer plasma it becomes irregular, reflecting a mixture of
haematopoietic and tumour chromatin plus large-scale copy number changes. The
diagnostic signal is in the *variance and shape* of the profile, not its mean.

**Copy number becomes unbalanced.** Most solid tumours are aneuploid, so
tumour-derived fragments over-represent gained arms and under-represent lost
ones. At sufficient tumour fraction this is visible as arm-level coverage
deviation in low-coverage WGS.

**The nuclease repertoire shifts.** DNASE1L3 activity is reduced in several
cancers, and different nucleases have different sequence preferences at the cut
site. This shows up as an altered spectrum of 4-mers at fragment 5' ends and a
narrowed motif diversity score.

These four are not independent. Aneuploidy changes coverage, coverage changes
the counts underlying the short/long ratio, and a global size shift moves the
whole ratio distribution. The feature-block ablation exists partly to quantify
how much of that overlap is real redundancy.

## Confounders that are easy to mistake for biology

**Pre-analytical variation.** Time from venepuncture to plasma separation is
probably the single largest non-biological driver of cfDNA fragment length: as
white cells lyse in the tube, they release genomic DNA that is longer and less
nucleosome-protected. EDTA versus Streck tubes, storage temperature and
freeze-thaw cycles all move the size distribution. Cancer samples are usually
collected prospectively in a clinic; healthy controls are often convenience
samples from a biobank. If those two collection pathways differ, a size-based
classifier will separate them beautifully and will be measuring logistics.

**Library preparation and sequencing.** PCR cycle number changes GC bias.
Adapter chemistry changes the recoverable size range. Read length changes what
counts as a properly-paired fragment. When cases and controls are prepared in
different batches, all three become label-correlated.

**Clonal haematopoiesis.** Mutations from expanded haematopoietic clones appear
in plasma and are indistinguishable from tumour mutations without a matched
buffy coat. CHIP prevalence rises steeply with age. This is the main reason
mutation-based features are excluded from this project entirely: at low
coverage they would be noise, and at high coverage they would be confounded.

**Age and sex.** cfDNA concentration and composition both shift with age, and
cancer cohorts are almost always older than healthy control cohorts. Sex drives
obvious separation via chrX/chrY, which is why those chromosomes are excluded
here — but sex-linked differences in autosomal cfDNA composition exist too.

**Treatment and comorbidity.** Chemotherapy causes massive cell death and
raises cfDNA concentration. Inflammation, infection, recent surgery, exercise
and pregnancy all alter the cfDNA pool substantially, and some of them alter it
in the same direction cancer does.

The practical consequence for this project: whatever the classifier achieves,
the interesting question is which of these it might be reading. That is what the
technical-covariate negative control, the leave-one-study-out evaluation and the
metadata recorded in `data/metadata/` are for.

## Selected reading

- Snyder et al., *Cell* 2016 — nucleosome footprints in cfDNA reveal tissue of origin
- Cristiano et al., *Nature* 2019 — DELFI, genome-wide fragmentation profiles
- Jiang et al., *Cancer Discovery* 2020 — end-motif profiling
- Mouliere et al., *Science Translational Medicine* 2018 — size selection enriches ctDNA
- Adalsteinsson et al., *Nature Communications* 2017 — ichorCNA
- Lo, Han, Jiang & Chiu, *Science* 2021 — review of cfDNA biology and fragmentomics
