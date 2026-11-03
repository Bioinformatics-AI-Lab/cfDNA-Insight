# Data acquisition

## Why fragment files rather than BAMs

Almost all cfDNA WGS is deposited in dbGaP or EGA under controlled access,
because raw reads carry recoverable genotype. Access takes months and usually a
DTA between institutions.

Fragmentation analysis does not need the sequence — only where each fragment
started and ended. Several groups exploit this: they re-process published cfDNA
WGS through one uniform pipeline, discard the sequence, and release the
coordinates. That makes the data openly downloadable and, as a side benefit,
removes the aligner-choice component of between-study batch effects.

## Primary source: the CRAG validation cohort

**Zenodo record 6914806**, CC-BY-4.0, from Zhou, Zheng, Fu et al., *Genome
Medicine* 2022 (CRAG). This holds the independent validation cohort the CRAG
authors generated themselves, as `*.hg19.frag.bed.gz` files, genome-wide.

| Archive | Pairs | Size |
|---|---|---|
| `liver.tar` | 8 | 3.3 GB |
| `breast_1.tar` … `breast_4.tar` | 25 | 11.0 GB |

33 tumour plasma samples (25 breast, 8 liver) from the UC Biorepository, each
with a control drawn from the UC Fernald Community Cohort and matched on age,
gender, smoking history and alcohol use.

**The property that makes this cohort worth using:** cases and their matched
controls were processed in the same batch. Most public cfDNA case-control
comparisons collect patients prospectively in a clinic and controls from a
biobank, which aliases label with collection pathway and library prep. Here it
does not. The technical-covariate negative control in `02_train_evaluate.py`
should therefore sit near chance — and if it doesn't, that is a real finding
rather than the expected artefact.

Fetch with:

    bash scripts/00_fetch_crag.sh liver      # start here — 3.3 GB, 8 pairs
    bash scripts/00_fetch_crag.sh breast     # 11.0 GB, 25 pairs

### Open question: what the trailing `m` means

Every sample appears twice — `MAL0638.hg19.frag.bed.gz` and
`MAL0638m.hg19.frag.bed.gz`. The counts match the paper's cohort exactly
(25 + 25 breast, 8 + 8 liver), so the pair is almost certainly case and matched
control. **Which member is which is not documented** — not in the paper, not in
the Zenodo record, not in the file names.

Getting it backwards silently inverts every label in the study.
`00_build_samplesheet.py` therefore requires `--m-is control|cancer` and will
not guess.

Do **not** settle this by comparing fragment lengths and assuming the shorter
group is cancer. That is the hypothesis under test; using it to assign labels
guarantees a confirmation. Settle it from Additional file 1 Table S2 of the
paper, or by emailing the corresponding author (lyping1986@gmail.com).

**Resolved as of (date): _______ — `m` = _______ , source: _______**

### Study structure, and what can honestly be claimed

This is one cohort, one lab, one protocol. Leave-one-study-out does not apply,
and `config/samples.tsv` reflects that: `study` is `crag_validation` for every
row, so the pipeline will correctly report cross-study generalisation as
untested.

`config/samples_by_cancer_type.tsv` sets `study` to `crag_breast` /
`crag_liver`. Running against that sheet performs **leave-one-cancer-type-out**:
train on breast plus its controls, test on liver plus its controls. Weaker than
a true cross-cohort test, since protocol and lab are shared, but it does test
whether the signal transfers across tumour biology. Report it under that name,
never as cross-study validation.

## FinaleDB

[FinaleDB](http://finaledb.research.cchmc.org/) hosts uniformly processed cfDNA
WGS from several published studies and is the natural source for a genuine
second cohort.

**As of 2026-11-02 its backend is down.** The React frontend loads, but the only
two live API routes both return 500, so the portal renders blank:

    GET /api/v1/summary     -> 500
    GET /api/v1/publication -> 500

Worth re-checking periodically. If it comes back, pulling a second accession
with both cancer and healthy samples would upgrade the leave-one-cancer-type-out
analysis to a real leave-one-study-out.

Note for whenever it returns: the portal has two independent selection
mechanisms that do not share state. Row checkboxes feed the *Visualize* button
only; the **`+` at the right end of each row** builds the download list. Ticking
checkboxes leaves the Download button greyed out.

FinaleDB's processing, for the record, since it constrains cross-source
comparison: reads trimmed to 50 bp at the 3' end, Trimmomatic adapter/quality
trimming, BWA-MEM against both hg19 and hg38, samblaster duplicate marking, then
filtering to properly-paired primary alignments with MAPQ >= 30. The 3'-end
trimming preserves fragment length — length is inferred from the outer 5' ends
of the pair, and trimming removes inner bases.

## Reference files

    bash scripts/00_fetch_data.sh

Pulls UCSC hg19 (~1 GB compressed) and the ENCODE hg19 blacklist v2. **The build
is hg19**, matching the CRAG files. Do not mix builds in one sample sheet — bin
coordinates would not correspond, and the pipeline has no way to detect it.

Note for Apple Silicon: `samtools` and `pysam` both install cleanly from
conda-forge/bioconda as native arm64 builds. No Rosetta shim is needed anywhere
in this pipeline.

## Sample sheet

    python scripts/00_build_samplesheet.py --m-is control

Required columns are `sample_id`, `path`, `label` (`cancer` or `healthy`) and
`study`. The builder also emits `cancer_type`, `pair_id` and `size_mb`;
`pair_id` lets you keep a case and its matched control in the same CV fold if
you later decide that matters.

## Licensing

Code in this repository is MIT. The cfDNA data is not: Zenodo record 6914806 is
CC-BY-4.0 and must be cited as Zhou, Zheng & Liu (2022), with the CRAG paper
itself cited alongside. Do not redistribute downloaded fragment files through
this repo — `data/raw/` is gitignored for that reason as much as for size.
