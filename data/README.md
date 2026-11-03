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

### What the trailing `m` means — RESOLVED

`MAL0638` is the cancer case; `MAL0638m` is that patient's matched control.

Source: Supplementary Table 2 of Zhou et al., *Genome Medicine* 2022
(Additional file 1 on the journal page, **not** the similarly-named
`Supplemantary_Table 2_Hotspots.zip` on Zenodo, which holds per-cohort hotspot
intervals and no sample IDs). Table S2 is laid out one row per pair, with
`Subject_ID` the case and `Control Subject_ID` the matched control; the control
ID is the case ID with `m` appended, and the rule holds for all 33 pairs with
no exceptions. `00_extract_metadata.py` re-verifies it on every run and refuses
to emit labels if it ever fails.

Cohort as recorded in Table S2: 33 pairs — 25 breast, 8 liver. Case and control
ages match exactly, pair by pair.

**Batch structure is the part worth knowing.** Table S2 carries
`Library Preparation Batch` and `Sequencing Batch` (5 batches, Cabat1–5), and
every pair's case and control share both. Batch is therefore orthogonal to
label by construction, not merely claimed to be. That is what makes the
technical-covariate negative control in `02_train_evaluate.py` a real test: in
a cohort where cases and controls were prepared separately, that control fires
by construction and tells you nothing. Here it should sit near chance, and if
it does not, the cause is something other than case/control batch aliasing.

To regenerate:

    python scripts/00_extract_metadata.py --xlsx <Additional file 1>.xlsx

writes `data/metadata/crag_table_s2.tsv` (66 rows: label, cancer type, stage,
grade, age, gender, smoking, alcohol, batch, sequencing ID), which
`00_build_samplesheet.py` then uses as the authoritative label source.

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

    python scripts/00_extract_metadata.py --xlsx <Additional file 1>.xlsx
    python scripts/00_build_samplesheet.py

Required columns are `sample_id`, `path`, `label` (`cancer` or `healthy`) and
`study`. The builder also emits `cancer_type`, `pair_id` and `size_mb`;
`pair_id` lets you keep a case and its matched control in the same CV fold if
you later decide that matters, and the clinical columns from Table S2 feed the
confounder cross-tabs in notebook 01.

Before modelling, sanity-check the pairing with:

    python scripts/00_check_pairing.py --raw data/raw/crag/liver

It compares exact fragment coordinates (never lengths) between each `X`/`Xm`
pair against an unrelated-sample baseline. Table S2 already says they are
different people; this confirms the deposited files match that claim rather
than being two versions of one library.

## Licensing

Code in this repository is MIT. The cfDNA data is not: Zenodo record 6914806 is
CC-BY-4.0 and must be cited as Zhou, Zheng & Liu (2022), with the CRAG paper
itself cited alongside. Do not redistribute downloaded fragment files through
this repo — `data/raw/` is gitignored for that reason as much as for size.
