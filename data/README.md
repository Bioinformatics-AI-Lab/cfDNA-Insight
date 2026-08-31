# Data acquisition

## Why fragment files rather than BAMs

Almost all cfDNA WGS is deposited in dbGaP or EGA under controlled access,
because raw reads carry recoverable genotype. Access takes months and usually a
DTA between institutions.

Fragmentation analysis does not need the sequence — only where each fragment
started and ended. FinaleDB exploits this: it re-processed published cfDNA WGS
through one uniform pipeline, discarded the sequence, and released the
coordinates. That makes the data openly downloadable and, as a side benefit,
removes the aligner-choice component of between-study batch effects.

FinaleDB's processing, for the record, since it constrains what this repo can
claim: reads trimmed to 50 bp at the 3' end, Trimmomatic adapter/quality
trimming, BWA-MEM against both hg19 and hg38, samblaster duplicate marking,
then filtering to properly-paired primary alignments with MAPQ >= 30 and no
duplicate flag.

The 3'-end trimming is worth understanding rather than just accepting. Fragment
length is inferred from the outer 5' ends of the read pair; trimming the 3' end
of each read removes bases on the *inner* side, so the outer coordinates and
therefore the inferred fragment length are preserved. Trimming to a common
length across studies is what makes cross-study comparison defensible at all.

## Primary source

FinaleDB portal: http://finaledb.research.cchmc.org/

Query by pathological condition, select samples, add to the download list. If
the portal is unavailable, use the Zenodo mirrors below.

## Zenodo mirrors

    # Cristiano 2019 breast cancer + matched healthy controls
    mkdir -p data/raw/cristiano2019
    # https://zenodo.org/records/6914806  -> breast.tar

    # Snyder 2016 BH01: single deeply sequenced healthy plasma sample.
    # Useful as the healthy background for the dilution series because its
    # depth means mixing does not exhaust it.
    # https://zenodo.org/record/5227520

Download through the browser or with curl once you have the direct file URLs
from the record page.

## Reference files

    # hg38 primary assembly (needed for GC correction and end motifs)
    mkdir -p data/ref && cd data/ref
    curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
    gunzip GRCh38.primary_assembly.genome.fa.gz
    mv GRCh38.primary_assembly.genome.fa hg38.fa
    samtools faidx hg38.fa

    # ENCODE blacklist v2
    curl -LO https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
    gunzip hg38-blacklist.v2.bed.gz

Note for Apple Silicon: `samtools` and `pysam` both install cleanly from
conda-forge/bioconda as native arm64 builds. Unlike macs2, no Rosetta shim is
needed anywhere in this pipeline.

## Sample sheet

`config/samples.tsv`, tab-separated, four required columns:

    sample_id    path                                  label     study
    S001         data/raw/cristiano2019/S001.frag.gz   cancer    cristiano2019
    S002         data/raw/cristiano2019/S002.frag.gz   healthy   cristiano2019

`label` must be exactly `cancer` or `healthy`. `study` is not optional and is
not cosmetic: it defines the groups for leave-one-study-out evaluation. If every
sample shares one study value, the pipeline will say so and refuse to report a
cross-study number rather than quietly reporting a within-cohort one.

Record for each sample, in `data/metadata/`, whatever the source study
published: stage, age, sex, sequencing platform, collection site. Anything that
differs systematically between cancer and healthy samples is a candidate
confounder and belongs in the limitations section.

## Licensing

Code in this repository is MIT. The cfDNA datasets are not: each carries the
terms of its originating study and of FinaleDB. Cite the original studies, not
just FinaleDB. Do not redistribute downloaded fragment files through this repo —
`data/raw/` is gitignored for that reason as much as for size.
