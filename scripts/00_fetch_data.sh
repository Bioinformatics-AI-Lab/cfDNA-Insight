#!/usr/bin/env bash
# Fetch hg19 reference files.
#
# Build is hg19 because the CRAG validation-cohort fragment files are
# *.hg19.frag.bed.gz. The cfDNA fragment tars themselves are fetched by
# 00_fetch_crag.sh.
set -euo pipefail

REF_DIR="data/ref"
mkdir -p "$REF_DIR" data/raw

if [[ ! -f "$REF_DIR/hg19.fa" ]]; then
  echo "[fetch] hg19 primary assembly (~1 GB compressed, several minutes)"
  curl -L -o "$REF_DIR/hg19.fa.gz" \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
  gunzip "$REF_DIR/hg19.fa.gz"
  if command -v samtools >/dev/null 2>&1; then
    samtools faidx "$REF_DIR/hg19.fa"
  else
    echo "[note] samtools not found; pysam will build the index on first use"
  fi
else
  echo "[skip] hg19.fa present"
fi

if [[ ! -f "$REF_DIR/hg19-blacklist.v2.bed" ]]; then
  echo "[fetch] ENCODE blacklist v2 (hg19)"
  curl -L -o "$REF_DIR/hg19-blacklist.v2.bed.gz" \
    https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
  gunzip "$REF_DIR/hg19-blacklist.v2.bed.gz"
else
  echo "[skip] blacklist present"
fi

echo
echo "Reference files ready (hg19)."
echo "Next: bash scripts/00_fetch_crag.sh to pull the cfDNA fragment files."
