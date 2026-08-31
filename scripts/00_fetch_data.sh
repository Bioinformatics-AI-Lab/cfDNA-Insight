#!/usr/bin/env bash
# Fetch reference files. cfDNA fragment files must be downloaded manually from
# FinaleDB or Zenodo (see data/README.md) — there is no stable bulk API.
set -euo pipefail

REF_DIR="data/ref"
mkdir -p "$REF_DIR" data/raw

if [[ ! -f "$REF_DIR/hg38.fa" ]]; then
  echo "[fetch] hg38 primary assembly (~3 GB, several minutes)"
  curl -L -o "$REF_DIR/hg38.fa.gz" \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
  gunzip "$REF_DIR/hg38.fa.gz"
  samtools faidx "$REF_DIR/hg38.fa"
else
  echo "[skip] hg38.fa present"
fi

if [[ ! -f "$REF_DIR/hg38-blacklist.v2.bed" ]]; then
  echo "[fetch] ENCODE blacklist v2"
  curl -L -o "$REF_DIR/hg38-blacklist.v2.bed.gz" \
    https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
  gunzip "$REF_DIR/hg38-blacklist.v2.bed.gz"
else
  echo "[skip] blacklist present"
fi

echo
echo "Reference files ready."
echo "Now download cfDNA fragment files (see data/README.md) and write"
echo "config/samples.tsv from config/samples_template.tsv."
