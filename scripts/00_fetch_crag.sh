#!/usr/bin/env bash
# Fetch the CRAG independent-validation cfDNA fragment files from Zenodo.
#
#   bash scripts/00_fetch_crag.sh liver      # 3.3 GB, 8 pairs  (start here)
#   bash scripts/00_fetch_crag.sh breast     # 11.0 GB, 25 pairs
#   bash scripts/00_fetch_crag.sh all        # 14.3 GB
#
# Zenodo record 6914806 (CRAG, Zhou/Zheng/Liu, Genome Medicine 2022),
# CC-BY-4.0. Files are *.hg19.frag.bed.gz, genome-wide.
#
# Per the paper, this cohort is 33 tumour plasma samples (25 breast, 8 liver)
# from the UC Biorepository, each with a control from the UC Fernald Community
# Cohort matched on age, gender, smoking and alcohol use — and cases and their
# matched controls were processed in the same batch. That last property is why
# this cohort is worth using: it removes most of the case/batch aliasing that
# the technical-covariate control exists to detect.
set -euo pipefail

TARGET="${1:-liver}"
DEST="data/raw/crag"
BASE="https://zenodo.org/records/6914806/files"

case "$TARGET" in
  liver)  FILES=(liver.tar) ;;
  breast) FILES=(breast_1.tar breast_2.tar breast_3.tar breast_4.tar) ;;
  all)    FILES=(liver.tar breast_1.tar breast_2.tar breast_3.tar breast_4.tar) ;;
  *) echo "usage: $0 [liver|breast|all]" >&2; exit 1 ;;
esac

mkdir -p "$DEST"
cd "$DEST"

for f in "${FILES[@]}"; do
  if [[ -f "${f%.tar}/.done" ]]; then
    echo "[skip] $f already extracted"
    continue
  fi
  echo "[fetch] $f"
  curl -L -C - -o "$f" "$BASE/$f?download=1"
  echo "[extract] $f"
  tar xf "$f"
  touch "${f%.tar}/.done"
  rm -f "$f"          # the tar is redundant once extracted; saves ~14 GB
done

echo
echo "Fragment files in $DEST:"
find . -name '*.frag.bed.gz' | wc -l
echo "Next: python scripts/00_build_samplesheet.py"
