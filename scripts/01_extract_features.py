#!/usr/bin/env python3
"""Extract per-sample feature blocks from fragment files.

    python scripts/01_extract_features.py \
        --samples config/samples.tsv \
        --outdir results/features \
        --fasta data/ref/hg38.fa

The FASTA is optional. Without it, GC correction and end-motif features are
skipped and the run is logged as `gc_corrected=False` so the provenance of any
downstream result is unambiguous.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from cfdna_insight.features import (  # noqa: E402
    compute_bin_gc, extract_sample_features, save_features, stack_blocks,
)
from cfdna_insight.io import SampleSheet  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--outdir", default="results/features")
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--genome", default="hg38", choices=["hg38", "hg19"])
    ap.add_argument("--bin-size", type=int, default=100_000)
    ap.add_argument("--window-factor", type=int, default=50)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--gc-cache", default="results/cache/bin_gc.parquet")
    args = ap.parse_args()

    sheet = SampleSheet.load(args.samples).existing()
    print(f"[info] {len(sheet)} samples")
    print(sheet.table["label"].value_counts().to_string())

    gc_track = None
    if args.fasta:
        print("[info] computing per-bin GC (cached after first run)")
        gc_track = compute_bin_gc(args.fasta, bin_size=args.bin_size,
                                  genome=args.genome, cache=args.gc_cache)
    else:
        print("[warn] no --fasta: skipping GC correction and end motifs")

    per_sample, failures = {}, []
    t0 = time.time()
    for i, row in sheet.table.iterrows():
        try:
            per_sample[row["sample_id"]] = extract_sample_features(
                row["path"], row["sample_id"], gc_track=gc_track,
                fasta_path=args.fasta, bin_size=args.bin_size,
                window_factor=args.window_factor, genome=args.genome,
                min_mapq=args.min_mapq,
            )
            print(f"  [{i + 1}/{len(sheet)}] {row['sample_id']} ok")
        except Exception as exc:  # noqa: BLE001
            print(f"  [{i + 1}/{len(sheet)}] {row['sample_id']} FAILED: {exc}")
            failures.append({"sample_id": row["sample_id"], "error": str(exc)})

    if not per_sample:
        print("[error] no samples processed")
        return 1

    blocks = stack_blocks(per_sample)
    save_features(blocks, args.outdir)

    manifest = {
        "n_samples": len(per_sample),
        "n_failed": len(failures),
        "failures": failures,
        "gc_corrected": bool(args.fasta),
        "end_motifs": bool(args.fasta),
        "genome": args.genome,
        "bin_size": args.bin_size,
        "window_size": args.bin_size * args.window_factor,
        "min_mapq": args.min_mapq,
        "blocks": {k: list(v.shape) for k, v in blocks.items()},
        "elapsed_sec": round(time.time() - t0, 1),
    }
    with open(os.path.join(args.outdir, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)

    sheet.table.to_csv(os.path.join(args.outdir, "samples_used.tsv"),
                       sep="\t", index=False)
    print(f"[done] {json.dumps(manifest['blocks'])} -> {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
