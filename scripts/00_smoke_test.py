#!/usr/bin/env python3
"""End-to-end smoke test on simulated data. No downloads, ~2 minutes.

Run this first. If it completes, the environment is correct and the pipeline
is wired up; only then is it worth downloading 50 real samples.
"""
import os, subprocess, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cfdna_insight.simulate import simulate_cohort

root = os.path.join(os.path.dirname(__file__), "..")
sim = os.path.join(root, "data", "simulated")
print("[1/3] simulating a 24-sample cohort ...")
sheet = simulate_cohort(sim, n_per_group=12, n_fragments=120_000, effect_size=1.0)
print(f"      sample sheet: {sheet}")

print("[2/3] extracting features ...")
subprocess.run([sys.executable, os.path.join(root, "scripts", "01_extract_features.py"),
                "--samples", sheet, "--outdir", os.path.join(root, "results", "smoke", "features")],
               check=True)

print("[3/3] training and evaluating ...")
subprocess.run([sys.executable, os.path.join(root, "scripts", "02_train_evaluate.py"),
                "--features", os.path.join(root, "results", "smoke", "features"),
                "--outdir", os.path.join(root, "results", "smoke"),
                "--n-repeats", "3", "--n-perm", "30"], check=True)
print("\nSmoke test complete. Results in results/smoke/")
print("REMINDER: these numbers are from simulated data and mean nothing.")
