#!/usr/bin/env python3
"""Chai-1 smoke test: one peptide, one seed.

Purpose is to trigger the weight download, confirm the API call shape, record
the output layout and a timing figure, before committing to the full
5-seed x N-peptide run. Not part of the pipeline.
"""
import pathlib
import time

import torch
from chai_lab.chai1 import run_inference

ROOT = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp")
OUT = ROOT / "results/02_cofold/_smoke"
OUT.mkdir(parents=True, exist_ok=True)

receptor = (ROOT / "results/00_receptor/mOR_clean.fasta").read_text().split("\n")[1].strip()
peptide = "YGGFM"  # Met-enkephalin

fa = OUT / "smoke.fasta"
fa.write_text(
    f">protein|name=receptor\n{receptor}\n"
    f">protein|name=peptide\n{peptide}\n"
)
print(f"receptor {len(receptor)} aa + peptide {len(peptide)} aa "
      f"= {len(receptor)+len(peptide)} tokens")
print("GPU:", torch.cuda.get_device_name(0),
      "| bf16:", torch.cuda.is_bf16_supported())

t0 = time.time()
out = run_inference(
    fasta_file=fa,
    output_dir=OUT / "run_seed0",
    use_msa_server=False,        # single-sequence mode, per the instructions
    use_esm_embeddings=True,
    num_trunk_recycles=3,
    num_diffn_timesteps=200,
    seed=0,
    device="cuda:0",
)
dt = time.time() - t0
print(f"\n=== inference finished in {dt:.1f} s ===")
print("returned type:", type(out))
try:
    print("attrs:", [a for a in dir(out) if not a.startswith('_')][:40])
except Exception as e:
    print("introspect failed:", e)

print("\n=== output tree ===")
for p in sorted((OUT / "run_seed0").rglob("*")):
    if p.is_file():
        print(f"  {p.relative_to(OUT / 'run_seed0')}  ({p.stat().st_size} B)")
