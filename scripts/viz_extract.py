#!/usr/bin/env python3
"""Extract compact structures and profiles for visualisation."""
import json
import pathlib

import numpy as np

ROOT = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp")
OUT = ROOT / "viz"
OUT.mkdir(exist_ok=True)

# ---------------------------------------------------------- 1. complex only
rec = (ROOT / "results/00_receptor/mOR_clean.pdb").read_text().splitlines()
dam = (ROOT / "results/00_receptor/damgo_ref.pdb").read_text().splitlines()
keep = [l for l in rec if l.startswith("ATOM") and (l[76:78].strip() != "H")]
damk = [l for l in dam if l.startswith(("ATOM", "HETATM"))]
# mark DAMGO as HETATM chain L so viewers style it separately
damk = ["HETATM" + l[6:17] + "DAM" + " L" + l[22:] for l in damk]
(OUT / "complex.pdb").write_text("\n".join(keep) + "\nTER\n" + "\n".join(damk) + "\nEND\n")
print(f"complex.pdb: {len(keep)} receptor heavy atoms + {len(damk)} DAMGO")

# ---------------------------------------------------------- 2. membrane slab
sysf = ROOT / "results/03_membrane_prod/membrane_system.pdb"
lines = [l for l in sysf.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]

def z(l):
    return float(l[46:54])

prot, lip, wat, ion = [], [], [], []
for l in lines:
    r = l[17:20].strip()
    if r in ("PC", "PA", "OL", "CHL"):
        lip.append(l)
    elif r in ("WAT", "HOH"):
        wat.append(l)
    elif r.strip() in ("K+", "Cl-"):
        ion.append(l)
    else:
        prot.append(l)
print(f"system: protein {len(prot)}  lipid {len(lip)}  water {len(wat)}  ion {len(ion)}")

# A y-slab through the middle: shows the bilayer in cross-section.
ys = np.array([float(l[38:46]) for l in lines])
lo, hi = -12.0, 12.0
def in_slab(l):
    return lo <= float(l[38:46]) <= hi

slab = [l for l in prot]                                  # all protein
slab += [l for l in lip if in_slab(l)]                    # lipids in slab
slab += [l for l in wat if in_slab(l)][::4]               # thin the water
slab += [l for l in ion if in_slab(l)]
slab = [l for l in slab if l[76:78].strip() != "H"]       # heavy atoms only
(OUT / "slab.pdb").write_text("\n".join(slab) + "\nEND\n")
print(f"slab.pdb: {len(slab)} atoms (y in [{lo},{hi}], heavy atoms, water thinned 4x)")

# ---------------------------------------------------------- 3. density profile
bins = np.arange(-70, 60.5, 2.0)
centres = (bins[:-1] + bins[1:]) / 2
def hist(sel):
    if not sel:
        return [0] * len(centres)
    h, _ = np.histogram([z(l) for l in sel], bins=bins)
    return h.tolist()

prof = {
    "z": centres.tolist(),
    "protein": hist(prot),
    "lipid_tails": hist([l for l in lip if l[17:20].strip() in ("PA", "OL")]),
    "lipid_heads": hist([l for l in lip if l[17:20].strip() == "PC"]),
    "cholesterol": hist([l for l in lip if l[17:20].strip() == "CHL"]),
    "water": hist(wat),
    "ions": hist(ion),
}
dz = [z(l) for l in lines if l[17:20].strip() == "DAM"]
prof["damgo_z"] = [min(dz), max(dz)] if dz else None
(OUT / "density.json").write_text(json.dumps(prof))
print(f"density.json written; DAMGO z-range {prof['damgo_z']}")

# ---------------------------------------------------------- 4. Stage 2 poses
pep = ROOT / "results/02_cofold"
for name in ("real__met_enkephalin", "real__casoxin_C"):
    got = []
    for sd in sorted((pep / name).glob("seed*")):
        c = sorted(sd.glob("pred.model_idx_*.cif"))
        if c:
            got.append((sd.name, c[0]))
    if got:
        print(f"{name}: {len(got)} seed models available")
        (OUT / f"{name}_seeds.txt").write_text(
            "\n".join(str(p) for _, p in got) + "\n")
print("\nWrote", OUT)
