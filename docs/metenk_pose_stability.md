# Met-enkephalin pose-stability campaign — parsed state, 2026-10-05

Five membrane MD systems, one peptide (YGGFM), four independent starting poses
plus a scramble, launched by `run_metenk100.sh` against `config/prod100.mdp`
(100 ns, dt 2 fs, 50 M steps).

**The runs did not finish. All five hit the SLURM wall clock at 20 h on
2026-10-02 07:58 and were cancelled at ~27.6 ns of the 100 ns target.** Nothing
below is a 100 ns result; it is the first quarter of one, from a single replicate
per arm.

Run root: `kuhfeldr@orca:/scratch/kuhfeldr-Kuhfeld_temp/results/04_md/metenk_*`

## Status per arm

| arm | pose source | SLURM | state | reached | target |
|---|---|---|---|---:|---:|
| `metenk_exp` | experimental 8F7Q pose (positive reference) | 193175 | TIMEOUT | 27.66 ns | 100 ns |
| `metenk_hdk` | HADDOCK3 blind, `emref_67` | 193176 | TIMEOUT | 27.60 ns | 100 ns |
| `metenk_chai` | Chai-1 seed0 model1 | 193177 | TIMEOUT | 27.65 ns | 100 ns |
| `metenk_ptx` | Protenix + MSA, seed1 sample0 | 193178 | TIMEOUT | 27.68 ns | 100 ns |
| `metenk_scr` | scrambled GMYGF (negative control) | 193179 | TIMEOUT | 27.66 ns | 100 ns |

Systems are ~102,000 atoms each (receptor 281 res + peptide 5 res, POPC+30 %
chol, OPC water, Na⁺/Cl⁻); the five differ only in the starting pose, as
intended (atom counts 102,007–102,039).

### Why they stopped, and it is not the wall clock alone

Two independent faults compound, and both are in the repo rather than on the
cluster:

1. **mdrun ran on one CPU thread.** `prod.log` for every arm:
   `The number of OpenMP threads was set by environment variable
   OMP_NUM_THREADS to 1`. The profile does request `cpus_per_task: 8`, but
   `rule md_production` in `workflow/rules/04_md.smk` declares no `threads:`,
   so the Snakemake SLURM executor exported `OMP_NUM_THREADS=1`, and
   `04_gmx_run.sh`'s `-ntomp ${OMP_NUM_THREADS:-8}` resolved to 1, not 8.
   Measured throughput was **32.3 ns/day** (from the `npt_10` performance line
   on the same node and thread count), against **137.7 ns/day** for the
   validated 50 ns DAMGO run. The L40S was allocated as asked
   (`NVIDIA L40S, compute cap. 8.9`), so the card was not the problem — the
   L40S-vs-A30 comparison this campaign was also meant to settle is still
   unmeasured.
2. **The walltime request was capped at 20 h.** `rule md_production` asks for
   `60 * production_ns` = 6000 min, but `config/slurm/config.yaml` sets
   `md_production: runtime: 1200`, and the profile wins. 20 h × 32.3 ns/day
   ≈ 27 ns, which is exactly where all five stopped. At the intended 8 threads,
   100 ns would still need ~18 h, i.e. it only just fits in 20 h and has no
   margin.

**The trajectories were deleted.** `prod.xtc` and `prod.tpr` are declared
outputs of a rule that failed, so Snakemake removed them. What survives per arm
is `prod.cpt` (step ~13.83 M), `prod.edr` (50 ps sampling, 553 points) and
`prod.trr` — the full-precision trajectory at `full_traj_ps: 5000`, i.e. **6
frames at a 5 ns stride**. The 100 ps `prod.xtc` that every downstream analysis
expects does not exist.

**A plain rerun will discard the 27.6 ns.** `04_gmx_run.sh` passes
`npt_final.cpt` to `grompp -t` and never passes `-cpi prod.cpt` to `mdrun`, so
re-running the target restarts from the end of equilibration at t=0.

## System stability — all five are healthy

From `prod.edr` over the 27.6 ns that ran (mean ± SD, drift by OLS slope):

| arm | T (K) | P (bar) | density (kg/m³) | box-Z (nm) | E_pot drift (kJ/mol per 10 ns) | LINCS rmsd |
|---|---|---|---|---|---:|---:|
| `metenk_exp` | 310.00 ± 1.24 | 9.3 ± 130.5 | 1022.20 ± 1.90 | 12.431 ± 0.142 | −651 | 4.0e−6 |
| `metenk_hdk` | 310.07 ± 1.32 | 2.3 ± 131.5 | 1021.94 ± 1.79 | 12.439 ± 0.241 | −709 | 4.0e−6 |
| `metenk_chai` | 309.97 ± 1.31 | −0.1 ± 135.5 | 1021.84 ± 1.68 | 12.379 ± 0.189 | −796 | 4.0e−6 |
| `metenk_ptx` | 310.05 ± 1.28 | 2.5 ± 135.8 | 1021.96 ± 1.81 | 12.375 ± 0.258 | −900 | 4.0e−6 |
| `metenk_scr` | 310.02 ± 1.33 | −9.4 ± 145.6 | 1022.38 ± 1.85 | 12.549 ± 0.305 | −708 | 4.0e−6 |

Thermostat on target, barostat centred on zero, constraint rmsd at the
integrator floor, density flat. The residual potential drift is ~0.07 % of
|E_pot| per 10 ns — still relaxing, but nothing is going wrong. **The physics
setup is not what failed here; only the scheduling did.**

## Pose behaviour over the first 25 ns

Peptide heavy-atom RMSD and COM displacement after superposing on receptor CA,
PBC-corrected to the peptide image nearest the receptor centre. Contact =
receptor/peptide heavy atoms within 4.0 Å; *Jaccard* compares the contacted
residue set against t=0.

| arm | RMSD final (Å) | RMSD max (Å) | COM disp. (Å) | receptor CA RMSD (Å) | contact res. at t=0 | contact persistence | Jaccard final |
|---|---:|---:|---:|---:|---:|---:|---:|
| `metenk_exp` | 1.39 | 1.39 | 0.70 | 2.57 | 18 | 0.89 | 0.73 |
| `metenk_ptx` | 1.78 | 1.78 | 0.62 | 1.56 | 15 | 0.84 | 0.76 |
| `metenk_hdk` | 1.89 | 2.02 | 1.13 | 2.29 | 17 | 0.79 | 0.70 |
| `metenk_chai` | 2.06 | 2.14 | 0.71 | 2.00 | 15 | 0.86 | 0.60 |
| `metenk_scr` | 2.76 | 3.34 | 1.80 | 1.87 | 11 | 0.76 | 0.45 |

Per-frame values: `results/04_md_metenk/pose_stability_per_frame.tsv`.
Per-residue contact persistence (packed and auth numbering):
`results/04_md_metenk/pose_contact_persistence.tsv`.

### What this does and does not show

- **Every pose stayed in the pocket.** No arm loses the peptide; even the
  scramble keeps 11→21 contacting residues and a COM displacement under 2 Å.
  "The peptide stayed bound" again fails to discriminate — the same thing the
  pharmacophore work and the MM/GBSA AUROC of 0.483 already said.
- **The ordering is the one you would hope for** (exp < ptx < hdk < chai < scr
  on final RMSD), and all three predicted poses sit within ~0.7 Å of the
  experimental reference. That is consistent with independent pose generators
  converging on the same binding mode under unrestrained physics.
- **But it will not carry weight as written.** n = 1 per arm, 6 frames, 5 ns
  stride, 25 ns of sampling, no error bars. The scramble is separated by about
  1.4 Å of RMSD from the reference — with a single trajectory and a 5 ns stride
  there is no way to say that exceeds run-to-run noise. The receptor itself is
  still drifting (CA RMSD rising monotonically to 2.6 Å in `metenk_exp`), so
  the frame-0 reference is not yet a settled one.
- The Jaccard contact-set overlap is the metric that separates the scramble most
  cleanly (0.33–0.45 across frames, vs 0.52–0.85 for the four real poses). Worth
  keeping as the primary readout if the campaign is re-run.

## Stability is not quality — distance to the experimental pose

The RMSDs above are each arm against **its own** first frame, so they measure
*did this pose move*, not *is this pose right*. Measuring the second question
needs all arms in one reference frame, which is what `md_viz/metenk` provides
(everything superposed on the experimental arm's frame-0 receptor CA).
Peptide heavy-atom RMSD against the experimental pose, no further fitting:

| arm | RMSD vs exp at t=0 | at t=25 ns | COM dist at t=0 | at t=25 ns |
|---|---:|---:|---:|---:|
| `metenk_ptx` | 2.95 Å | 2.99 Å | 1.51 Å | 1.25 Å |
| `metenk_chai` | 4.09 Å | 4.17 Å | 1.99 Å | 1.12 Å |
| `metenk_hdk` | 5.06 Å | 4.97 Å | 1.77 Å | 1.75 Å |
| `metenk_scr` | n/a (GMYGF vs YGGFM) | n/a | 2.09 Å | 2.74 Å |

**No predicted pose is the experimental pose, and 25 ns of unrestrained MD moved
none of them towards it** — 2.95→2.99, 4.09→4.17, 5.06→4.97 Å. Centres of mass
agree to 1–2 Å, so all three are in the orthosteric pocket; the backbone and
side-chain placement inside it differ by 3–5 Å and each arm is sitting stably in
its own distinct local minimum.

This is what the stability numbers actually mean. Low RMSD-vs-own-start is a
statement about how deep a well the pose sits in, not about whether it is the
right well. On this evidence the three generators do **not** converge on a common
binding mode, which is the question the campaign was built to ask.

## What the 3–5 Å gap actually is

Two distances explain the discrepancy, both in
`results/04_md_metenk/anchor_distances.tsv` (packed numbering, auth = +64):

**The message salt bridge is not the problem.** The peptide N-terminal ammonium
sits on the D147 (D3.32) carboxylate at 2.6–3.0 Å in all four real-peptide arms,
for every frame. Every generator gets this right and MD holds it.

| arm | N-term → D147, 0 → 25 ns |
|---|---|
| `metenk_exp` | 2.89 → 2.72 Å (locked throughout) |
| `metenk_ptx` | 2.87 → 2.64 Å (locked throughout) |
| `metenk_chai` | 2.59 → 2.72 Å (locked throughout) |
| `metenk_hdk` | 2.79 → 2.87 Å (locked throughout) |
| `metenk_scr` | **7.35 → 2.84 Å** (forms one it did not start with) |

**The tyrosine phenol is.** Distance from the Tyr hydroxyl to the nearest of
W293 / H297 / Y326:

| arm | Tyr OH → W293/H297/Y326 |
|---|---|
| `metenk_exp` | 5.7–7.0 Å — *not* engaged |
| `metenk_ptx` | 2.8–3.6 Å — H-bonded |
| `metenk_chai` | 2.8–3.6 Å — H-bonded |
| `metenk_hdk` | 3.2–4.3 Å — H-bonded |

**All three generators make the same systematic error**: they anchor the
ammonium correctly and then rotate the phenol into an H-bond with the
W293/H297/Y326 cluster that the experimental structure does not make. That one
difference is most of the 3–5 Å RMSD, and 25 ns of unrestrained MD did not
relax any of them out of it — it is a shared bias of the predictors, not three
independent errors.

**And the scramble is the cautionary note.** It starts with its N-terminus 7.4 Å
from D147 — its free amine belongs to a glycine, and its tyrosine is buried at
position 3 — then swings in over ~20 ns and forms a salt bridge indistinguishable
by distance (2.84 Å) from the native one, with its Tyr already H-bonded to the
anchor cluster from t=0. The pocket will coordinate any free ammonium presented
to it. A contact-distance criterion cannot tell the scramble from the real
peptide here; only the contact-set *identity* (Jaccard) does.

## To re-run

1. Add `threads: 8` to `rule md_production` so the executor stops exporting
   `OMP_NUM_THREADS=1` (`cpus_per_task: 8` alone does not do it).
2. Raise `md_production: runtime` in `config/slurm/config.yaml` above the
   rule's own request, or drop the override; at 8 threads 100 ns needs ~18 h
   plus margin, so `long` (7 d, A30-only) or a chained restart is the safer
   home.
3. Decide whether to resume or restart. Resuming needs `mdrun -cpi prod.cpt`
   wired into `04_gmx_run.sh` and a re-`grompp`ed `prod.tpr` (deterministic from
   the same inputs); without it the 27.6 ns already on disk is thrown away.
4. `full_traj_ps: 5000` gives 20 frames over 100 ns. The 100 ps `reduced_traj_ps`
   xtc is what `04c_md_analysis.py` needs and is what was lost — the frame
   budget should not depend on the run completing.
5. `04c_md_analysis.py` is still DAMGO-specific (`resname DAM`, `reduced.tpr`).
   It cannot analyse these systems; the numbers above come from
   `scripts/contrib/parse_metenk_prod.py`, which selects the peptide by segment
   (`seg_7_system3`).

## Viewing the trajectories

`md_viz/metenk/` holds a receptor-aligned, PBC-corrected viewing set per arm
(receptor + peptide, 4685 atoms, the 6 surviving frames), all superposed on one
common reference so arms can be compared directly.

```
pymol md_viz/metenk/view_ptx_vs_scr.pml   # Protenix vs scramble
pymol md_viz/metenk/view_all5.pml         # all five arms together
```

`mplay` animates; `set all_states, on` overlays all six frames at once, which is
the more honest display — at 5 ns spacing an animation implies a continuity the
data does not have. `overlay_view1.png` / `_view2.png` are rendered overlays of
the Protenix and scrambled arms in the pocket, and
`results/04_md_metenk/metenk_pose_stability.png` plots the three readouts.

Receptor is `resi 1-281`, peptide is `resi 898-902`; auth = packed + 64 in 6DDF
numbering, so the salt-bridge Asp is D147 (D3.32) at packed 83, not the 8F7Q
D149. Both land in chain X through the PDB writer, so `resi` is the only safe
selection handle.
