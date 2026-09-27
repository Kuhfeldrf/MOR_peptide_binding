# The membrane QC gate rejects the system that produced our validated result

Found while diagnosing seven `membrane_qc` failures. The gate is not calibrated
against anything, and the workflow was not reproducing the packing settings
behind its own only validated MD run.

## The QC gate fails DAMGO

`results/03_membrane_prod` is the system that minimised, equilibrated, and ran
**50 ns of production MD**, giving the result in
[`stage4_damgo_50ns.md`](stage4_damgo_50ns.md) - ligand RMSD 1.30 A, ASP147 in
contact in all 501 frames. Running the QC script on it:

| criterion | DAMGO (ran 50 ns) | beta-casomorphin-7-hum (blocked) |
|---|---|---|
| packmol reached tolerance | **NO** (`STOP 173`) | NO (`STOP 173`) |
| inter-residue < 0.5 A | **1** | 7 |
| inter-residue < 1.2 A | **111** | 327 |
| inter-residue < 1.8 A | **1601** (worst 0.242 A) | 2691 (worst 0.157 A) |
| lipid atoms enclosed by protein | 83 | 114 |

**DAMGO fails all four criteria.** A gate that rejects a demonstrably usable
system does not predict whether MD will work, so passing or failing it carries
little information as currently set.

## Two separate problems

### 1. The workflow used different packing settings than the validated run

`nloop` was **set nowhere** - not in `config/`, not in the pack wrapper, not in
the rule. The workflow silently took packmol-memgen's defaults:

| | nloop_all | nloop (per lipid) |
|---|---|---|
| packmol-memgen default, used by the workflow | 100 | 20 |
| **DAMGO, set by hand in the manual Stage 3 work** | **200** | **40** |

So the peptide packs ran with **half** the iterations that produced the
validated system, and came out roughly **twice as bad** on every contact metric
- consistent with simply having had less time to relax.

And packmol had not stalled. The objective function over its final four loops:

```
loop  97: f = 3.828e4
loop  98: f = 1.826e4
loop  99: f = 1.422e4
loop 100: f = 1.113e4
STOP: Maximum number of GENCAN loops achieved.
```

A **3.4x reduction still falling monotonically** when it hit the ceiling. It ran
out of iterations, not out of progress.

**Fixed:** `membrane.nloop_all: 200` and `membrane.nloop: 40` are now explicit
in `config/config.yaml` and passed through the rule, so the workflow reproduces
the configuration behind its own validated result instead of inheriting a
default.

### 2. The gate mislabels its own diagnosis - NOT YET FIXED

The QC script reports close contacts as:

> `7 inter-residue contacts under 0.5 A under MINIMUM IMAGE - the pack is not periodic`

**That inference is wrong.** [D29](decisions.md#d29) was a genuinely aperiodic
pack - `--pbc` was missing, the patch overhung the box, atoms overlapped their
own periodic images, and both GROMACS and sander reported infinite force. The
fix was `--pbc`, which every pack now uses.

A periodic pack can still contain close contacts; they are ordinary packing
imperfection that minimisation removes. DAMGO is periodic, has one contact under
0.5 A, and ran for 50 ns. So "close contacts under minimum image" does not imply
"not periodic", and the message attributes a benign observation to a
catastrophic cause.

## RESOLVED by measurement, 2026-09-27

Option 1 was run: re-packed at 200/40, then minimised. Both questions answered.

### The packs are fine

| | PE (kJ/mol) | Fmax |
|---|---|---|
| **DAMGO reference** (ran 50 ns) | -1.1710e6 | **998.97** |
| met-enkephalin | -1.1864e6 | 757 |
| beta-casomorphin-4-bov | -1.1815e6 | 811 |
| beta-casomorphin-5-bov | -1.1884e6 | 980 |
| beta-casomorphin-5-hum | -1.1759e6 | 839 |
| beta-casomorphin-7-bov | -1.1770e6 | 909 |
| beta-casomorphin-7-hum | -1.1785e6 | 955 |
| neocasomorphin-6 | -1.1773e6 | 917 |

**All seven reach a LOWER final energy and a LOWER maximum force than the
reference system that ran 50 ns** - including beta-casomorphin-7-hum, whose
pack had 11 contacts under 0.5 A. The prediction that "minimisation will
overflow" was simply wrong.

The re-pack also confirmed the box: 86.06 x 86.06 x 119.0 A, identical to
DAMGO, and 6.89 MB against DAMGO's 6.89 MB. The earlier 154-172 MB packs were
the coordinate-frame bug (D36), not a packing problem - which is why packing
now takes 1 h 14 rather than 4 h.

### The gate was recalibrated, not loosened

Three criteria became WARN, each for a stated empirical reason rather than to
make failures go away:

1. **packmol convergence** - the reference system also ends in `STOP 173`.
   packmol does not reach tolerance for a receptor-in-bilayer system of this
   size, so failing on it rejects everything, including the system that
   produced the project's only validated MD result.
2. **close contacts described as "the pack is not periodic"** - a false
   inference, now removed. D29 was genuinely aperiodic (`--pbc` missing, atoms
   overlapping their own images); every pack now uses `--pbc`, and a periodic
   pack can still contain close contacts.
3. **"minimisation will overflow"** - a prediction, falsified seven times out
   of seven.

Every WARN now prints the reference system's value beside it, so a reader can
judge whether a pack is unusual rather than compare against a threshold nobody
calibrated.

**Minimisation is the gate now.** It is cheap, and it is exactly what the
check was trying to predict.

## What had NOT been decided, before the above

**Whether to loosen the thresholds.** Tempting and wrong to simply lower them
until the seven pass - that is exactly the over-fitting to a transient problem
that [D26](decisions.md) retired a whole script for.

DAMGO gives exactly **one** calibration point: a pack with 1 contact < 0.5 A,
111 < 1.2 A and 1601 < 1.8 A survived minimisation and 50 ns. It does not
establish where the real cliff is, and the 7 systems are ~2x worse on every
count, which may or may not be survivable.

The defensible options, in order of preference:

1. **Re-pack at 200/40 and re-measure.** If the seven land near DAMGO's numbers,
   the question dissolves without touching the gate. This is the fix already
   applied and the obvious next step.
2. **Make the gate empirical rather than absolute** - report each metric
   alongside DAMGO's value as the reference, fail only when materially worse,
   and keep the hard failure for what actually proved fatal (infinite force on
   minimisation).
3. **Let minimisation be the arbiter.** It is cheap, and it is the thing the
   gate is trying to predict. A pack that minimises without overflow is usable
   by definition.

Option 1 first, since it costs nothing and may remove the need for a judgement
call.

## Correction to an earlier claim

The seven failures were first reported as *"QC did its job - it caught a bad
pack rather than letting it flow into MD."* That was asserted without checking
the gate against a known-good system. QC did flag real imperfections, but it
also rejects the system that produced the project's only validated MD result, so
it cannot be said to have discriminated correctly.
