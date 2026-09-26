# Is ABFE the right method for these peptides?

Asked before implementing Stage 6, because the answer changes what gets built.
Short version: **ABFE is defensible for DAMGO and marginal for the
casomorphins. It is the wrong tool for the 10-20+ residue peptides this
project ultimately wants to screen, and scaling it there would be a
methodological error a reviewer would catch.**

## The DAMGO reference value - NOT YET USABLE

A computed ΔG needs an experimental number to be validated against. The
reference table currently has `ki_status = NOT_REPORTED` for every row, and
DAMGO is not in it at all.

| | |
|---|---|
| Value found | **Ki = 1.23 nM** (→ ΔG ≈ -12.2 kcal/mol at 310 K) |
| Attributed to | "Binding affinity and selectivity of opioids at mu, delta and kappa receptors in monkey brain membranes" (J Pharmacol Exp Ther) |
| Status | **UNCONFIRMED** |

Two reasons it is not yet usable as the yardstick:

1. **Not read in full text.** The value came from a search summary, not from
   the paper. Under the project's own sourcing rule that is not a source.
2. **Wrong preparation.** It is monkey brain membrane, whereas the receptor
   modelled here is **human** MOR (6DDF). Opioid affinities differ between
   native tissue and cloned receptor preparations, so the assay system has to
   match the claim.

The right citation is likely Raynor et al. 1994, *Mol Pharmacol* 45:330-334
([PMID 8114680](https://pubmed.ncbi.nlm.nih.gov/8114680/)), which characterised
the **cloned** receptors - but the DAMGO value was not extracted from it, so it
is recorded here as a lead, not a source.

**Action required before any ΔG is called validated:** obtain the full text,
extract the value with its assay conditions, and record it with its PMID.

## Why ABFE degrades with peptide length

Alchemical ABFE was designed for small, rigid, drug-like ligands. Five
specific things break as the ligand grows:

1. **Conformational entropy.** The free peptide samples an enormous ensemble;
   the bound one is restricted. ΔG_bind contains that entropy change, so BOTH
   ensembles must be sampled. Cost grows roughly exponentially with the number
   of rotatable torsions - about 2 backbone dihedrals per residue plus side
   chains.

2. **Boresch restraints assume ONE bound pose.** The six restrained degrees of
   freedom carry an *analytical* standard-state correction that is valid only
   if the ligand occupies a single well-defined basin. **Our β-casomorphins
   show 3.5-7.5 Å RMSD across co-folding seeds** - there is no single pose. The
   correction would then be wrong, not merely imprecise.

3. **Charged-solute finite-size artifacts.** Decoupling a net-charged ligand
   under PBC needs explicit corrections. Peptides are frequently charged, and a
   20-mer can carry several formal charges.

4. **Timescale mismatch.** Peptide conformational rearrangement in a pocket
   runs to µs-ms. A 5 ns window samples a sliver of that.

5. **The decoupled state is pathological.** At λ→0 the peptide is effectively a
   gas-phase chain, free to collapse or extend into conformations it would
   never adopt in either real end state. Recoupling that into a binding site is
   a severe sampling problem that barely exists for a rigid small molecule.

Published peptide ABFE is largely confined to **short peptides (roughly ≤5-8
residues)** and usually requires enhanced sampling - Hamiltonian replica
exchange, or similar - rather than plain windows. For 10-20+ residues it is not
the standard tool.

## What to use instead, by ligand class

| Ligand | Method | Why |
|---|---|---|
| ≤8 res, single well-defined pose | **ABFE** | in its domain; DAMGO is the validation case |
| Congeneric series | **RBFE (ΔΔG)** | common scaffold cancels; converges far better |
| 10-20+ res, flexible | **PMF / umbrella sampling** | never annihilates the ligand, so no decoupled-state pathology |
| Everything, as triage | **MM/GBSA** | cheap ranking within a series; poor absolute accuracy |

Two of these deserve emphasis.

**RBFE is the strongest use of the data already in hand.** The benchmark set is
a congeneric series - YPFP / YPFPG / YPFPGPI and YPFVE / YPFVEPI differ by
appended residues. Computing ΔΔG *between* them is far better conditioned than
computing ΔG for each independently, because the shared scaffold's
contributions cancel instead of having to be sampled and then subtracted.

**PMF is the scalable path to long peptides.** Pulling the peptide out of the
pocket along a physical coordinate and integrating the mean force gives
ΔG_bind with a standard-state correction, without ever creating a gas-phase
chain. It handles size and flexibility far better than alchemical decoupling.

## Decision

- **Build ABFE, but scope it as method validation on DAMGO**, not as the
  screening engine. DAMGO is a modified pentapeptide with an experimental
  pose - about as favourable as peptide ABFE gets. If it cannot reproduce a
  known ΔG there, nothing further is worth running.
- **Do not plan to scale ABFE to 10-20mers.** Add PMF/umbrella sampling as the
  route for those.
- **Build Stage 6 so the equilibrated systems feed either engine.** The
  expensive, reusable asset is the equilibrated membrane system, not the
  alchemical machinery on top of it.

For a grant, selecting the method by ligand class is a **stronger** claim than
applying ABFE uniformly. Uniform ABFE on a 20-mer invites the reviewer to ask
whether the authors know its domain of validity; an explicit tiering says they
do.

## Revisit trigger

If the co-folded poses for a peptide converge tightly across seeds AND it is
≤8 residues, it becomes an ABFE candidate. Otherwise it goes to PMF. Seed
divergence, already computed in Stage 2, is the gate.
