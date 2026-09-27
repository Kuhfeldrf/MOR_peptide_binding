#!/usr/bin/env python3
"""Shared ligand-identity helpers.

WHY THIS MODULE EXISTS
----------------------
The same error has now occurred four times: code that needs the ligand was
handed something that was the wrong molecule, in the wrong frame, or under the
wrong atom names.

  1. ligand_params received the whole receptor+peptide complex. Caught by a
     charge assertion.
  2. membrane_complex received the whole complex as its --ligand, duplicating
     the receptor into a 2294-atom "LIG" (D35).
  3. Stage 5 MM/GBSA: the LIG in complex_opm.pdb has the predictor's protein
     atom names, repeated across residues, while the mol2 has 75 uniquely named
     atoms - tleap declared 78 atoms missing and built them from template.
  4. membrane_topology: the same naming mismatch, on the packed system.

Each was fixed at its own call site, and the next site repeated it. The fix
belongs in one place that every caller goes through.

THE INVARIANT
-------------
Parameters live in the mol2. Coordinates live wherever the pipeline last put
the ligand. Those two carry DIFFERENT atom names and often different atom
counts, because parameterisation completes the chemistry (PDBFixer adds the
C-terminal OXT and all hydrogens) while the pipeline's copy is whatever the
structure predictor emitted.

So the ligand handed to tleap must always be built by taking the
mol2-consistent file and MOVING it onto the coordinates in question - never by
renaming or patching the coordinates themselves.
"""
from __future__ import annotations

import pathlib

import numpy as np


def pdb_atoms(path: pathlib.Path) -> list[str]:
    return [l for l in path.read_text().splitlines()
            if l.startswith(("ATOM", "HETATM"))]


def coords(lines: list[str]) -> np.ndarray:
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in lines])


def is_hydrogen(line: str) -> bool:
    el = (line[76:78].strip() or line[12:16].strip()[0]).upper()
    return el == "H"


def kabsch(P: np.ndarray, Q: np.ndarray):
    """Rigid transform taking P onto Q. Rotation only, no scaling."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def ligand_in_frame(param_pdb: pathlib.Path, target_lines: list[str],
                    resname: str = "LIG", chain: str = "L", resid: int = 1,
                    tol: float = 0.1):
    """Return the mol2-consistent ligand moved onto `target_lines`' frame.

    `param_pdb` is params/ligand_unique.pdb - the file whose atom names and
    count match the mol2. `target_lines` are the ligand's ATOM records wherever
    the pipeline currently has it (a prepared complex, a packed membrane).

    Both derive from the same predicted ligand and preserve heavy-atom order,
    so the transform between them is rigid and exact. It is VERIFIED rather
    than assumed: a silent mis-superposition would place a correctly named
    ligand in the wrong position, which nothing downstream would flag.

    Returns (lines, rmsd, n_extra).
    """
    param = pdb_atoms(param_pdb)
    heavy = [l for l in param if not is_hydrogen(l)]
    n = len(target_lines)
    if n == 0:
        raise ValueError("no target ligand atoms given")
    if len(heavy) < n:
        raise ValueError(
            f"{param_pdb} has {len(heavy)} heavy atoms, fewer than the {n} "
            f"in the target - these are not the same molecule")

    R, t = kabsch(coords(heavy[:n]), coords(target_lines))
    fit = float(np.sqrt(np.mean(np.sum(
        ((R @ coords(heavy[:n]).T).T + t - coords(target_lines)) ** 2, axis=1))))
    if fit > tol:
        raise ValueError(
            f"superposition RMSD {fit:.3f} A exceeds {tol} A - these should be "
            f"the same conformation in two frames, so the atom order does not "
            f"correspond")

    out = []
    for l in param:
        v = R @ np.array([float(l[30:38]), float(l[38:46]),
                          float(l[46:54])]) + t
        out.append(l[:17] + f"{resname:>3} {chain}{resid:>4}" + l[26:30]
                   + f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}" + l[54:])
    return out, fit, len(heavy) - n
