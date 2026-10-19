#!/usr/bin/env python3
"""Small multiples: the three pose-stability readouts for all five arms.

Six points per arm at a 5 ns stride - every point is drawn as a marker, because
drawing these as smooth lines alone would imply sampling the run does not have.

Palette: validated categorical slots 1/2/3/7 for the three predicted poses and
the scramble; the experimental reference wears text ink, not a series hue, since
it is the yardstick rather than one of the things being compared. Both
non-predicted arms are dashed, so identity never rests on color alone.
"""
import csv
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path(__file__).resolve().parents[2]
SRC = REPO / "results/04_md_metenk/pose_stability_per_frame.tsv"
OUT = REPO / "results/04_md_metenk"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#dcdbd6"

# arm -> (label, colour, linestyle, z)
STYLE = {
    "metenk_ptx":  ("Protenix",     "#2a78d6", "-",  5),
    "metenk_chai": ("Chai-1",       "#eb6834", "-",  4),
    "metenk_hdk":  ("HADDOCK3",     "#1baf7a", "-",  4),
    "metenk_scr":  ("scrambled",    "#4a3aa7", "--", 3),
    "metenk_exp":  ("experimental", INK_2,     "--", 3),
}
ORDER = ["metenk_exp", "metenk_ptx", "metenk_chai", "metenk_hdk", "metenk_scr"]

PANELS = [
    ("pep_rmsd_A",    "Peptide RMSD", "Å", None),
    ("com_disp_A",    "Centre-of-mass displacement", "Å", None),
    ("jaccard_vs_t0", "Contact-set overlap with t=0", "Jaccard", (0, 1.05)),
]


def load():
    rows = {}
    with open(SRC) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.setdefault(r["arm"], []).append(
                {k: float(v) for k, v in r.items() if k != "arm"})
    return rows


def main():
    data = load()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), facecolor=SURFACE)
    fig.subplots_adjust(left=0.055, right=0.86, top=0.80, bottom=0.145, wspace=0.28)

    for ax, (key, title, unit, ylim) in zip(axes, PANELS):
        ax.set_facecolor(SURFACE)
        for arm in ORDER:
            label, colour, ls, z = STYLE[arm]
            pts = data[arm]
            x = [p["time_ns"] for p in pts]
            y = [p[key] for p in pts]
            ax.plot(x, y, ls=ls, color=colour, lw=2.0, zorder=z,
                    marker="o", ms=5.0, mew=1.6, mec=SURFACE, clip_on=False)
        ax.set_title(title, color=INK, fontsize=11.5, loc="left", pad=9)
        ax.set_xlabel("production time (ns)", color=INK_2, fontsize=9.5)
        ax.set_ylabel(unit, color=INK_2, fontsize=9.5)
        ax.set_xlim(0, 25)
        ax.set_xticks([0, 5, 10, 15, 20, 25])
        if ylim:
            ax.set_ylim(*ylim)
        else:
            ax.set_ylim(0, None)
        ax.grid(True, color=GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=INK_2, labelsize=9, length=3)

    # Direct labels on the last panel, which is bounded and has room; this is
    # also the relief the contrast WARN on the aqua slot obliges.
    ax = axes[-1]
    last = {a: data[a][-1]["jaccard_vs_t0"] for a in ORDER}
    # Push labels apart where two arms end within 0.045 Jaccard of each other,
    # so the text never overlaps even though the marks legitimately do.
    placed = []
    for arm, y in sorted(last.items(), key=lambda kv: -kv[1]):
        if placed and placed[-1][1] - y < 0.045:
            y = placed[-1][1] - 0.045
        placed.append((arm, y))
    for arm, y in placed:
        label, colour, _, _ = STYLE[arm]
        ax.annotate(label, xy=(25, y), xytext=(6, 0), textcoords="offset points",
                    color=colour if arm != "metenk_exp" else INK_2,
                    fontsize=9.5, va="center", annotation_clip=False)

    fig.suptitle("Met-enkephalin pose stability by starting pose — "
                 "first 25 ns of a run cancelled at 27.6 ns",
                 color=INK, fontsize=13, x=0.055, ha="left", y=0.955)
    fig.text(0.055, 0.885,
             "6 frames per arm at a 5 ns stride · n = 1 trajectory per arm · "
             "receptor-CA superposed, PBC-corrected",
             color=INK_2, fontsize=9.5, ha="left")

    handles = [plt.Line2D([], [], color=STYLE[a][1], ls=STYLE[a][2], lw=2.0,
                          marker="o", ms=5.0, mec=SURFACE, mew=1.6,
                          label=STYLE[a][0]) for a in ORDER]
    axes[0].legend(handles=handles, loc="upper left", fontsize=9,
                   frameon=False, labelcolor=INK_2, handlelength=2.4)

    for ext in ("png", "svg"):
        fig.savefig(OUT / f"metenk_pose_stability.{ext}", dpi=150,
                    facecolor=SURFACE)
    print("wrote", OUT / "metenk_pose_stability.png")


if __name__ == "__main__":
    main()
