# Protenix pose vs scrambled control, met-enkephalin at MOR.
#
# WHAT YOU ARE WATCHING: 6 frames, 5 ns apart, from a 100 ns production run that
# was cancelled at 27.6 ns by the wall clock. This is a flipbook, not a movie -
# anything happening faster than 5 ns is invisible here. See
# docs/metenk_pose_stability.md.
#
# Every arm is superposed on the SAME reference (the experimental arm's frame-0
# receptor CA), so the two peptides sit in a common frame and can be compared
# directly. Peptide motion is receptor-relative by construction.
#
#   pymol md_viz/metenk/view_ptx_vs_scr.pml     (from anywhere)

# Run it from anywhere: cd to the script's own directory first, so the relative
# data paths below resolve. pymol.__script__ is the path of the running .pml.
python
import os, pymol
try:
    cmd.cd(os.path.dirname(os.path.abspath(pymol.__script__)))
except Exception:
    print("NOTE: could not locate the script directory; "
          "cd to md_viz/metenk and re-run if the loads below fail.")
python end

load metenk_ptx_frame0.pdb, ptx
load_traj metenk_ptx_prod.xtc, ptx, state=1
load metenk_scr_frame0.pdb, scr
load_traj metenk_scr_prod.xtc, scr, state=1

hide everything

# Receptor is resi 1-281, peptide is resi 898-902. Both are protein and both
# land in chain X through the PDB writer, so resi is the only safe handle.
select rec, ptx and resi 1-281
show cartoon, rec
color grey80, rec
set cartoon_transparency, 0.55, rec

select pep_ptx, ptx and resi 898-902
select pep_scr, scr and resi 898-902
show sticks, pep_ptx or pep_scr
color marine, pep_ptx and elem C
color firebrick, pep_scr and elem C
util.cnc("pep_ptx or pep_scr")
set stick_radius, 0.18

# Orthosteric pocket landmarks. The receptor was built from 6DDF and packmol
# renumbered it from 1, so auth = packed + 64 and the Asp of the salt bridge is
# D147 (D3.32), not the 8F7Q D149. Verified against the coordinates:
#   D147 -> 83, W293 -> 229, H297 -> 233, Y326 -> 262
select pocket, ptx and resi 83+229+233+262
show sticks, pocket
color palegreen, pocket and elem C
set stick_radius, 0.12, pocket

orient pep_ptx or pocket
zoom pep_ptx or pep_scr, 4

set movie_fps, 2
set all_states, off
bg_color white
set ray_shadows, 0

print ""
print "  blue   = Protenix pose (metenk_ptx)"
print "  red    = scrambled control (metenk_scr)"
print "  green  = D147 / W293 / H297 / Y326 (packed 83/229/233/262)"
print ""
print "  6 states = 0, 5, 10, 15, 20, 25 ns.  mplay to animate, mstop to halt."
print "  'set all_states, on' overlays all six at once - more honest than the"
print "  animation at this frame spacing."
print ""

mplay
