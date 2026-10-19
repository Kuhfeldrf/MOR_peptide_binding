# All five met-enkephalin arms in one common reference frame.
#
# 6 frames, 5 ns apart, from a 100 ns run cancelled at 27.6 ns. A flipbook, not
# a movie. See docs/metenk_pose_stability.md.
#
# Receptor is resi 1-281, peptide is resi 898-902. auth = packed + 64 (6DDF).
# (No semicolons in comments: PyMOL splits on ';' even inside a comment line and
# tries to execute what follows.)
#
#   pymol md_viz/metenk/view_all5.pml

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

python
arms = [("exp",   "metenk_exp",  "white"),
        ("ptx",   "metenk_ptx",  "marine"),
        ("hdk",   "metenk_hdk",  "orange"),
        ("chai",  "metenk_chai", "purple"),
        ("scr",   "metenk_scr",  "firebrick")]
for name, stem, colour in arms:
    cmd.load(f"{stem}_frame0.pdb", name)
    cmd.load_traj(f"{stem}_prod.xtc", name, state=1)
    cmd.hide("everything", name)
    cmd.show("sticks", f"{name} and resi 898-902")
    cmd.color(colour, f"{name} and resi 898-902 and elem C")
python end

# Receptor context once, from the experimental arm.
show cartoon, exp and resi 1-281
color grey80, exp and resi 1-281
set cartoon_transparency, 0.6, exp and resi 1-281

select pocket, exp and resi 83+229+233+262
show sticks, pocket
color palegreen, pocket and elem C
set stick_radius, 0.12, pocket
set stick_radius, 0.18, resi 898-902

orient pocket
zoom resi 898-902, 4
bg_color white
set ray_shadows, 0
set movie_fps, 2

print ""
print "  white=exp  blue=ptx  orange=hdk  purple=chai  red=scr"
print "  mplay animates all five together; 'set all_states, on' overlays."
print ""
