# 50 ns MOR + DAMGO production run (Reduced group: receptor, DAMGO, proximal lipids)
load prod_frame0.pdb, md
load_traj prod_viz.xtc, md, state=1
intra_fit md and polymer.protein and name CA
hide everything
show cartoon, polymer.protein
spectrum count, rainbow, md and polymer.protein and name CA
show sticks, resn DAM
color magenta, resn DAM and elem C
show lines, resn PA+PC+OL+CHL
color grey70, resn PA+PC+OL+CHL and elem C
orient resn DAM or polymer.protein
set movie_fps, 20
mplay
