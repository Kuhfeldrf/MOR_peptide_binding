# Headless render: the six surviving production frames of the Protenix pose and
# the scrambled control, overlaid in the pocket.
#
# An overlay, not an animation: at 5 ns spacing an animation implies a
# continuity the data does not have.
#
# Receptor context is loaded WITHOUT its trajectory, so it stays a single state
# while the peptides carry six - six ray-traced transparent cartoons is minutes
# of render for no information.
load metenk_exp_frame0.pdb, context
load metenk_ptx_frame0.pdb, ptx
load_traj metenk_ptx_prod.xtc, ptx, state=1
load metenk_scr_frame0.pdb, scr
load_traj metenk_scr_prod.xtc, scr, state=1

hide everything
remove hydrogens          # 6 overlaid states of explicit H is unreadable
show cartoon, context and resi 1-281
color grey90, context and resi 1-281
set cartoon_transparency, 0.45, context   # cheap: context is a single state

show sticks, (ptx or scr) and resi 898-902
color marine, ptx and resi 898-902 and elem C
color firebrick, scr and resi 898-902 and elem C
show sticks, context and resi 898-902
color grey50, context and resi 898-902 and elem C
set stick_radius, 0.13

select pocket, context and resi 83+229+233+262
show sticks, pocket
color palegreen, pocket and elem C
set stick_radius, 0.10, pocket

set all_states, on
bg_color white
set ray_shadows, 0
set ray_opaque_background, 1
set antialias, 1
orient (ptx or scr) and resi 898-902
zoom (ptx or scr) and resi 898-902, 7
png overlay_view1.png, width=1100, height=850, ray=1
turn y, 90
png overlay_view2.png, width=1100, height=850, ray=1
