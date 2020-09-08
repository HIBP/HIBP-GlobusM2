import numpy as np
import hibplib as hb
import hibpplotlib as hbplot
import copy

# %%
''' pass trajectories to different slits
'''
Ebeam = 40.
UA2 = 0.0

n_slits = 5
# add slits to Geometry
geomGlob.add_slits(n_slits=n_slits, slit_dist=0.01, slit_w=5e-3,
                   slit_l=0.1, slit_gamma=0.)
r_slits = geomGlob.slits_edges
rs = geomGlob.r_dict['slit']
# calculate normal to slit plane
slit_plane_n = geomGlob.slit_plane_n

# %%
traj_list_copy = copy.deepcopy(traj_list_oct)
# traj_list_copy = copy.deepcopy(traj_list_passed)

# %%
print('*** Passing fan to {} slits'.format(n_slits))
for tr in traj_list_copy:
    if tr.Ebeam == Ebeam and tr.U[0] == UA2:
        print('\nEb = {}, UA2 = {}'.format(tr.Ebeam, tr.U[0]))
    else:
        continue

    tr = hb.pass_to_slits(tr, dt, E, B, geomGlob, timestep_divider=15)
    break

# %% plot trajectories
hbplot.plot_traj_toslits(tr, geomGlob, Btor, Ipl,
                         plot_fan=True, plot_flux=False)
