#!/usr/bin/env python
# coding: utf-8

import numpy as np

import matplotlib.pyplot as plt


x_L = 10.
z_L = 10.   # think of these lengths as domain size in KM

x_loc = 10.0
z_loc = 3.0  # think of these as localization half-widths for radar

ob_loc = [0.5*x_L, 3.]

npx, npz = 51, 51

#--------------------------------------------------------------------------
#
def compute_localization(coords, ob_loc=0.5, vert_norm=0.2, cutoff=1.0):

    def gc_lt1(r):
        return (((-0.25*r +0.5)*r + 0.625)*r - (5.0/3.0))*r**2 + 1.0

    def gc_lt2(r):
        return (((((r/12.0) - 0.5)*r +0.625)*r + (5.0/3.0))*r -5.0)*r + 4.0 - 2.0 / (3.0 * r)


    print('Location of Ob:  {0}'.format(ob_loc))
    location_diff = (coords - ob_loc)/vert_norm

    local = np.zeros((location_diff.shape))

    location_diff  = np.abs( location_diff )

    r = location_diff / cutoff

    local = np.where(r < 2.0*cutoff, gc_lt2(r), local)

    local = np.where(r <     cutoff, gc_lt1(r), local)

    return local
    
#--------------------------------------------------------------------------
#
def compute_localization2(c_x, c_z, ob_loc=[0.5,0.5], vert_norm=[0.2,0.2], cutoff=1.0):

    def gc_lt1(r):
        return (((-0.25*r +0.5)*r + 0.625)*r - (5.0/3.0))*r**2 + 1.0

    def gc_lt2(r):
        return (((((r/12.0) - 0.5)*r +0.625)*r + (5.0/3.0))*r -5.0)*r + 4.0 - 2.0 / (3.0 * r)

    print('Location of Ob:  {0}'.format(ob_loc))
    
    location_x = (c_x - ob_loc[0])/vert_norm[0]

    location_z = (c_z - ob_loc[1])/vert_norm[1]

    location_diff  = np.sqrt(location_x**2 + location_z**2)

    local = np.zeros((location_diff.shape))

    r = location_diff / cutoff

    local = np.where(r < 2.0*cutoff, gc_lt2(r), local)

    local = np.where(r <     cutoff, gc_lt1(r), local)

    return local



#--------------------------------------------------------------------------
#
# Main program

xg = np.linspace(0.0, x_L, npx)
zg = np.linspace(0.0, z_L, npz)

xx, zz = np.meshgrid(xg, zg)

fig, ax = plt.subplots(1,3, figsize=(12,6))


locx = compute_localization(xx, ob_loc=ob_loc[0], vert_norm=9., cutoff=1.0);
locz = compute_localization(zz, ob_loc=ob_loc[1], vert_norm=3., cutoff=1.0);
gsi_loc = locx * locz

ax[0].contour(xx, zz, gsi_loc, 10, colors='k')
ax[0].set_title("GSI Localization")
ax[0].set_xlabel('X (km)')


dart_loc = compute_localization2(xx, zz, ob_loc=ob_loc, vert_norm=[9.,3.], cutoff=1.0)
ax[1].contour(xx, zz, dart_loc, 10, colors='k')
ax[1].set_title("Dart Localization")
ax[1].set_xlabel('X (km)')


diff = dart_loc-gsi_loc
ax[2].contour(xx, zz, diff, 6, colors='k')
ax[2].set_title("Difference between \n Dart-GSI Localization\n Max: {:.2e}".format(diff.max()))
ax[2].set_xlabel('X (km)')

plt.savefig("Hloc.png")



