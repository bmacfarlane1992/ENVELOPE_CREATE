'''

gridmake.py

Routine to take SPH distribution of particles, and construct grid based on
particle location limits for world size. Method grids based on main.py variable
grid_method, for either regular, or octree grid (dev.).

Code outputs to amr_grid.inp as per required RADMC-3D IO, and returns grid
to main.py to build density grid.

Author: Benjamin MacFarlane
Date: 19/07/2017
Contact: bmacfarlane@uclan.ac.uk

'''
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
nbins = 64     # For 'reg' grid_make selection, how many bins in each dimension
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import time
import sys
from hyperion.importers.sph import construct_octree
#
# To open the hyperion from the source, open a shell then inputs the following:
#
# >>> import os
# >>> import hyperion.importers.sph
# >>> os.system('atom '+hyperion.importers.sph)
#
try:
    import numpy as np
except ImportError:
    np = None
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def reg(arch_dir, pos, r_dust):
#
    # Define world size and bin size based on nbins variable
#
    binlims = [-1.005*np.amax(pos), 1.005*np.amax(pos)]
    bin_it = (binlims[1] - binlims[0]) / float(nbins)
#
    # Write to amr_grid.inp for RADMC-3D the locations of bin edges in all
    # dimensions. As grid is cubic, r_it = x_it = y_it = z_it
#
    f = open(arch_dir+'amr_grid.inp','w')
    f.write('1 \n')
    f.write('0 \n')
    f.write('0 \n')
    f.write('0 \n')
    f.write('1 1 1 \n')
    f.write(str(nbins)+' '+str(nbins)+' '+str(nbins)+'\n')
    for i in range(0,nbins+1):
        x = ( binlims[0] + (i * bin_it) )
        f.write(str(x)+'\n')
    for i in range(0,nbins+1):
        y = ( binlims[0] + (i*bin_it) )
        f.write(str(y)+'\n')
    for i in range(0,nbins+1):
        z = ( binlims[0] + (i*bin_it) )
        f.write(str(z)+'\n')
    f.close()
    count = 0
#
    # Now compute the radial values for each [x,y,z] bin, to eventually feed
    # into profiles.py
#
    rgrid = []
    for zz in range(0, nbins):
        for yy in range(0, nbins):
            for xx in range(0, nbins):
                x = binlims[0] + ( xx * bin_int ) + ( 0.5 * bin_int )
                y = binlims[0] + ( yy * bin_int ) + ( 0.5 * bin_int )
                z = binlims[0] + ( zz * bin_int ) + ( 0.5 * bin_int )
                r[count] = np.sqrt(x**2. + y**2. + z**2.) / r_dust
                count = count + 1
    rgrid = np.array(rgrid)
#
    # Set up a 3D float array for all bins in cubic grid to be used in
    # densgrid.py
#
    grid =  np.array(  \
       [ [ [[0.] for i in range(nbins)] for j in range(nbins) ] for k in range(nbins) ] \
       )
#
    return nbins, bin_it, binlims, grid, rgrid




### ------------------------------------------------------------------------ ###
### ------------------------------------------------------------------------ ###
### ------------------------------------------------------------------------ ###
### ------------------------------------------------------------------------ ###
### ------------------------------------------------------------------------ ###




def oct(arch_dir, pos, rho, h, m_part, r_dust):

### ------------------------------------------------------------------------ ###
    ### Set up octree global variables, modifying envlope.py outputs as necessary  ###
### ------------------------------------------------------------------------ ###

    plims = [-np.amax(pos), np.amax(pos)]
    WORLD_SIZE = 1.1*(plims[1]-plims[0])
    dx = 0.5*WORLD_SIZE ; dy = dx ; dz = dx
    n_octbase = 2 ; octbase = dx ; base0 = 0. - dx
    m_part = np.array( [m_part]*len(pos[0]) )

### ------------------------------------------------------------------------ ###
    ### Run Hyperion octree code, and loop over outputs to find nleafsmax and nbranchmax  ###
### ------------------------------------------------------------------------ ###

    levelmax = 15
#
    def stop(x, y, z, dx, dy, dz, px, py, pz, sigma):
        return len(px) <= 10
#
    o = construct_octree(0., 0., 0., dx, dy, dz, \
       np.array(pos[0]), np.array(pos[1]), np.array(pos[2]), np.array(h), m_part, \
       n_levels = levelmax, stopping_criterion=stop)
#
    n = len(o.refined)
    nleafsmax = 0 ; nbranchmax = 0
    for i in range(0, n):
        if o.refined[i] == True:
            nbranchmax += 1
        elif o.refined[i] == False:
            nleafsmax += 1
    ncells = nleafsmax
#
#    levelmax = raw_input("Input the maximum recursion level that has been reached: ")
#
### ------------------------------------------------------------------------ ###
    ### Use results to ouput grid data into amr_grid.inp for RADMC-3D  ###
### ------------------------------------------------------------------------ ###

    f = open(arch_dir+'amr_grid.inp','w')
    f.write('1 \n')
    f.write('1 \n')
    f.write('0 \n')
    f.write('0 \n')
    f.write('1 1 1 \n')
    f.write(str(1)+' '+str(1)+' '+str(1)+'\n')
    f.write(str(levelmax)+' '+str(n)+' '+str(n)+'\n')
    f.write(str(base0)+' '+str(-base0)+'\n')
    f.write(str(base0)+' '+str(-base0)+'\n')
    f.write(str(base0)+' '+str(-base0)+'\n')
    for i in range(0, n):
        if o.refined[i] == True:
            f.write('1\n')
        else:
            f.write('0\n')
    f.close()
#
### ------------------------------------------------------------------------ ###
    ### Use results to ouput grid to indicate radial locations of leaf nodes  ###
### ------------------------------------------------------------------------ ###
#
    rgrid = []
    for i in range(0, n):
        if o.refined[i] == False:
            rgrid.append( np.sqrt( (o['xcen'][0].array[i])**2. + \
               (o['ycen'][0].array[i])**2. + (o['zcen'][0].array[i])**2. ) / r_dust )
        else:
            continue
    rgrid = np.array(rgrid)
#
    return o, n, ncells, rgrid
