'''

densgrid.py

Routine to take SPH distribution of particles, populating each grid bin with
an associated density.
Code ouputs grid results to density.inp for required RADMC-3D IO. Also returns
rhogrid to main.py for anaylysis of profiles.

Author: Benjamin MacFarlane
Date: 12/05/2017
Contact: bmacfarlane@uclan.ac.uk

'''
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
find_kern = True
gamma_correct = False
#
dim = 3         # Considering a 3D model - variable for binary search exit condition
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import numpy as np
import math
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def particle(arch_dir, pos, rho, m_part, h, grid, bin_it, limits, nbins):

    print "\nUsing finders-keepers method to compute kernel contrbutions to bin location\n"

    rhogrid = np.array(grid)

    if (math.log(nbins,2) % 1 != 0):
        print "\nRevise nbins, to be used in binary search"
        print "(Must be n, where nbins = 2^{n})\n"
        return

    for part in range(0, len(pos[0])):

        if ( float(part+1) % ( len(pos[0]) / 10 ) == 0):
            print "Particle "+str(part+1)+" of "+str(len(pos[0]))+" being processed"

        ploc = [0.,0.,0.] ; prho = 0.

        ploc[0] = pos[0][part] ; ploc[1] = pos[1][part] ; ploc[2] = pos[2][part]
        prho = rho[part]

        dimpoint = 0
        ibin = [0,0,0]
        pbin = [[0,0],[0,0],[0,0]]
        pc = [0,0,0]

        while dimpoint < dim:

            levels = math.log(nbins,2)
            bottom = limits[0] ; top = limits[1]
            refine = 0

            while refine < levels:
                middle = (bottom + top) / 2.
                if ploc[dimpoint] < middle and ploc[dimpoint] > bottom:
                    pbin[dimpoint][0] = bottom ; pbin[dimpoint][1] = middle
                    bottom = bottom ; top = middle
                    refine = refine + 1
                elif ploc[dimpoint] > middle and ploc[dimpoint] < top:
                    pbin[dimpoint][0] = middle ; pbin[dimpoint][1] = top
                    bottom = middle ; top = top
                    refine = refine + 1
                else:
                    "Particle outside of binning domain"
                    return
            pc[dimpoint] = (pbin[dimpoint][0]+pbin[dimpoint][1]) / 2.
            ibin[dimpoint] = int((pc[dimpoint] + 0.5*bin_it - limits[0] ) / bin_it ) - 1
            dimpoint = dimpoint + 1
#
    # Find kernel contribution to centre of bin, to sum densities in a consistent manner
#
        if find_kern is True:
#
            r_sph = np.sqrt( \
            (pos[2][part] - pc[2])**2. + (pos[1][part] - pc[1])**2. + (pos[0][part] - pc[0])**2.
               )
            q = r_sph / h[part]
#
            w = 0.
            if ((q < 1.) and (q > 0.)):
                w = (1./(4.*math.pi)) * ( (2. - q)**3. - 4*(1 - q)**3. )
            elif ((q < 2.) and (q > 1.)):
                w = (1./(4.*math.pi)) * ( (2. - q)**3. )
            else:
                w = 0.
#
            w = w / h[part]**(3.)
#
    # Apply density correction based on the ratio of bin and particle volumes
#
            if gamma_correct is True:
#
                volbin = deltp**(3.)
                volsph = (4./3.) * math.pi * (2.*h[part])**(3.)
                gamma = volsph / volbin
#
                if (gamma < 1.):
                    prho = prho * gamma
#
#
#            rhogrid[ibin[0]][ibin[1]][ibin[2]] = \
#               (rhogrid[ibin[0]][ibin[1]][ibin[2]] + (w * prho) )
            rhogrid[ibin[0]][ibin[1]][ibin[2]] = \
               (rhogrid[ibin[0]][ibin[1]][ibin[2]] + (w * m_part) )
#
        elif find_kern is False:
#
            rhogrid[ibin[0]][ibin[1]][ibin[2]] = \
               (rhogrid[ibin[0]][ibin[1]][ibin[2]] + prho )
#
            if gamma_correct is True:
#
                volbin = deltp**(3.)
                volsph = (4./3.) * math.pi * (2.*h[part])**(3.)
                gamma = volsph / volbin
#
                if (gamma < 1.):
                    prho = prho * gamma
#
#
    # Write data to density.inp for RADMC-3D IO
#
    print "\nWriting to dust_density.inp\n"
    ngrid = int(nbins**(3.))
    f = open(arch_dir+'dust_density.inp','w')
    f.write('1 \n')
    f.write(str(ngrid)+'\n')
    f.write('1 \n')
    for i in range(nbins):
        for j in range(nbins):
            for k in range(nbins):
                f.write(str(rhogrid[i][j][k][0])+'\n')
    f.close()
#
    return rhogrid

### ------------------------------------------------------------------------ ###

def kernel(arch_dir, pos, rho, m_part, h, grid, bin_it, limits, nbins):
#
    # First, compute number of surrounding bins contributing to bin under
    # inspection, adding 1 to account for limiting contributions
#
    sph_inf = int(  ( 2. * max(h) ) / bin_it ) + 1
#
    print "\nUsing sphere-of-influence method to compute "
    print "kernel contrbutions to bin location,"
    print "with ", str(sph_inf), " surrounding bin(s)"
#
    # Initialise arrays to be filled, and do sanity checks on the grid adopted
    # (blank grid used to account for any out-of-domain particle finds)
#
    blank_grid = [ [ [[] for i in range(nbins+sph_inf)] \
       for j in range(nbins+sph_inf) ] for k in range(nbins+sph_inf) ]
    pgrid = blank_grid ; rhogrid = grid
#
    if (math.log(nbins,2) % 1 != 0):
        print "\nRevise nbins, to be used in binary search"
        print "(Must be n, where nbins = 2^{n})\n"
        return
    print "\n"
#
    # Now use bianry search algoithm to identify the particles belonging
    # in each bin, saving index of particles for kernel contribution work
#
    for part in range(0, len(pos[0])):
#
        if ( float(part+1) % ( float(len(pos[0])) / float(10) ) == 0):
            print "Particle "+str(part+1)+" of "+str(len(pos[0]))+" being binned"
#
        ploc = [0.,0.,0.]
        ploc[0] = pos[0][part] ; ploc[1] = pos[1][part] ; ploc[2] = pos[2][part]
#
        dimpoint = 0
        pc = [0,0,0] ; pbin = [[0,0],[0,0],[0,0]] ; ibin = [0,0,0]
#
        while dimpoint < dim:
            levels = math.log(nbins,2)
            bottom = limits[0] ; top = limits[1]
            refine = 0
            while refine < levels:
                middle = (bottom + top) / 2.
                if ploc[dimpoint] < middle and ploc[dimpoint] > bottom:
                    pbin[dimpoint][0] = bottom ; pbin[dimpoint][1] = middle
                    bottom = bottom ; top = middle
                    refine = refine + 1
                elif ploc[dimpoint] > middle and ploc[dimpoint] < top:
                    pbin[dimpoint][0] = middle ; pbin[dimpoint][1] = top
                    bottom = middle ; top = top
                    refine = refine + 1
                else:
                    "Particle outside of binning domain"
                    return
#
            pc[dimpoint] = (pbin[dimpoint][0]+pbin[dimpoint][1]) / 2.
            ibin[dimpoint] = int((pc[dimpoint] + 0.5*bin_it - limits[0] ) \
                / bin_it ) - 1 + sph_inf
            dimpoint = dimpoint + 1
#
        pgrid[ibin[0]][ibin[1]][ibin[2]].append(part)
#
    # With grid filled with particle indices and using the sphere of influence,
    # loop over particles in bins in sphere of influence (including central bin)
    # and add kernel contribution to the density a the bin centre
#
    print "Finding kernel contributions for base bin"
#
    count = [0]*(nbins**3) ; gridn = 0
#
    for igrid in range(0, nbins):                               # main grid - z array
#
        if (float(igrid+1) / math.log(nbins,2) % 1. == 0.):
            print "\t"+str(igrid+1)+" of "+str(nbins)
#
        for jgrid in range(0, nbins):                           # ... y array
            for kgrid in range(0, nbins):                       # ... x array
#
                pc[2] = limits[0] + (igrid * bin_it) + ( 0.5 * bin_it )
                pc[1] = limits[0] + (jgrid * bin_it) + ( 0.5 * bin_it )
                pc[0] = limits[0] + (kgrid * bin_it) + ( 0.5 * bin_it )
#
                ibin_inf = [0,0,0]
#
                for iinf in range(-sph_inf, sph_inf):                # Influence grid -  z array
                    for jinf in range(-sph_inf, sph_inf):            # ... y array
                        for kinf in range(-sph_inf,sph_inf):         # ... x array
#
                            ibin_inf[2] = igrid + sph_inf + iinf
                            ibin_inf[1] = jgrid + sph_inf + jinf
                            ibin_inf[0] = kgrid + sph_inf + kinf
#
    # If influencing bin locatin is outside of grid domain, continue
#
                            if (ibin_inf[2] < 0) or (ibin_inf[2] >= nbins):
                                continue
                            if (ibin_inf[1] < 0) or (ibin_inf[1] >= nbins):
                                continue
                            if (ibin_inf[0] < 0) or (ibin_inf[0] >= nbins):
                                continue
#
    # Count to find number of particle contributions to each bin
#
                            count[gridn] = count[gridn] + len(pgrid[ibin_inf[0]][ibin_inf[1]][ibin_inf[2]][:])
#
                            for part in range(0, len(pgrid[ibin_inf[0]][ibin_inf[1]][ibin_inf[2]] ) ):
#
                                r_sph = np.sqrt( \
                                (pos[2][pgrid[ibin_inf[0]][ibin_inf[1]][ibin_inf[2]][part]] - pc[2])**2. + \
                                (pos[1][pgrid[ibin_inf[0]][ibin_inf[1]][ibin_inf[2]][part]] - pc[1])**2. + \
                                (pos[0][pgrid[ibin_inf[0]][ibin_inf[1]][ibin_inf[2]][part]] - pc[0])**2.
                                )
                                q = r_sph / h[pgrid[ibin_inf[0]][ibin_inf[1]][ibin_inf[2]][part]]
#
                                w = 0.
                                if ((q < 1.) and (q > 0.)):
                                    w = (1./(4.*math.pi)) * ( (2. - q)**3. - 4*(1 - q)**3. )
                                elif ((q < 2.) and (q > 1.)):
                                    w = (1./(4.*math.pi)) * ( (2. - q)**3. )
                                else:
                                    w = 0.
#
                                w = w / h[pgrid[ibin_inf[0]][ibin_inf[1]][ibin_inf[2]][part]]**(3.)
#
                                rhogrid[ibin_inf[0]][ibin_inf[1]][ibin_inf[2]] = \
                                ( rhogrid[ibin_inf[0]][ibin_inf[1]][ibin_inf[2]] + \
                                (w*m_part)

                                )
#
                gridn = gridn+1
    print "Average number of particles contributing to bin density is:", np.mean(count)
#
#
    # Write data to density.inp for RADMC-3D IO
#
    print "\nWriting to dust_density.inp\n"
    ngrid = int(nbins**(3.))
    f = open(arch_dir+'dust_density.inp','w')
    f.write('1 \n')
    f.write(str(ngrid)+'\n')
    f.write('1 \n')
    for i in range(nbins):
        for j in range(nbins):
            for k in range(nbins):
                f.write(str(rhogrid[i][j][k][0])+'\n')
    f.close()
#
    return rhogrid

### ------------------------------------------------------------------------ ###

def brute(arch_dir, pos, rho, m_part, h, grid, bin_it, limits, nbins):

    print "\nUsing brute force method to compute kernel contrbutions to bin location\n"
#
    rhogrid = grid
#
    for ii in range(0, nbins):
        z0 = limits[0] + (ii * bin_it) + ( 0.5 * bin_it )

        if (float(ii+1) / 4. % 1. == 0.):
            print "On bin "+str(ii+1)+" of "+str(nbins)

        for jj in range(0, nbins):
            y0 = limits[0] + (jj * bin_it) + ( 0.5 * bin_it )

            for kk in range(0, nbins):
                x0 = limits[0] + (kk * bin_it) + ( 0.5 * bin_it )
#
                rho_contrib = 0. ; bin_rho = 0.
#
                for part in range (0, len(pos[0]) ):

                    w = 0.
                    r =  np.sqrt((x0 - pos[0][part])**2. + (y0 - pos[1][part])**2. + (z0-pos[2][part])**2.)
                    q = r / h[part]
#
                    if ((q < 1.) and (q > 0.)):
                        w = (1./(4.*math.pi)) * ( (2. - q)**3. - 4*(1 - q)**3. )
                    if ((q < 2.) and (q > 1.)):
                        w = (1./(4.*math.pi)) * ( (2. - q)**3. )
                    else:
                        q = 0.
#
                    w = w / h[part]**(3.)
#
                    rho_contrib = (w * m_part)
                    rhogrid[ii][jj][kk] = rhogrid[ii][jj][kk] + rho_contrib
#
#
    # Write data to density.inp for RADMC-3D IO
#
    print "\nWriting to dust_density.inp\n"
    ngrid = int(nbins**(3.))
    f = open(arch_dir+'dust_density.inp','w')
    f.write('1 \n')
    f.write(str(ngrid)+'\n')
    f.write('1 \n')
    for i in range(nbins):
        for j in range(nbins):
            for k in range(nbins):
                f.write(str(rhogrid[i][j][k][0])+'\n')
    f.close()
#
    return rhogrid

### ------------------------------------------------------------------------ ###

def byhand(arch_dir, grid, bin_it, limits, nbins, rho_0, r_dust, p):

    print "\nSetting all bin density values by hand\n"

    rhogrid = grid

    if (p == 0):
#
        for ii in range(0, nbins):
            for jj in range(0, nbins):
                for kk in range(0, nbins):
                    rhogrid[ii][jj][kk] = rho_0
#
#
    elif (p == 2):
#
        for ii in range(0, nbins):
            x = limits[0] + ( ii * bin_it ) + ( 0.5 * bin_it )
#
            for jj in range(0, nbins):
                y = limits[0] + ( jj * bin_it ) + ( 0.5 * bin_it )
#
                for kk in range(0, nbins):
                    z = limits[0] + ( kk * bin_it ) + ( 0.5 * bin_it )
                    r = np.sqrt(x**2. + y**2. + z**2.) / r_dust
#
                    rhogrid[ii][jj][kk] = rho_0 * r**(-p)
#
#
    # Write data to density.inp for RADMC-3D IO
#
    print "\nWriting to dust_density.inp\n"
    ngrid = int(nbins**(3.))
    f = open(arch_dir+'dust_density.inp','w')
    f.write('1 \n')
    f.write(str(ngrid)+'\n')
    f.write('1 \n')
    for i in range(nbins):
        for j in range(nbins):
            for k in range(nbins):
                f.write(str(rhogrid[i][j][k][0])+'\n')
    f.close()
#
    return rhogrid
