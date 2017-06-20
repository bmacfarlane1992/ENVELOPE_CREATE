'''

profiles.py

Routine to take grid outputted density, and RADMC-3D outputs, and plot in
radial distributions, to check both the inputs and outputs to radiative transfer
code.

Author: Benjamin MacFarlane
Date: 04/05/2017
Contact: bmacfarlane@uclan.ac.uk

'''
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
profile_comp = True
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import numpy as np
import math
import matplotlib.pyplot as plt
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def plotrho(arch_dir, bin_int, binlims, nbins, r_sph, r_dust, rho_sph, rho_0, ds_dir, \
   p, tau, plt_dir):
#
    # Set up arrays/variables to feed in data
#
    ngrid = int( float(nbins)**(3.) )
    r = [0.] * ngrid ; rho = [] ; temp = []
#
    f_rho = open(arch_dir+'dust_density.inp')
#
    # Now read in data (note that rho is artificially added to to prevent
    # log-log plot issues later)
#
    for i in range(0, 3):
        trash = f_rho.readline()
    for lines in f_rho:
        lines = lines.strip() ; columns = lines.split()
        rho.append(float(columns[0])+1.e-99)
    f_rho.close()
#
    # Now fill arrays, comuputing locations
#
    count = 0
    for zz in range(0, nbins):
        for yy in range(0, nbins):
            for xx in range(0, nbins):
                x = binlims[0] + ( xx * bin_int ) + ( 0.5 * bin_int )
                y = binlims[0] + ( yy * bin_int ) + ( 0.5 * bin_int )
                z = binlims[0] + ( zz * bin_int ) + ( 0.5 * bin_int )
                r[count] = np.sqrt(x**2. + y**2. + z**2.) / r_dust
                count = count + 1
    r = np.array(r) ; rho = np.array(rho) ; temp = np.array(temp)
#
    # Short procedure to find number of unique radii from r array, using indices
    # of these to produce azimuthally averaged density and temperature profiles
#
    r_plt, indices = np.unique(r, return_inverse=True)
    rho_plt = [0.]*len(r_plt) ; count = [0.]*len(r_plt)

    for i in range(0, len(indices)):
        count[indices[i]] = count[indices[i]] + 1.
        rho_plt[indices[i]] = rho_plt[indices[i]] + rho[i]
    for i in range(0, len(count)):
        if (count[i] != 0):
            rho_plt[i] = rho_plt[i] / count[i]
        else:
            continue
#
    # If called for, read in DS output densities to plot in comparison
#
    if profile_comp is True:
        f_ds = open(ds_dir+'outcells.dat')
        trash = f_ds.readline()
        r_ds = [] ; rho_ds = []
        for lines in f_ds:
            lines = lines.strip() ; columns = lines.split()
            r_ds.append(float(columns[0])/r_dust)
            rho_ds.append(float(columns[7]))
        r_ds = np.array(r_ds) ; rho_ds = np.array(rho_ds)
        f_ds.close()
#
    # Generate the analytic distribution of density, to compare to grid outputs
#
    rho_an = [0.]*len(r_plt)
    for i in range(0, len(r_plt)):
        rho_an[i] = rho_0 * (r_plt[i])**(-p)
#
    # Now plot azimuthally averaged profile
#
    rho_l = np.mean(rho_plt) / 1e5
    plt.scatter(r_plt, rho_plt, label='Grid Output', color = 'b')
    plt.plot(r_plt, rho_an, linestyle = 'solid', color = 'g', label='Analytic')
    plt.legend(loc = 'lower left', fontsize=16, scatterpoints=1)
    plt.xlim(8e-1, 4e3)
    plt.ylim(rho_l,max(rho_plt)*3.)
    plt.yscale('log') ; plt.xscale('log')
    plt.xlabel('Radius (AU)') ; plt.ylabel('Volume density '+r'(g cm$^{-3}$)')
#
    if profile_comp is True:
        plt.plot(r_ds, rho_ds, color = 'r', label = 'DS test')
        for i in range(0, len(r_sph)):
            r_sph[i] = r_sph[i] / r_dust
        plt.scatter(r_sph, rho_sph, color = '0.75', label = 'SPH Output')
        plt.legend(loc = 'lower left', fontsize = 16, scatterpoints = 1)
#
    plt.savefig(plt_dir+'rho_r.png') ; plt.clf()
#
    return


### ------------------------------------------------------------------------ ###


def plottemp(arch_dir, bin_int, binlims, nbins, r_dust, ds_dir, p, tau, plt_dir):
#
    # Set up arrays/variables to feed in data from RADMC-3D outputs
#
    ngrid = int( float(nbins)**(3.) )
    r = [0.] * ngrid ; temp = []
#
    f_temp = open(arch_dir+'dust_temperature.dat')
#
    # Now read in data
#
    for i in range(0,3):
        trash = f_temp.readline()
    for lines in f_temp:
        lines = lines.strip() ; columns = lines.split()
        temp.append(float(columns[0]))
    f_temp.close()
#
    # Now fill arrays, comuputing locations
#
    count = 0
    for zz in range(0, nbins):
        for yy in range(0, nbins):
            for xx in range(0, nbins):
                x = binlims[0] + ( xx * bin_int ) + ( 0.5 * bin_int )
                y = binlims[0] + ( yy * bin_int ) + ( 0.5 * bin_int )
                z = binlims[0] + ( zz * bin_int ) + ( 0.5 * bin_int )
                r[count] = np.sqrt(x**2. + y**2. + z**2.) / r_dust
                count = count + 1
    r = np.array(r) ; temp = np.array(temp)
#
    # Short procedure to find number of unique radii from r array, using indices
    # of these to produce azimuthally averaged density and temperature profiles
#
    r_plt, indices = np.unique(r, return_inverse=True)
    temp_plt = [0.]*len(r_plt) ; count = [0.]*len(r_plt)

    for i in range(0, len(indices)):
        count[indices[i]] = count[indices[i]] + 1.
        temp_plt[indices[i]] = temp_plt[indices[i]] + temp[i]
    for i in range(0, len(count)):
        if (count[i] != 0):
            temp_plt[i] = temp_plt[i] / count[i]
        else:
            continue
#
    # If called for, read in DS output temperatures to plot in comparison
#
    if profile_comp is True:
        f_ds = open(ds_dir+'outcells.dat')
        trash = f_ds.readline()
        r_ds = [] ; temp_ds = []
        for lines in f_ds:
            lines = lines.strip() ; columns = lines.split()
            r_ds.append(float(columns[0])/r_dust)
            temp_ds.append(float(columns[6]))
        r_ds = np.array(r_ds) ; temp_ds = np.array(temp_ds)
        f_ds.close()
#
    # Now plot azimuthally averaged profiles
#
    plt.scatter(r_plt, temp_plt, label = 'Grid Output')
    plt.xlim(8e-1, 4e3)
    plt.ylim(2e1,3e3)
    plt.yscale('log') ; plt.xscale('log')
    plt.xlabel('Radius (AU)') ; plt.ylabel('Temperature (K)')
#
    if profile_comp is True:
        plt.plot(r_ds, temp_ds, color = 'r', label = 'DS test')
        plt.legend(loc = 'lower left', fontsize = 16, scatterpoints = 1)
#
    plt.savefig(plt_dir+'temp_r.png') ; plt.clf()
#
    return
