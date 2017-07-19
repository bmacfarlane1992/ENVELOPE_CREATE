'''

 main.py

 Procedure written to construct a distribution of SPH particles
 into an envelope distribution, as per benchmarking tests of Ivezic et al. (1997)
 and Bjorkmann & Wood (2001)

 Program outputs SPH properties as regular grid, as per requirements of RADMC-3d,
 using kernel information to smooth SPH properties consistently over the space.

 Author: Benjamin MacFarlane
 Date: 19/07/2017
 Contact: bmacfarlane@uclan.ac.uk

'''
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
    # System Variables
#
arch_dir = "/Users/bmacfarlane1992/PHD/PROJECTS/ENVELOPE_CREATE/"
#
    # Variables for envelope creation
#
envelope_create = True      # Choose whether or not to create new SPH -> grid translated model
#
tau = [1.,10.,100.]                  # Optical depth at 1 micron
p=[0,2]                       # Exponent on density profile
#
    # Variables for gridding
#
grid_method = 'reg'         # ['reg', 'oct']
dens_method = 'sphinf'     # For 'reg' grid_method only: ['part','sphinf','brute','byhand']
#
    # Variables for RADMC-3D outputs
#
therm = True        # Decide whether run RADMC-3D MC therm program
sed = True          # Decide whether SED is generated and plotted
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - CONSTANTS AND CONVERSIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
Rsol_to_cm = 6.957e10
r_star = 1
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import os
import sys
#
import envelope
import gridmake
import densgrid
import inp_gen
import radmc3dPy
import profiles
import sedplt
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
    # Loop over both tau and p values, to batch produce data and outputs
#
for tau_ind in range(0, len(tau)):
    for p_ind in range(0, len(p)):
#
    # First, create p, tau and dens_method specific plot directory for ouptuts,
    # and isolate where DS data is for profile/SED comparisons
#
        plt_dir = arch_dir+'PLOTS/'+str(grid_method)+'_p'+\
           str(int(p[p_ind]))+'_tau'+str(int(tau[tau_ind]))+'/'
        os.system('mkdir '+plt_dir)
#
        ds_dir = arch_dir+'/DS_outputs/p'+str(int(p[p_ind]))+'/tau' \
           +str(int(tau[tau_ind]))+'/'

### ------------------------------------------------------------------------ ###

        code_dir = arch_dir+"/Code/"
        os.chdir(code_dir)

### ------------------------------------------------------------------------ ###

#
    # Create envelope SPH setup
#
        if envelope_create is True:
#
            pos, rho, rho_0, h, m_part, r, r_dust \
               = envelope.create(arch_dir, r_star, tau[tau_ind], p[p_ind],
               plt_dir, Rsol_to_cm)

### ------------------------------------------------------------------------ ###

#
    # Create grid to which density is translated from SPH distribution.
    # Method of gridding either:
    #           'reg' = Regular gridding with n_bins in each dimension defined in gridmake
    #           'oct' = Build octree through hyperion modules
#
            if grid_method == 'reg':
                nbins, bin_it, binlims, grid = gridmake.reg(arch_dir, pos)
#
            elif grid_method == 'oct':
                o, n, ncells = gridmake.oct(arch_dir, pos, rho, h, m_part, r_dust)
#
            else:
                print "Incorrect selection of grid_make variable in main.py"
                exit()

### ------------------------------------------------------------------------ ###

#
    # For populate bin density contributions:
    #
    #   For 'reg' gridding, work must be done, with method based on dens_method
    #
    #       'part' = When particle found in bin, add central density to bin
    #       'sphinf' = Add kernel contribution of particles in relevant sphere of influence
    #       'brute' = Add kernel contributions for all particles that intersect bin centre
    #       'byhand' = By-hand set the density for every bin in p=0 case, where density must be rho_0
    #
    #   For 'oct' gridding, work is done in gridmake.py through hyperion modules
    #   so work is only in identifying leaves, and writing to file
#
            if grid_method == 'reg':
                if dens_method == 'part':
                    rhogrid = densgrid.particle(arch_dir, pos, rho, m_part, h, \
                       grid, bin_it, binlims, nbins)
                elif dens_method == 'sphinf':
                      rhogrid = densgrid.kernel(arch_dir, pos, rho, m_part, h, \
                      grid, bin_it, binlims, nbins)
                elif dens_method == 'brute':
                    rhogrid = densgrid.brute(arch_dir, pos, rho, m_part, h, \
                       grid, bin_it, binlims, nbins)
                elif (dens_method == 'byhand'):
                    rhogrid = densgrid.byhand(arch_dir, grid, binlims, nbins, \
                       rho_0, r_dust, p[p_ind])
                else:
                    print "\nIncorrect selection of density gridding (dens_method)\n"
                    exit()
#
            elif grid_method == 'oct':
                densgrid.oct(arch_dir, o, n, ncells)
#
            else:
                print "\nIncorrect selection of gridding geometry (grid_method)\n"
#
#
### ------------------------------------------------------------------------ ###
#
#
    # Check the radial density profile output from SPH -> grid method
#
            if grid_method == 'reg':
                profiles.plotrho(arch_dir, bin_it, binlims, nbins, r, r_dust, rho, \
                   rho_0, ds_dir, p[p_ind], tau[tau_ind], plt_dir)
#
#
### ------------------------------------------------------------------------ ###
#
#
        if therm is True:
#
    # Generate the dustkappa_silicate.inp and wavelength_micron.inp files as
    # per Ivezic et al. (1997) benchmark test requirements
    # Also generate grid-independent input files for RADMC-3D, before mctherm run:
    #           dustopac.inp
    #           wavelength_micron.inp
    #           stars.inp
    #           lines.inp
    #           radmc3d.inp

#
            t_star = inp_gen.all(arch_dir, code_dir, r_star, Rsol_to_cm)
#
#
### ------------------------------------------------------------------------ ###
#
#
    # Generate the temperature structure
#
            os.chdir(arch_dir)
            os.system('radmc3d mctherm')
#
#
### ------------------------------------------------------------------------ ###
#
#
    # Check the radial profiles of the temperature that is outputted
    # from SPH -> grid -> MC therm test
#
            if grid_method == 'reg':
                profiles.plottemp(arch_dir, bin_it, binlims, nbins, r_dust, \
                   ds_dir, p[p_ind], tau[tau_ind], plt_dir)
#
#
### ------------------------------------------------------------------------ ###
#
#
    # Generate SED, then read in data from spectrum.out file and plot
#
        if sed is True:
            os.system('rm -r '+arch_dir+'lines.inp')
            os.system('radmc3d sed incl 0')
#
            sedplt.plot(arch_dir, r_star, Rsol_to_cm, t_star, \
               p[p_ind], tau[tau_ind], plt_dir, ds_dir)
#
#
### ------------------------------------------------------------------------ ###

os.chdir(arch_dir+"Code/")
os.system("rm -r *.pyc") ; os.system("rm -r *~")

### ------------------------------------------------------------------------ ###

exit()
