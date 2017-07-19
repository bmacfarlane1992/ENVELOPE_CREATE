'''

envelope.py

Routine to create envelope is inputs to RADMC-3D benchmarking test.
Envelope distribution set by Ivezic et al. (1997) parameters, noted in
Bjorkmann & Wood (2001).

Algorithm to generate envelope closely follows method of Stamatellos (2003) disc
construction, using random numbers to distribute particles and associate SPH
properties.

Author: Benjamin MacFarlane
Date: 04/05/2017
Contact: bmacfarlane@uclan.ac.uk

'''
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
nsph = int(1e6)         # Number of SPH particles in model
#
N_neigh = 50.     # Number of neighbours that SPH distributes particle properties over
#
rdust_rstar = [[8.44, 8.46, 8.60],[9.11, 11.37, 17.67]] # rdust/rstar as a function of p and tau
#                                                       # [ [ (p_1,tau_i) , ... ], [ (p_2,tau_i) , ... ]]
seren_check = False         # Output DRAGON formatted file to check rho(r) profile in SPLASH
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - CONSTANTS AND CONVERSIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
cgsdens_to_sidens = 1e3
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
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def create(arch_dir, r_star, tau, p, plt_dir, Rsol_to_cm):
#
    print "\nCreating envelope for p="+str(int(p))+", tau="+str(int(tau))+" model\n"
#
    # Initialise arrays both for particle data and random seeds
#
    pos = [ [[0.] for i in range(nsph)] for i in range(3)]
    r = [0.] * nsph ; theta = [0.] * nsph ; phi = [0.] * nsph
    rho = [0.] * nsph ; h = [0.] * nsph
    Rs = [[[0.] for i in range(nsph)] for i in range(3)]
#
    # Generate random seeds for spatial distributions
#
    for i in range(0, nsph):
        Rs[0][i] = random.random()
        Rs[1][i] = random.random()
        Rs[2][i] = random.random()
#
        # Using definitions of Ivezic et al. (1997) define spatial limits
        # of envelope
#
    if ( p == 0. ):
        if (tau == 1.):
            r_dust = rdust_rstar[0][0] * r_star * Rsol_to_cm
        elif (tau == 10.):
            r_dust = rdust_rstar[0][1] * r_star * Rsol_to_cm
        elif (tau == 100.):
            r_dust = rdust_rstar[0][2] * r_star * Rsol_to_cm
        else:
            print "incorrect tau/p combination entered"
            exit()
    if ( p == 2. ):
        if (tau == 1.):
            r_dust = rdust_rstar[1][0] * r_star * Rsol_to_cm
        elif (tau == 10.):
            r_dust = rdust_rstar[1][1] * r_star * Rsol_to_cm
        elif (tau == 100.):
            r_dust = rdust_rstar[1][2] * r_star * Rsol_to_cm
        else:
            print "incorrect tau/p combination entered"
            exit()
#
    r_env = 1.e3 * r_dust
#
    rho_0 = tau / (1. - (r_dust / r_env))
    if (p == 0):
        rho_0 = rho_0 / r_env
    elif (p == 2):
        rho_0 = rho_0 / r_dust
#
    # Define omega value for envelope from which randomised spatial distribution
    # can be generated with selected radial density profile. Also compute
    # starting density and envelope mass
#
    exp = 1. - (p / 3.)
    w_env = r_env**(3.) / r_dust**(3.)
#
    for i in range(0, nsph):
#
    # Use omega value computed above, and generate radial specfic omega value.
    # With this value, use random seed to find particle location
#
        w_r = ( ( ( w_env**exp - 1.) * Rs[0][i] ) + 1. )**(1./exp)
        r[i]= r_dust * ( (w_r)**(1./3.) )
#
    # Use Rs_2 to distribute the particles through theta and phi, taking
    # relevent angular ranges into account in random seeds
#
        theta[i] = math.acos(1. - (2. * Rs[1][i]))
        phi[i] = 2. * math.pi * Rs[2][i]
#
    # Convert from spherical to cartesian co-ordinates
#
        pos[0][i] = ( r[i] * (math.sin(theta[i]) * math.cos( phi[i] ) ) )
        pos[1][i] = ( r[i] * ( math.sin(theta[i]) * math.sin( phi[i] ) ) )
        pos[2][i] = ( r[i] * math.cos(theta[i]))
#
    # Using neighbours, calculate smoothing length of each particle, using
    # particle density and mass information
#
        m_env = ( (4. * math.pi * rho_0 * r_dust**(3.)) / (3. * exp) ) \
           * (w_env**(exp) - 1.)
        m_part = m_env / nsph
#
        rho[i] = rho_0 * (r_dust / r[i])**(p)
        h[i] = ( (3. * N_neigh * m_part) / (32. * math.pi * rho[i]))**(0.333)
#
    # If called for, output DRAGON file to evaluate SPH distribution with SEREN
#
    if seren_check is True:
        f = open(arch_dir+'/seren/envelope/envelope.dat','w')
        f.write(str(nsph)+'\n')
        for i in range(0, 69):
            f.write('0\n')
        for i in range(0, nsph):
            f.write(str(pos[0][i])+' '+str(pos[1][i])+' '+str(pos[2][i])+'\n')
        for i in range(0, nsph):
            f.write('1 1 1 \n')         # Mask for v arrays
        for i in range(0, nsph):
            f.write('10 \n')             # Mask for T array
        for i in range(0, nsph):
            f.write(str(h[i])+'\n')
        for i in range(0, nsph):
            f.write(str(rho[i])+'\n')
        for i in range(0, nsph):
            f.write(str(m_part)+'\n')
        for i in range(0, nsph):
            f.write('1 \n')             # Mask for itype array
        f.close()
#
    # Set up arrays to, and then plot smoothing length radial profile
#
    r_plt = [] ; h_plt = []
    for i in range(0, len(r)):
        r_plt.append(r[i]/r_dust)
        h_plt.append(2.*h[i]/r_dust)
    r_plt = np.array(r_plt) ; h_plt = np.array(h_plt)
#
    fig = plt.figure(1)
    ax1 = plt.subplot(111)
    plt.scatter(r_plt, h_plt, color = 'b')
#
    plt.xlabel('r/'+(r'R$_{dust}$'), fontsize = 18, labelpad=0.5)
    plt.xticks(fontsize = 15) ;   ax1.set_xscale('log') ; ax1.set_xlim(1,1200.)
    plt.ylabel('2h(r)/'+(r'R$_{dust}$'), fontsize = 18, labelpad=0.5)
    ax1.set_yscale('log')
    plt.savefig(plt_dir+'h_r.png') ; plt.clf()
#
    return pos, rho, rho_0, h, m_part, r, r_dust
