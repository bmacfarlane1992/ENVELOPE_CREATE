'''

sedplt.py

Routine to plot spectral energy distributions output as spectrum.out from
RADMC-3D. Can also compare outputs to DS tests, based on the Ivzic et al. (1997)
benchmark tests.

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
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - CONSTANTS AND CONVERSIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
sb_const_cgs = 5.6705e-5
pc_cgs = 3.086e18

#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def plot(arch_dir, r_star, Rsol_to_cm, t_star, p, tau, plt_dir, ds_dir):
#
    wav = [] ; flux = [] ; flam = []
    f = open(arch_dir+'spectrum.out')
#
    # Read in data, converting F_nu to lambda F_lambda
#
    for j in range(0, 3):
        header = f.readline()
    for lines in f:
        lines = lines.strip() ; columns = lines.split()
        wav.append(float(columns[0]))
        flam.append( float(columns[1]) * ( 3e10 / (float(columns[0])*1e-4) ) )
#
    # Compute total flux from central source, to normalise SED
#
    fnorm = ( (r_star*Rsol_to_cm)**2. * sb_const_cgs * (t_star)**4. ) / (pc_cgs)**2.
    for j in range(0, len(flam)):
        flam[j] = flam[j] / fnorm
#
    # Now plot, using normalised SED, comparing to DS outputs if called
#
    if profile_comp is True:
        f_ds = open(ds_dir+'sed0.0.0.dat')
        trash = f_ds.readline()
        wav_ds = [] ; flam_ds = []
        for lines in f_ds:
            lines = lines.strip() ; columns = lines.split()
            wav_ds.append(float(columns[0]))
            flam_ds.append(float(columns[1]))
        wav_ds = np.array(wav_ds) ; flam_ds = np.array(flam_ds)
        f_ds.close()
#
    fig = plt.figure(1)
    ax1 = plt.subplot(111)
    plt.scatter(wav, flam, label = 'BM test')
#
    if profile_comp is True:
        plt.plot(wav_ds, flam_ds, color = 'r', label = 'DS test')
        plt.legend(loc = 'lower center', fontsize = 16, scatterpoints = 1)
#
    plt.xlabel('Wavelength ('+(r'$\mu$m')+')', fontsize = 18, labelpad=0.5)
    plt.xticks(fontsize = 15) ;   ax1.set_xscale('log') ; ax1.set_xlim(0.2,2500.)
    plt.ylabel(r'$\lambda$ F$_{\lambda}$ / F', fontsize = 18, labelpad=0.5)
    ax1.set_yscale('log') ; ax1.set_ylim( 1e-5 , 1)
    plt.savefig(plt_dir+'SED.png') ; plt.clf()

    return
