'''

opac_wav.py

Routine to generate both the dustkappa_silicate.dat opacity table and
wavelength_micron.inp files, using the definitions for scattering/absorption
as per Bjorkmann & Wood (2001) with user-defined wavelength range/binning.

Author: Benjamin MacFarlane
Date: 13/03/2017
Contact: bmacfarlane@uclan.ac.uk

'''
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import os
import numpy as np
import math
import radmc3dPy
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
    # Stellar variables
#
m_star = 1
t_star = 2500.
#
    # Opacity file specific variables
#
select = 1          # Defined routine to select wavelength: 1 - User defined, logarithmic spaced
#                                                           2 - As per Bjorkmann & Wood (2001)
wavlims = np.array( [0.15, 2500] )   # Wavelength (units of microns) limits
wavbins = 256.          # No. of bins in wavelenght_micron.inp file
#
    # Other variables
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - CONSTANTS AND CONVERSIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
c = 3e8
h = 6.626e-34
k_sb = 1.38e-23
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def all(arch_dir, code_dir, r_star, Rsol_to_cm):

### ------------------------------------------------------------------------ ###

#
    # wavelength_micorn.inp and dustkappa_silicate.inp
#
    f1 = open(arch_dir+'wavelength_micron.inp', 'w')
    f1.write(str(int(wavbins))+'\n')

    f2 = open(arch_dir+'dustkappa_silicate.inp','w')
    f2.write('2\n')
    f2.write(str(int(wavbins))+'\n')

    if (select == 1):
        d0 = math.log(wavlims[0],10) ; d1 = math.log(wavlims[1],10)
        dinc = (d1 - d0) / wavbins
    elif (select == 2):
        nu_max = (16. * k_sb * t_star) / h
        del_nu = nu_max / wavbins
    llambda = []
    for i in range(0, int(wavbins)):

        if (select == 1):
            wav = 10**( d0 + dinc*float(i) )
        elif (select == 2):
            wav = c / (nu_max - (float(i) * del_nu))
            wav = wav / 1e-6
        llambda.append(wav)

        f1.write(str(wav)+'\n')

        if (wav <= 1.):
            absorp = 0.5
            scatt = 0.5
        elif (wav > 1.):
            absorp = 0.5 * (1. / wav)
            scatt = 0.5 * (1. / (wav)**4. )

        f2.write(str(wav)+' '+str(absorp)+' '+str(scatt)+'\n')

    f1.close()
    f2.close()

### ------------------------------------------------------------------------ ###

#
    # dustopac.inp
#
    radmc3dPy.analyze.writeDefaultParfile('ppdisk')
    radmc3dPy.setup.problemSetupDust('ppdisk')
    os.system('rm -r amr_grid.inp dust_density.binp \
       problem_params.inp radmc3d.inp stars.inp wavelength_micron.inp')
    os.system('mv '+code_dir+'dustopac.inp '+arch_dir)

### ------------------------------------------------------------------------ ###

#
    # stars.inp
#
    f = open(arch_dir+'stars.inp','w')
    f.write('2 \n')
    f.write('1 '+str(int(wavbins))+'\n')
    f.write(str(r_star*Rsol_to_cm)+' '+str(m_star)+' 0. 0. 0. \n')
    for i in range(0, int(wavbins)):
        f.write(str(llambda[i])+'\n')
    f.write('-'+str(t_star)+'\n')
    f.close()

### ------------------------------------------------------------------------ ###

#
    # lines.inp
#
    f = open(arch_dir+'lines.inp','w')
    f.write('2 \n')
    f.write('1 \n')
    f.write('co leiden 0 0 0 \n')
    f.close()

### ------------------------------------------------------------------------ ###

#
    # radmc3d.inp
#
    f = open(arch_dir+'radmc3d.inp','w')
    f.write('nphot = 1000000 \n')
    f.write('scattering_mode_max = 1 \n')
    f.write('scattering_mode= 0\n')
    f.write('lines_mode = 1\n')
    f.write('modified_random_walk = 0 \n')
    f.write('tgas_eq_tdust = 1 \n')
    f.close()

### ------------------------------------------------------------------------ ###

    return t_star
