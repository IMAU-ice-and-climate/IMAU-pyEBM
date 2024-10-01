import numpy as np
from numba import njit

from globals import *
from info import *
from routines import irreducible

# !===============================================================================
# !
# ! All routines related to the calculation of the (initial) density and temperature profiles
# !
# !===============================================================================

@njit
def snowtemp(z):
# !===============================================================================
    if chstation == 'ant_domeA':
        if (z < 0.82):
            temp = -48.56 + 273.15
        elif (z <= 10):
            temp = -1.44e-05*z**6. + 0.0002271*z**5 + 0.008398*z**4 - 0.2618*z**3. + 2.703*z**2. -11.85*z - 40.52 + 273.15
        else:
            temp =   T10m  + 273.15
    else:
        temp = T10m + 273.15
    if(temp >= 273.0): temp = 273.0

    return temp


@njit
def snowdens(z,id):
# !===============================================================================

    if chstation == 'ant_domeA':
        if (z <= 10):
            dens =  -0.1292 * z**4 + 3.2804* z**3 -28.5735* z**2 + 102.9793*z + 287.4564 
        else:
            dens = 450.0
    else:
        dens = densfirn +(rhosninit-densfirn)*np.exp(-z/10.)

    if ((dens > densice) | (id == 0)):
        dens = densice
    elif(dens < densnow):
        dens = densnow

    return dens

@njit
def resetdens(z, lid, dens, mass, dz, nl): 
    # !===============================================================================
    # ! Resets the sub-surface profiles of density and mass profiles
    # !===============================================================================
    # USE SNOW_EBM , ONLY : nl, z, dz, dens, mass , lid
    # USE INPUT_EBM , ONLY : chstation

    # IMPLICIT NONE
    # !Local
    # INTEGER :: il
    # REAL :: snowdens, snowdens_STGL

    for il in range(0,nl):
        if (chstation != "STGL_1999"): 
            dens[il] = snowdens(z[il],lid[il])

        elif (chstation == "STGL_1999"):
            print('snowdens_STGL not implementend')
            # if (il > 0) :
            #     dens[il] = snowdens_STGL(z[il],dz[il],z(il-1),dz(il-1),il,lid[il])
            # else:
            #     dens[il] = snowdens_STGL(z[il],dz[il],0.,0.,il,lid[il])

        mass[il] = dens[il]*dz[il]
    
    irrwater = irreducible(dens,nl)

    return dens, mass, dz, irrwater

