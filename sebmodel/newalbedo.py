import numpy as np
from numba import njit

from globals import *
from info import *


'''
===============================================================================

This file contains an albedo parameterisation along with its relevant subroutines

===============================================================================
'''

@njit
def initgrains(lid,nl):
    ''' 
    Subroutine that initializes the grainsizes 
    For now, only for Neumayer a spin has been performed to start with a realistic profile
    For all other stations, the initial grainsize is equal to the fresh snow grain size
    
    Parameters
    ----------
    lid : ndarray (1D float)
        layer ID (0 = ice, 1 = snow, 2 and larger = firn)
    nl : int
        Number of subsurface layers
        
    Returns
    -------
    grainsize: ndarray
        subsurface grainsize profile
    radfresh: float
        fresh snow grain radius
    '''

    grainsize = np.zeros(nlmax)
    
    if((chstation == "ant_neuma") & (SSAfresh) == 60):
        radfresh = 0.00025
    else: 
        radfresh = 3.0/(SSAfresh*densice)
    for il in range(0,nl):
        if (lid[il] == .0):
            grainsize[il] = radrefr
        else:
            grainsize[il] = radfresh

    return grainsize, radfresh


@njit
def newalbedo(sbuf,szenith,cloud,z,dz,lid,nlsnow,dtdz,dens,temp,water, mass, refrfrac, freshfrac,
              grainsize,radfresh,DENSVals,TVals,DTDZVals,TAUMAT,KAPMAT,DR0MAT):
    '''
    Routine that contains the albedo parameterisation as described by Kuipers Munneke et al., 2011 (PKM)
    which is in turn based on Gardner & Sharp, 2010 (G&S)

    Parameters
    ----------
    sbuf : ndarray (1D float)
        Array of interpolated forcing values at time t
    szenith : float
        solar zenith angle in radians
    cloud : float
        cloud cover
    z : ndarray (1D float)
        depth of layer centers, positive downwards
    dz : ndarray (1D float)
        layer thickness
    lid : ndarray (1D float)
        layer ID (0 = ice, 1 = snow, 2 and larger = firn)
    nlsnow : int
        Number of snow layers
    dtdz : ndarray (1D float)
        Vertical temperature gradient
    dens : ndarray (1D float)
        subsurface density profile
    temp : ndarray (1D float)
        subsurface temperature profile
    water : ndarray (1D float)
        subsurface liquid water content profile
    mass : ndarray (1D float)
        subsurface mass of layer profile
    refrfrac : ndarray (1D float)
        subsurface profile of refrozen snow
    freshfrac : ndarray (1D float)
        subsurface profile of fresh snow fraction 
    grainsize : ndarray (1D float)
        subsurface grainsize profile
    radfresh : float
        fresh snow grain radius
    DENSVals : ndarray (1D float)
        Values of density in table    
    TVals : ndarray (1D float)
        Values of temperature in table    
    DTDZVals : ndarray (1D float)
        Values of temperature gradient in table    
    TAUMAT : ndarray (3D float)
        Tabulated values of TAU
    KAPMAT : ndarray (3D float)
        Tabulated values of KAPPA
    DR0MAT : ndarray (3D float)
        Tabulated values of DR0
        
    Returns
    -------
    albedo: float
        shortwave surface albedo
    '''

    layeralbedo = np.zeros(100)

    grainsize = metamorphism(nlsnow,dtdz,dens,temp,water, mass, refrfrac, freshfrac, 
                             grainsize,radfresh,DENSVals,TVals,DTDZVals,TAUMAT,KAPMAT,DR0MAT)

    factor = (30./densice)**(-0.07)

    coszenith = np.cos(szenith)
    # !	tau = 1.14*(EXP(3.2*cloud)-1.0) !Kuipers Munneke et al., 2011, International Journal of Climatology, Eq .3
    # !	tau = 8.36*(EXP(1.548*cloud)-1.0)

    if(chstation == "alp_aws02"):
        tau = 0.0535*(np.exp(6.812*cloud)-1.0)
        
    elif(chstation == "ant_aws04"):
        tau = 6.541*(np.exp(2.065*cloud)-1.0)
        
    else:
        # These are actually Neumayer values
        # !	 tau = 4.58*(EXP(2.054*cloud)-1.0) !Lower bound
        # !	 tau = 6.228*(EXP(2.359*cloud)-1.0) !Upper bound
        tau = 5.404*(np.exp(2.207*cloud)-1.0)

    # Correction for zenith angle (G&S)
    if (coszenith < 1E-9):
        x = 0.
    else:
        x = min(np.sqrt(tau/(3.0*coszenith)),1.0) # Definition of x in G&S Eq. 10 / PKM Eq. 6
    tempzenith = 0.53*((1-0.64*x-(1-x)*coszenith)**1.2) # Part of PKM Eq. 6

    # Correction for clear sky optical thickness (PKM Eq. 11, converted to Pa)
    dalbclear = min(0.0,0.03247*np.log(sbuf[4]/153880.0)) #sbuf[4] = pressure

    il=0
    while ((z[il]+0.5*dz[il]) <= 0.5): #We go only 10 cm deep
        if(lid[il] == 0):
            layeralbedo[il] = albice
        else:
            basealb = 1.48 - factor*(grainsize[il]**0.07) #PKM Eq. 5
        
            #Correction due to loading by light-absorbing carbon (G&S Eq. 8, modified for grain size in metres)
            dalbcarbon = max(0.04-basealb,(-(soot**0.55))/(0.16+np.sqrt((0.012*densice)/grainsize[il])+\
                1.8*(soot**0.6)*((densice/30.)**0.25)*(grainsize[il]**0.25)))
            
            tempsum = basealb + dalbcarbon
            
            dalbzenith = basealb*(1.-tempsum)*tempzenith #Rest of PKM Eq. 6
            
            dalbtau = (0.1*tau*(tempsum**1.3))/((1.0+1.5*tau)**basealb) #PKM Eq. 8
            
            layeralbedo[il] = basealb + dalbzenith + dalbtau + dalbclear #PKM Eq. 4 + Eq. 11

        if (il==0):
            albedo = layeralbedo[il]
        else:
            albedo = albedo + ((layeralbedo[il]-layeralbedo[il-1])*np.exp(-(z[il]+0.5*dz[il])/0.01))
        il = il + 1
    if (albedo > albmax):
        albedo = albmax
    elif (albedo < albmin):
        albedo = albmin

    return albedo

@njit
def metamorphism(nlsnow,dtdz,dens,temp,water,mass,refrfrac,freshfrac,
                 grainsize,radfresh,DENSVals,TVals,DTDZVals,TAUMAT,KAPMAT,DR0MAT):
    '''
    Routine that calculates the metamorphism of grainsizes due to old snow, refrozen snow and new snow

    Parameters
    ----------
    nlsnow : int
        Number of snow layers
    dtdz : ndarray (1D float)
        Vertical temperature gradient
    dens : ndarray (1D float)
        subsurface density profile
    temp : ndarray (1D float)
        subsurface temperature profile
    water : ndarray (1D float)
        subsurface liquid water content profile
    mass : ndarray (1D float)
        subsurface mass of layer profile
    refrfrac : ndarray (1D float)
        subsurface profile of refrozen snow
    freshfrac : ndarray (1D float)
        subsurface profile of fresh snow fraction 
    grainsize : ndarray (1D float)
        subsurface grainsize profile
    radfresh : float
        fresh snow grain radius
    DENSVals : ndarray (1D float)
        Values of density in table    
    TVals : ndarray (1D float)
        Values of temperature in table    
    DTDZVals : ndarray (1D float)
        Values of temperature gradient in table    
    TAUMAT : ndarray (3D float)
        Tabulated values of TAU
    KAPMAT : ndarray (3D float)
        Tabulated values of KAPPA
    DR0MAT : ndarray (3D float)
        Tabulated values of DR0
        
    Returns
    -------
    grainsize: ndarray (1D float)
        vertical profile of grainsize
    '''

    drdry = drysnow(nlsnow,dtdz,dens,temp,grainsize,radfresh,DENSVals,TVals,DTDZVals,TAUMAT,KAPMAT,DR0MAT)

    for il in range (0,nlsnow):
        watcont = water[il]/(mass[il]+water[il])
        if (watcont > 1):
            raise ValueError("Water content greater than one...")

        drwet = (4.22E-13*(watcont**3.))/(4.0*np.pi*grainsize[il]**2)*tstep # PKM Eq. 2

        oldfrac = 1 - refrfrac[il] - freshfrac[il]
        if ((oldfrac > 1.0) | (refrfrac[il] > 1.0) | (freshfrac[il]> 1.0)):
            print("Some fraction is greater than 1", oldfrac, refrfrac[il], freshfrac[il])
            # raise ValueError("Some fraction is greater than 1", oldfrac, refrfrac[il], freshfrac[il])
        elif ((oldfrac < 0.0) | (refrfrac[il] < 0.0) | (freshfrac[il] < 0.0)):
            print("Some fraction is smaller than 0", oldfrac, refrfrac[il], freshfrac[il])
            # raise ValueError("Some fraction is smaller than 0", oldfrac, refrfrac[il], freshfrac[il])

        grainsize[il] = (grainsize[il]+drdry[il]+drwet)*oldfrac+refrfrac[il]*radrefr+freshfrac[il]*radfresh

        # Reset fractions of refrozen and fresh snow
        refrfrac[il] = 0.0
        freshfrac[il] = 0.0

        if (grainsize[il] < radfresh):
            if ((grainsize[il]-radfresh) <= (-1.0E-5)):
                # raise ValueError("Grainsize too small!")
                print("Grainsize too small!")
            else:
                if (lcomment == 1):
                    print("Watch out, grainsize slightly smaller than radfresh:")
                    print("Difference ", (radfresh-grainsize[il]))
                    
    return grainsize

@njit
def drysnow(nlsnow,dtdz,dens,temp,grainsize,radfresh,DENSVals,TVals,DTDZVals,TAUMAT,KAPMAT,DR0MAT):
    '''
    Routine that determines the dry snow metamorphism
    It uses the look-up tables as initialized in inittables.py

    Parameters
    ----------
    nlsnow : int
        Number of snow layers
    dtdz : ndarray (1D float)
        Vertical temperature gradient
    dens : ndarray (1D float)
        subsurface density profile
    temp : ndarray (1D float)
        subsurface temperature profile
    grainsize : ndarray (1D float)
        subsurface grainsize profile
    radfresh : float
        fresh snow grain radius
    DENSVals : ndarray (1D float)
        Values of density in table    
    TVals : ndarray (1D float)
        Values of temperature in table    
    DTDZVals : ndarray (1D float)
        Values of temperature gradient in table    
    TAUMAT : ndarray (3D float)
        Tabulated values of TAU
    KAPMAT : ndarray (3D float)
        Tabulated values of KAPPA
    DR0MAT : ndarray (3D float)
        Tabulated values of DR0
        
    Returns
    -------
    drdry: ndarray (1D float)
        vertical profile of change in grainsize of dry snow
    '''

    absdtdz = np.abs(dtdz)
    
    drdry = np.zeros(nlsnow)

    for il in range(0,nlsnow):
        if(dens[il] < DENSVals[1]):
            idens = 1
            tempdens = DENSVals[1]
        elif (dens[il] >= DENSVals[8]):
            idens = 7
            tempdens = DENSVals[8]
        else:
            idens = 1
            while (DENSVals[idens] <= dens[il]):
                idens = idens + 1
            idens = idens - 1

        if(temp[il] < TVals[1]):
            itemp = 1
            temptemp = TVals[1]
        elif(temp[il] >= TVals[11]):
            itemp = 10
            temptemp = TVals[11]
        else:
            itemp = 1
            while(TVals[itemp] <= temp[il]):
                itemp = itemp + 1
            itemp = itemp - 1

        if(absdtdz[il] < DTDZVals[1]):
            idtdz = 1
            tempdtdz = DTDZVals[1]
        elif(absdtdz[il] >= DTDZVals[31]):
            idtdz = 30
            tempdtdz = DTDZVals[31]
        else:
            idtdz = 1
            while (DTDZVals[idtdz] <= absdtdz[il]):
                idtdz = idtdz + 1
            idtdz = idtdz - 1

        fracdens = (dens[il]-DENSVals[idens])/(DENSVals[idens+1]-DENSVals[idens])
        fractemp = (temp[il]-TVals[itemp])/(TVals[itemp+1]-TVals[itemp])
        fracdtdz = (absdtdz[il]-DTDZVals[idtdz])/(DTDZVals[idtdz+1]-DTDZVals[idtdz])

        # Now retrieve the values for tau, kap and dr0 from the look-up tables
        TAU = 0.0
        KAP = 0.0
        DR0 = 0.0

        TAU = TAU + (1.0 - fracdens) * (1.0 - fractemp) * \
        (1.0 - fracdtdz) * TAUMAT[itemp  ,idtdz  ,idens  ]
        TAU = TAU + (1.0 - fracdens) * (1.0 - fractemp) * \
        fracdtdz           * TAUMAT[itemp  ,idtdz+1,idens]
        TAU = TAU + (1.0 - fracdens) * fractemp           * \
        (1.0 - fracdtdz) * TAUMAT[itemp+1,idtdz  ,idens]
        TAU = TAU + (1.0 - fracdens) * fractemp           * \
        fracdtdz           * TAUMAT[itemp+1,idtdz+1,idens]
        TAU = TAU + fracdens           * (1.0 - fractemp) * \
        (1.0 - fracdtdz) * TAUMAT[itemp  ,idtdz  ,idens+1]
        TAU = TAU + fracdens           * (1.0 - fractemp) * \
        fracdtdz           * TAUMAT[itemp  ,idtdz+1,idens+1]
        TAU = TAU + fracdens           * fractemp           * \
        (1.0 - fracdtdz) * TAUMAT[itemp+1,idtdz  ,idens+1]
        TAU = TAU + fracdens           * fractemp           * \
        fracdtdz           * TAUMAT[itemp+1,idtdz+1,idens+1]

        KAP = KAP + (1.0 - fracdens) * (1.0 - fractemp) * \
        (1.0 - fracdtdz) * KAPMAT[itemp  ,idtdz  ,idens  ]
        KAP = KAP + (1.0 - fracdens) * (1.0 - fractemp) * \
        fracdtdz           * KAPMAT[itemp  ,idtdz+1,idens]
        KAP = KAP + (1.0 - fracdens) * fractemp           * \
        (1.0 - fracdtdz) * KAPMAT[itemp+1,idtdz  ,idens]
        KAP = KAP + (1.0 - fracdens) * fractemp           * \
        fracdtdz           * KAPMAT[itemp+1,idtdz+1,idens]
        KAP = KAP + fracdens           * (1.0 - fractemp) * \
        (1.0 - fracdtdz) * KAPMAT[itemp  ,idtdz  ,idens+1]
        KAP = KAP + fracdens           * (1.0 - fractemp) * \
        fracdtdz           * KAPMAT[itemp  ,idtdz+1,idens+1]
        KAP = KAP + fracdens           * fractemp           * \
        (1.0 - fracdtdz) * KAPMAT[itemp+1,idtdz  ,idens+1]
        KAP = KAP + fracdens           * fractemp           * \
        fracdtdz           * KAPMAT[itemp+1,idtdz+1,idens+1]

        DR0 = DR0 + (1.0 - fracdens) * (1.0 - fractemp) * \
        (1.0 - fracdtdz) * DR0MAT[itemp  ,idtdz  ,idens  ]
        DR0 = DR0 + (1.0 - fracdens) * (1.0 - fractemp) * \
        fracdtdz           * DR0MAT[itemp  ,idtdz+1,idens]
        DR0 = DR0 + (1.0 - fracdens) * fractemp           * \
        (1.0 - fracdtdz) * DR0MAT[itemp+1,idtdz  ,idens]
        DR0 = DR0 + (1.0 - fracdens) * fractemp           * \
        fracdtdz           * DR0MAT[itemp+1,idtdz+1,idens]
        DR0 = DR0 + fracdens           * (1.0 - fractemp) * \
        (1.0 - fracdtdz) * DR0MAT[itemp  ,idtdz  ,idens+1]
        DR0 = DR0 + fracdens           * (1.0 - fractemp) * \
        fracdtdz           * DR0MAT[itemp  ,idtdz+1,idens+1]
        DR0 = DR0 + fracdens           * fractemp           * \
        (1.0 - fracdtdz) * DR0MAT[itemp+1,idtdz  ,idens+1]
        DR0 = DR0 + fracdens           * fractemp           * \
        fracdtdz           * DR0MAT[itemp+1,idtdz+1,idens+1]

        if((KAP <= 0.0) | (TAU <= 0.0)):
            drdry[il] = 0.0
        else:
            if (np.abs(grainsize[il] - radfresh) < 1.0E-5):
                drdry[il] = (DR0*tstep*(1E-6)*((TAU/(TAU+1.0))**(1./KAP)))/3600.
            else:
                drdry[il] = (DR0*tstep*(1E-6)*((TAU/(TAU+1E6*(grainsize[il]-radfresh)))**(1./KAP)))/3600.
        
    return drdry

