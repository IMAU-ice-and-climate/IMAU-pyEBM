from numba import njit
import numpy as np

from globals import *
from info import *

# !===============================================================================
# !
# ! All routines related to the correction and calculation of radiation 
# ! components of the energy balance
# !
# !===============================================================================
@njit
def calcalbedo(szenith,lid,hmass,dsnowr,rhosn,water,precip,t0,snowdays,alb_old):
# !===============================================================================
# ! Routine that calculates the surface albedo based on:
# ! 1= Oerlemans and Knap, 1998, J. Glaciol. 
# ! 2= Douville et al 1995. 
# ! 3= combination of Oerlemans and Knap, 1998, J. Glaciol. and Zuo and Oerlemans, 1996, J. Glaciol. taken from Bougamont et al 2005
# ! including a dependence on snowdepth following Oerlemans and Knap, 1998 
# ! including solar zenith angle dependency: Segal et al., 1991 J. Atm. Sci., 48(8), 1024-1042. 
# !===============================================================================
# USE INPUT_EBM , ONLY : tstep , lalbedo , albsnow , albice , albfirn , snowstar , &
# &               tstardry0 , tstardry10  , tstarwet , solzenyes 
# USE GLOBALS_EBM , ONLY : snowdays , t0 , szenith , alb_old , albedo
# USE SNOW_EBM , ONLY : lid , dsnowh , water, nl , precip , dsnowr , hmass , rhosn
# USE CONSTANTS_EBM , ONLY : Tkel , secday , pi , errorval

# IMPLICIT NONE

# REAL , PARAMETER :: efoldcoef = 0.24			! Douville: exponential e-folding coefficient
# REAL , PARAMETER :: reducrate=0.008		! Douville: constant reduction ratio
# INTEGER :: il , firnrest
# REAL :: snowrest, snowwet
# REAL :: albmin, albmax
# REAL :: snowalb,tstar,constant
# REAL :: coszen, factor

    efoldcoef = 0.24 # Douville: exponential e-folding coefficient 
    reducrate=0.008 # Douville: constant reduction ratio
    # solzenyes  1 = yes, 0 = no correction albedo for low solar zenith angle
    if (solzenyes == 1):	
        coszen = np.cos(szenith)
        if (np.cos(szenith) < np.cos(80.*np.pi/180.)): coszen = np.cos(80.*np.pi/180.)
        factor = 0.32*0.5*((3./(1.+4.*coszen))-1.)
        if (factor < 0.): factor = 0.
    else:
        factor = 0.

    il = 0
    while (lid[il] == 1): 		#snow cover
        il = il + 1

    firnrest = 0
    if (lid[il]>= 1): firnrest = 1	#snow or firn

    albmax = albsnow
    snowrest = hmass		# snowrest in mm w.e.
    if ((snowrest <=0.) & (dsnowr > 0.)): snowrest = dsnowr*rhosn		# snowrest in mm w.e.
    if (snowrest < 0.): snowrest = 0.
    if (lid[0] >= 2): snowrest = 0.001		# firn / snow more than 1 year old
    snowwet = water[0]
    constant = (tstardry10-tstardry0)/10.

    if (precip <= 0.):		# in case of no precipitation snow albedo decreases in time due to aging of snow
        snowdays = snowdays + tstep/secday				# time since last snowfall in days.
        albmin = albfirn

        if (snowrest > 0.):		# there is still snow on the surface 

            if (lalbedo <= 1):					# Oerlemans and Knap
                snowalb = albmin + (albmax - albmin) * np.exp((-snowdays)/tstarwet)

            elif (lalbedo == 2):		# in Bougamont et al., 2005 using Oerlemans and Knap, 1998, Zuo and Oerlemans, 1996
                if (snowwet > 0.):
                    tstar = tstarwet
                elif ((t0 <= Tkel) & (t0 >= Tkel-10.)):
                    tstar = abs(t0-Tkel)*constant+tstardry0
                else:
                    tstar = tstardry10
                snowalb = alb_old - (alb_old - albmin)*(tstep/(tstar*secday))

            elif (lalbedo == 3):		# Douville et al, 1995 modified by VdHurk and Viterbo 2003
                if (t0 >= Tkel):
                    snowalb = albmin + (alb_old - albmin) *np.exp(-efoldcoef*tstep/secday)
                else:
                    snowalb = alb_old - reducrate*(tstep/(secday*(1+0.1*(t0-Tkel)**4.)))		# Modified by vdHurk and Viterbo 2003

            
            # Correction albedo for low solar zenith angles, Segal et al, 1991
            if (solzenyes == 1):
                snowalb = snowalb + factor
                if (snowalb > albmax): snowalb = albmax
            
            albmin = albice
            if (firnrest == 1): albmin = albfirn
            albedo = snowalb + (albmin - snowalb) * np.exp(-snowrest/snowstar)
            
        else:		#in case of no precipitation and ice at the surface
            albedo = albice		

    else: 		# in case of fresh snowfall
        albedo = albsnow
        snowdays = 0.

    alb_old = albedo
    if ((solzenyes == 1) & (precip == 0) & (snowrest != 0.)): alb_old = alb_old - factor

    return albedo, snowdays, alb_old

@njit
def solarspectrum (cc,pres,icount,SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear):
    # !===============================================================================
    # USE RADPEN_EBM , ONLY : SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear
    # IMPLICIT NONE
    # !Input
    # INTEGER, INTENT(IN)	:: icount
    # REAL, INTENT(IN)	:: cc, pres

    # !Output
    # REAL	:: solarspectrum

    # ! local
    # REAL :: presplat, pressea
    # REAL :: cloudf,dsdp,SolarawsClear


    presplat = 60900.0		# pressure plateau for standard solar spectrum
    pressea = 98000.0		# presure sea level for standard solar spectrum

    if (SolarPlateauClear[icount] > 0.):
        cloudf = 1.0 - cc * (1.0 - SolarPlateauCloudy[icount]/SolarPlateauClear[icount])
    else:
        cloudf = 0.


    dsdp = (SolarSeaClear[icount] - SolarPlateauClear[icount])/( pressea - presplat)
    SolarawsClear = SolarPlateauClear[icount]+ dsdp * (pres - presplat)
    solarspectrum = SolarawsClear*cloudf

    return solarspectrum

@njit
def cloudcover(temp,lwin):
    # !===============================================================================
    # USE INPUT_EBM, ONLY		: chstation, lwmax, lwmin
    # USE CONSTANTS_EBM, ONLY	: Tkel

    # IMPLICIT NONE

    # !Input
    # REAL,INTENT(IN) :: temp, lwin
    # REAL :: cloudcover
    # !local
    # REAL :: max, min, reltemp

    reltemp = temp - Tkel

    max = lwmax[0] + lwmax[1]*reltemp + lwmax[2]*reltemp**2.
    min = lwmin[0] + lwmin[1]*reltemp + lwmin[2]*reltemp**2.

    cloudcover = (lwin - min)/(max - min)
    if (cloudcover > 1.): cloudcover = 1.
    if (cloudcover < 0.): cloudcover = 0.

    return cloudcover


@njit
def RadiationPenetration(z,dz,grainsize,dens,sbuf,Snet,lid,cloud,
        Lambda, SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear, 
        lambdaAll, dlambdaAll, asymAll, qextAll, cosinglescatAll,
        dlambda,asymsn, qextsn, cosinglescatsn,
        asymice, qextice, cosinglescatice):
    # !===============================================================================
    # ! This routine calculates the fluxes of shortwave radiation that penetrate into
    # ! the snow layer following Brand and Warren J. Glac. 1993
    # ! calculations are done in terms of the radiation flux divergence dqdz
    # !===============================================================================
    # USE CONSTANTS_EBM, ONLY: nlmax
    # USE GLOBALS_EBM, ONLY: cloud, sbuf, Snet
    # USE INPUT_EBM, ONLY: densice, dzrad, lcomment, lwcloud, radiusice, radiussn, zradmax, lalbedo
    # USE RADPEN_EBM, ONLY:	asymice, asymsn, bandmax, cosinglescatice, cosinglescatsn, dlambda,&
    # &						dsdz, qextice, qextsn, sumdivs,&
    # &						asymAll, cosinglescatAll, qextAll, dlambdaAll
    # USE SNOW_EBM, ONLY: dens, dz, lid, nlsnow, z, grainsize
    # USE FILES, ONLY: uo1

    # IMPLICIT NONE

    # !Function
    # REAL :: solarspectrum

    # !Local
    # INTEGER :: i, il, k
    # INTEGER :: kmax, kmin, nlradnd, nlradst
    # REAL :: aext, radius, rcoef, rhorad, scaling, sumdo, sumup, sz0, sigext
    # REAL, DIMENSION(bandmax) :: asym, cosinglescat, klam, mlam, qext, solarscaled
    # REAL, DIMENSION(nlmax)  :: sz

    # REAL :: aext1, aext2, sigext1, sigext2, rcoef1, rcoef2, sz01, sz02, sumup1, sumup2,&
    # & sumdo1, sumdo2, sz1, sz2, scaling1, scaling2
    # REAL, DIMENSION(bandmax) :: asym1, asym2, cosinglescat1, cosinglescat2, dlambda1, dlambda2,&
    # & klam1, klam2, mlam1, mlam2, qext1, qext2
    # REAL :: fraction

    # INTEGER :: igrainbnds
    # REAL , DIMENSION(7) :: grainbnds

    grainbnds = [0.00005,0.0001,0.0002,0.00035,0.0005,0.001,0.0025]

    sz = np.zeros(nlmax)
    dsdz = np.zeros(nlmax)
    sumdivs = 0.
 
    klam = np.zeros(bandmax)
    mlam = np.zeros(bandmax)
    solarscaled = np.zeros(bandmax)
    klam1 = np.zeros(bandmax)
    klam2 = np.zeros(bandmax)
    mlam1 = np.zeros(bandmax)
    mlam2 = np.zeros(bandmax)

    
    if((lalbedo == 4) | (lalbedo == 5)):	# Kuipers Munneke albedo scheme with snow grain size calculation
        for il in range(0,nlmax):
            if((z[il]+0.5*dz[il]) > zradmax):
                nlradnd = il
                break
                
        for il in range(0,nlradnd):
            igrainbnds = 7
            for i in range(6,0,-1):
                if (grainsize[il] <= grainbnds[i]): igrainbnds = i
            if (grainsize[il] == grainbnds[6]): igrainbnds = 7

            if ((igrainbnds == 7) | (igrainbnds == 0)):
                if (igrainbnds == 7): igrainbnds = 6
                dlambda = dlambdaAll[igrainbnds,:]
                asym = asymAll[igrainbnds,:]
                cosinglescat = cosinglescatAll[igrainbnds,:]
                qext = qextAll[igrainbnds,:]

                for i in range(0,bandmax):
                    sigext = (qext[i]*3.*dens[il])/(4.*densice*grainsize[il])
                    aext = sigext * cosinglescat[i]
                    rcoef = 0.5*sigext*(1.0 - asym[i])*(1.0 - cosinglescat[i])
                    klam[i] = np.sqrt(aext*aext + 2.0*aext*rcoef)
                    mlam[i] = (aext + rcoef - klam[i])/rcoef		# note density dependence falls out of mlam

                
                sz0 = 0.
                for i in range(0, bandmax):
                    if (lwcloud != 1):
                        solarscaled[i] = (sbuf[6]/413.8)*solarspectrum(0.0,sbuf[4],i,SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear)		#413 = �sum of spectrum
                    elif (lwcloud == 1):
                        solarscaled[i] = (sbuf[6]/413.8)*solarspectrum(cloud,sbuf[4],i,SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear)		#413 = �sum of spectrum
                    sz0 = sz0 + 1.0E-3 * dlambda[i] * solarscaled[i] * (1.0 - mlam[i])
                
                scaling = (Snet-sumdivs)/sz0
                
                if (il == 0):
                    sumup = 0.0
                    for i in range(0,bandmax):
                        sumup = sumup + 1.0E-3 * dlambda[i] * solarscaled[i] * (1.0-mlam[i]) * np.exp(-klam[i] * dzrad)

                sumdo = 0.0
                for i in range(0,bandmax):
                    sumdo = sumdo + 1.0E-3 * dlambda[i] * solarscaled[i] * (1.0-mlam[i]) * np.exp(-klam[i] * (z[i]+0.5*dz[i]))

                sz[i] = -abs(scaling*(sumup-sumdo))
                sumup = sumdo
                dsdz[i] = sz[i]/dz[i]
                sumdivs = sumdivs - sz[i]
                if ( sumdo*scaling < 0.01 ):
                    break # Go To 200
            else:

                dlambda1 = dlambdaAll[igrainbnds,:]
                asym1 = asymAll[igrainbnds,:]
                cosinglescat1 = cosinglescatAll[igrainbnds,:]
                qext1 = qextAll[igrainbnds,:]
                
                dlambda2 = dlambdaAll[igrainbnds-1,:]
                asym2 = asymAll[igrainbnds-1,:]
                cosinglescat2 = cosinglescatAll[igrainbnds-1,:]
                qext2 = qextAll[5,:]
                
                fraction = (grainsize[il]-grainbnds[igrainbnds-1])/(grainbnds[igrainbnds]-grainbnds[igrainbnds-1])
                
                for i in range(0,bandmax):
                    sigext1 = (qext1[i]*3.*dens[il])/(4.*densice*grainsize[il])
                    aext1 = sigext1 * cosinglescat1[i]
                    rcoef1 = 0.5*sigext1*(1.0 - asym1[i])*(1.0 - cosinglescat1[i])
                    klam1[i] = np.sqrt(aext1*aext1 + 2.0*aext1*rcoef1)
                    mlam1[i] = (aext1 + rcoef1 - klam1[i])/rcoef1		# note density dependence falls out of mlam
                    
                    sigext2 = (qext2[i]*3.*dens[il])/(4.*densice*grainsize[il])
                    aext2 = sigext2 * cosinglescat2[i]
                    rcoef2 = 0.5*sigext2*(1.0 - asym2[i])*(1.0 - cosinglescat2[i])
                    klam2[i] = np.sqrt(aext2*aext2 + 2.0*aext2*rcoef2)
                    mlam2[i] = (aext2 + rcoef2 - klam2[i])/rcoef2		# note density dependence falls out of mlam

                    
                sz01 = 0.
                sz02 = 0.
                for i in range(0, bandmax):
                    if (lwcloud != 1):
                        solarscaled[i] = (sbuf[6]/413.8)*solarspectrum(0.0,sbuf[4],i,SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear)		#413 = �sum of spectrum
                    elif (lwcloud == 1):
                        solarscaled[i] = (sbuf[6]/413.8)*solarspectrum(cloud,sbuf[4],i,SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear)		#413 = �sum of spectrum
                    sz01 = sz01 + 1.0E-3 * dlambda1[i] * solarscaled[i] * (1.0 - mlam1[i])
                    sz02 = sz02 + 1.0E-3 * dlambda2[i] * solarscaled[i] * (1.0 - mlam2[i])
                
                scaling1 = (Snet-sumdivs)/sz01
                scaling2 = (Snet-sumdivs)/sz02
                
                if (il == 0):
                    sumup1 = 0.0
                    sumup2 = 0.0
                    for i in range(0,bandmax):
                        sumup1 = sumup1 + 1.0E-3 * dlambda1[i] * solarscaled[i] * (1.0-mlam1[i]) * np.exp(-klam1[i] * dzrad)
                        sumup2 = sumup2 + 1.0E-3 * dlambda2[i] * solarscaled[i] * (1.0-mlam2[i]) * np.exp(-klam2[i] * dzrad)

                sumdo1 = 0.0
                sumdo2 = 0.0
                for i in range(0,bandmax):
                    sumdo1 = sumdo1 + 1.0E-3 * dlambda1[i] * solarscaled[i] * (1.0-mlam1[i]) * np.exp(-klam1[i] * (z[il]+0.5*dz[il]))
                    sumdo2 = sumdo2 + 1.0E-3 * dlambda2[i] * solarscaled[i] * (1.0-mlam2[i]) * np.exp(-klam2[i] * (z[il]+0.5*dz[il]))

                sz1 = abs(scaling1*(sumup1-sumdo1))
                sz2 = abs(scaling2*(sumup2-sumdo2))
                sz[il] = -abs(min(sz1,sz2) + abs(sz1-sz2)*fraction)
                
                sumup1 = sumdo1
                sumup2 = sumdo2
                dsdz[il] = sz[il]/dz[il]
                sumdivs = sumdivs - sz[il]
                if ( (sumdo1*scaling1 < 0.01) | (sumdo2*scaling2 < 0.01) ): 
                    break
            
    # 200 CONTINUE
    else:  
        for il in range (0,nlmax):
            if((z[il]+0.5*dz[il])>zradmax):
                nlradnd = il
                break
                # GOTO 300
        # 300 CONTINUE
        
        for il in range(0,nlradnd):
            # Here read snow/ice values into the appropriate arrays
            if (lid[il] == 0):
                radius = radiusice
                asym = asymice
                qext = qextice
                cosinglescat = cosinglescatice
            else:
                radius = radiussn
                asym = asymsn
                qext = qextsn
                cosinglescat = cosinglescatsn
            
            for i in range(0,bandmax):
                sigext = (qext[i]*3.*dens[il])/(4.*densice*radius)
                aext = sigext * cosinglescat[i]
                rcoef = 0.5*sigext*(1.0 - asym[i])*(1.0 - cosinglescat[i])
                klam[i] = np.sqrt(aext*aext + 2.0*aext*rcoef)
                mlam[i] = (aext + rcoef - klam[i])/rcoef		#! note density dependence falls out of mlam

                
            sz0 = 0.
            for i in range(0, bandmax):
                if (lwcloud != 1):
                    solarscaled[i] = (sbuf[6]/413.8)*solarspectrum(0.0,sbuf[4],i,SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear)		#413 = �sum of spectrum
                elif (lwcloud == 1):
                    solarscaled[i] = (sbuf[6]/413.8)*solarspectrum(cloud,sbuf[4],i,SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear)		#413 = �sum of spectrum
                sz0 = sz0 + 1.0E-3 * dlambda[i] * solarscaled[i] * (1.0 - mlam[i])
            
            scaling = (Snet-sumdivs)/sz0
            if (il == 0):
                sumup = 0.0
                for i in range(0,bandmax):
                    sumup = sumup + 1.0E-3 * dlambda[i] * solarscaled[i] * (1.0-mlam[i]) * np.exp(-klam[i] * dzrad)

            sumdo = 0.0
            for i in range(0,bandmax):
                sumdo = sumdo + 1.0E-3 * dlambda[i] * solarscaled[i] * (1.0-mlam[i]) * np.exp(-klam[i] * (z[il]+0.5*dz[il]))

            sz[il] = -abs(scaling*(sumup-sumdo))
            sumup = sumdo
            dsdz[il] = sz[il]/dz[il]
            sumdivs = sumdivs - sz[il]
            if ( sumdo*scaling < 0.01 ): break #GOTO 400
        # k=kmin,kmax
        # 400 CONTINUE
        
    if(((il-2) == nlradnd)):
        print("Warning: shortwave radiation penetration has reached zradmax, remaining energy: ", sumdo*scaling)
        if(lcomment == 1):
            print("Warning: shortwave radiation penetration has reached zradmax, remaining energy: ", sumdo*scaling)

    return dsdz, sumdivs