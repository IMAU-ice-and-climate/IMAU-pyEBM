import numpy as np
from numba import njit

from globals import *
from info import *
from routines import densprec


# !===============================================================================
# !
# ! File with routines and functions related to using sonic altimeter data to limit
# ! the accumulation to observed values 
# ! and use it for evaluation of the modelled snow and ice melt. 
# !
# !===============================================================================
@njit
def getaccumandmelt(buf, nt, racc_old, acclevel):
# !===============================================================================
# ! Determines accumulation on ice and ice melt based in input sonic altimeter time series (17)
# ! and puts them in columns 19 = accum, 20 = melt , 21 = drift
# ! 1=year, 2=day, 3=hour, 4=T, 5=p, 6=WS, 7=Sin, 8=Sout, 9=albedo, 
# ! 10=qi, 11=Lin, 12=Lout, 13=z0m, 14=zt, 15=zm , 16=serie(running mean), 17=precip,  
# ! 18=accum, 19=melt , 20=drift 
# !===============================================================================
# USE GLOBALS_EBM , ONLY : ilast , buf 
# USE SNOW_EBM , ONLY : racc_old , acclevel , meltlevel 
# USE INPUT_EBM , ONLY : ibyear
# USE CONSTANTS_EBM , ONLY : errorval 

# IMPLICIT NONE

# ! input
# INTEGER :: iyr
# ! Local
# INTEGER :: ii
# INTEGER :: maxcnt , meltcnt
# REAL :: rserie , racc , rdrift , rsummelt , rmelt

    maxcnt = 30*24		# max number of samples delay with reset of meltlevel, prevents jumping behaviour		
    meltcnt = 0

    # if (iyr == ibyear):
    #     racc_old = buf(1,16)
    #     acclevel = 0		# raccm1
    #     meltlevel = 0		# rmelt1

    for ii in range(0,nt):	#  per calendar year
    
        rserie = buf[ii,15]
        
        if (rserie > errorval):
            if (rserie > acclevel):
                rsummelt = acclevel
                racc = rserie - acclevel      
                meltcnt = meltcnt + 1
                if (meltcnt > maxcnt): meltlevel=acclevel
            else:
                rsummelt = rserie
                acclevel = rserie
                racc = 0.
                meltcnt = 0

            rdrift = min(0.,(racc - racc_old))
            racc_old = racc
        #    rmelt = rsummelt-meltlevel

        buf[ii,17] = racc
        buf[ii,18] = rsummelt		# rmelt
        buf[ii,19] = rdrift


    return buf, racc_old, acclevel

@njit
def climaccuandmelt(buf, nt, climracc_old, climserie, climacclevel, climmeltlevel):
# !===============================================================================
# ! Determines accumulation on ice and ice melt based in input sonic altimeter time series (17)
# ! and puts them in columns 19 = accum, 20 = melt , 21 = drift
# ! 1=year, 2=day, 3=hour, 4=T, 5=p, 6=WS, 7=Sin, 8=Sout, 9=albedo, 
# ! 10=qi, 11=Lin, 12=Lout, 13=z0m, 14=zt, 15=zm , 16=serie(running mean), 17=precip,  
# ! 18=accum, 19=melt , 20=drift 
# !===============================================================================
    # USE GLOBALS_EBM , ONLY : ilast , buf 
    # USE SNOW_EBM , ONLY : rhosn , climserie , climracc_old , climacclevel , climmeltlevel 
    # USE INPUT_EBM , ONLY : ibyear, lclimtemp , trhoprec , rhosnprec , lclimprec , climprec
    # USE CONSTANTS_EBM , ONLY : errorval , Tkel 

    # IMPLICIT NONE

    # ! input
    # INTEGER :: iyr
    # ! Local
    # INTEGER :: ii
    # INTEGER :: maxcnt , meltcnt
    # REAL :: rserie , racc , rdrift , rsummelt , rmelt , rprecip
    # REAL :: dens , densprec
    # REAL :: sum

    maxcnt = 30*24		# max number of samples delay with reset of meltlevel, prevents jumping behaviour		
    meltcnt = 0

    # if (iyr == ibyear):
    # climracc_old = buf(1,15)
    # climserie = buf(1,15)
    # climacclevel = 0		# raccm1
    # climmeltlevel = 0		# rmelt1


    sum = climserie

    for ii  in range(0,nt):
    # !  rserie = buf(ii,16)
    # !  racc = buf(ii,18)
        rsummelt = buf[ii,18]
        rdrift = buf[ii,19]

        # modify serie, acc and drift in case climate sensitivity tests depending on temperature
        if (trhoprec == 0):		# constant value set in input
            dens = rhosnprec
        elif (trhoprec == 1):	# value depends on wind speed and surface temperature
            dens = densprec(buf[ii,5],buf[ii,3])
        rprecip = buf[ii,16]/dens		# 16 = precip

    # !  IF (lclimtemp == 1) THEN 
        if (buf[ii,3] >= Tkel + 2.5):		# precip turned into water
            rprecip = 0.
        elif (buf[ii,3] > Tkel + 0.5):
            rprecip = 0.5*(Tkel + 2.5 - buf[ii,3])*rprecip
        
        rprecip = rprecip*(1.+lclimprec*climprec*0.01)
        
        sum = sum + rprecip + rdrift
        climserie = sum + rsummelt
        if (climserie > errorval):
            if (climserie > climacclevel):
                rsummelt = climacclevel
                racc = climserie - climacclevel      
                meltcnt = meltcnt + 1
                if (meltcnt > maxcnt): climmeltlevel=climacclevel
            else:
                rsummelt = climserie
                climacclevel = climserie
                racc = 0.
                meltcnt = 0
            rdrift = min(0.,(racc - climracc_old))
            climracc_old = racc
        #   rmelt = rsummelt-meltlevel

        buf[ii,15] = climserie
        buf[ii,17] = racc
        buf[ii,16] = buf[ii,16]*(1.+lclimprec*climprec*0.01)
        buf[ii,18] = rsummelt		# rmelt
        buf[ii,19] = rdrift

    # !    WRITE(99,*) ii,rserie,racc,rsummelt,rdrift,rprecip,climserie
    
    return buf, climracc_old, climserie, climacclevel, climmeltlevel


@njit
def set_precip(t0,sbuf,water):
# !===============================================================================
# ! Routine that sets height of sensors -> station dependend in some cases !!!!!
# !===============================================================================
    precip_tmp = sbuf[16]
    rhosn = rhosninit
    # determine density of snow fall in order to convert precip in mm we to m snow
    if (tpdens > 0):
        if (trhoprec == 0):		# constant value set in input
            rhosn = rhosnprec
        elif (trhoprec == 1):	# value depends on wind speed and surface temperature
            rhosn = densprec(sbuf[5],t0)		# based on surf temp or 2m temp, or obs temp?
        #    rhosn = densprec(sbuf(6),t2m)
        #    rhosn = densprec(sbuf(6),sbuf(4))
        else:
            rhosn = rhosninit   
        # rhosn set at initialisation when dens profile is set 
    precip = precip_tmp/rhosn # convert to m snow

    # modify precipitation in case climate sensitivity tests depending on temperature
    if (lclimtemp == 1):  
        if (sbuf[3] >= Tkel + 2.5):		# precip turned into water added to first layer water content
            water[0] = water[0] + precip*rhosn
            precip = 0.
        elif ((sbuf[3] > Tkel + 0.5) & (sbuf[3] < Tkel + 2.5)):
            water[0] = water[0] + (1 - 0.5*(Tkel + 2.5 - sbuf[3]))*precip*rhosn    
            precip = 0.5*(Tkel + 2.5 - sbuf[3])*precip


    #precipsum = precipsum + precip*rhosn		!again in mm we
    return precip, water, rhosn

@njit
def set_height(sbuf, t0, dsnowh, dsnowr, rhosn, \
    corrsnow, temp, water, ice, mass, dens, dz, z, lid, nl, nlsnow, melt_l, icemelt,precip,drift):

# !===============================================================================
# ! Routine that corrects accumulation height in case acc restricted to obs !!!!!
# !===============================================================================
# USE GLOBALS_EBM , ONLY : t0 , sbuf
# USE SNOW_EBM , ONLY : acc , drift , dsnowh , dsnowr , precip , hsnow, corrsnow , &
# &                     rhosn, icemelt, melt_l, & !, precipsum
# &                     temp , water , mass , ice , dens , dz , z , lid , nl , nlsnow 
# USE INPUT_EBM , ONLY : tstep , dz0 , luseacc , lcomment
# USE CONSTANTS_EBM , ONLY : Tkel

    # IMPLICIT NONE

    # !local
    # INTEGER :: il , illost
    # REAL :: hu , llprec
    # REAL :: obsmelt
    # REAL :: corrdrift , corrmelt 

    # ! sbuf(16) = (running average) sonic height measurement in serie
    # ! sbuf(18) = (running average) acc measurement determined from sonic height
    # ! dsnow = snow thickness at start, equals series sonic height ranger at start
    # ! dsnowh = acc as moddeled
    # ! acc = acc observed
    # ! hsnow = sonic series observations, set in interp_data
    # ! hsnowmod_l = snow depth previous time step model

    hu = 0.01 #! estimated uncertainty in 1 sonic height observations
    llprec = hu*tstep/(24*3600.)  # equals hu in case tstep is 1 hour

    acc = sbuf[17]
    obsmelt = sbuf[18]
    #drift = sbuf[19] #	already defined in interp_data

    if (luseacc == 2): # accumulation is restricted to observed accumulation
        if ((dsnowh + dsnowr) - acc >= 0.):		# current model acc > observed acc, no use to add anymore snow
            if (t0 < 270.):
                drift = acc - (dsnowh + dsnowr)	# and later on remove the too much snow, assuming it to be drift
                precip = 0.     				# 16 = precip in mm we. 
            else:
                drift = 0.
                precip = 0.
        elif ((dsnowh + dsnowr) - acc < 0.):	# current model acc < observed acc, can add a bit more
            drift = acc - (dsnowh + dsnowr)		# tmp param > 0 because observed acc > modelled
            if (precip >= drift):
                precip = drift
                drift = 0.   # drift - precip
            else:
                drift = 0.
    
    # !! IF (dsnowh - acc > 0.01*dz0) THEN	!in case measured acc < modelled acc
    # ! IF ((dsnowh + dsnowr) - acc >= 0.) THEN		! current model acc > observed acc, not allowed to add anymore snow
    # !  IF (precip > 0.) THEN
    # !    drift = acc - (dsnowh + dsnowr)				! and later on remove the too much snow, assuming it to be drift
    # !    precip = 0.     							! 17 = precip in mm we. 
    # !  ELSE IF (drift < -2.*llprec .and. t0 < 270.) THEN	!Only adjust in case there is no melt and the initial drift islarge
    # !    drift = -0.5*(dsnowh-acc)
    # !  ELSE 
    # !    drift = 0.
    # !  ENDIF
    # ! ELSE
    # !  drift = 0.
    # ! ENDIF
    
    elif (luseacc == 3): # the snow accumulation is completely fixed by observed sonic ranger data
        if ((dsnowh + dsnowr) - acc >= 0.):		# current model acc > observed acc, no use to add anymore snow
            drift = acc - (dsnowh + dsnowr)				# and later on in snowmodel remove the too much snow, assuming it to be drift
            precip = 0.     							# 16 = precip in mm we. 
        elif ((dsnowh + dsnowr) - acc < 0.):	# current model acc < observed acc, can add a bit more
            drift = acc - (dsnowh + dsnowr)		# tmp param > 0 because observed acc > modelled
            if (precip >= drift):
                precip = drift
                drift = drift - precip
            elif((precip < drift) & (drift > 2.*llprec)):
                drift = 0.5*(acc - dsnowh)
    # !    IF (precip >= drift) THEN     
    # !      precip = drift
    # !      drift = 0
    # !    ELSE 
    # !      drift = 0.
    # !    ENDIF     

    # CORRECTION of hsnowmod after period with no data 
    if (corrsnow != 0.):
        corrdrift = dsnowh-acc
        corrmelt = obsmelt - melt_l		# > 0 then ice melt has occured
        if (corrdrift < 0.): precip = precip - corrdrift
        
        if (corrmelt < 0.): icemelt = icemelt + corrmelt
        if (corrsnow > 0.): dsnowh = dsnowh + corrsnow 

        corrdrift = -corrdrift*rhosn
        il = 0
        while (corrdrift < 0.):
            mass[il] = mass[il]+ corrdrift
            if (mass[il] < 0.):
                corrdrift = mass[il]
                il = il+1
            else:
                corrdrift = 0.
        if (il > 0):
            illost = il-1
            if (lcomment == 1): print(il-1,' Layer(s) removed to correct period with no observations')
            for il in range(0,nl-illost):
                temp[il] = temp[il+illost]
                water[il] = water[il+illost]
                ice[il] = ice[il+illost]
                mass[il] = mass[il+illost]
                dens[il] = dens[il+illost]
                dz[il] = dz[il+illost]
                if (il == 1): z[il] = 0.5*dz[il]
                if (il > 1): z[il] = 0.5*dz[il] + z[il-1]
                lid[il] = lid[il+illost]
            nl = nl-illost
            nlsnow = nlsnow - illost
            if (nlsnow < 0): nlsnow = 0
        corrsnow = 0.

    melt_l = obsmelt

    return temp, water, ice, mass, dens, dz, z, lid, nl, nlsnow, corrsnow, precip, drift, icemelt, dsnowh, melt_l