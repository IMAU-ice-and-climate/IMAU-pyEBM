#!/usr/bin/env python3
''' Standard library imports '''
import numpy as np
import pandas as pd
from numba import njit

''' Local imports '''
from info import * 
from routines import *
from globals import *
from initial import *
from massbalance import getaccumandmelt, climaccuandmelt

def read_info(file):
    """
    Read "awsid_ebm.txt" info file in previously used FORTRAN format
    """
    with open(file) as f:
        lines = f.readlines().replace(",","")

    lcomment = int(lines[0].split()[0])
    lerrorgap = int(lines[1].split()[0])
    lmc = int(lines[2].split()[0])
    lhourlysnowoutput = int(lines[3].split()[0])
    ibyear = int(lines[4].split()[0])
    ilyear = int(lines[5].split()[0])
    tstep = int(lines[6].split()[0])
    mbsumdy, mbwindy= float(lines[7].split()[0]), float(lines[7].split()[1])
    dz0 = float(lines[8].split()[0])
    dzdeep = float(lines[9].split()[0])
    zdeep = float(lines[10].split()[0])
    densice, densfirn= float(lines[11].split()[0]), float(lines[11].split()[1])
    rhosn, rhosnprec, trhoprec = float(lines[12].split()[0]), float(lines[12].split()[1]), int(lines[12].split()[2])
    tpcond = int(lines[13].split()[0])
    T10m = float(lines[14].split()[0])
    dsnow, dfirn  = float(lines[15].split()[0]), float(lines[15].split()[1])
    luseacc = int(lines[16].split()[0])
    tpdens, lrefr = int(lines[17].split()[0]), int(lines[17].split()[1])
    tpirre, cirre = int(lines[18].split()[0]), float(lines[18].split()[1])
    lslush = int(lines[19].split()[0])
    surfangle, tausteep, tauhor,tau1,slfact = float(lines[20].split()[0]), float(lines[20].split()[1]), float(lines[20].split()[2]), float(lines[20].split()[3]), float(lines[20].split()[4])
    accyear = float(lines[21].split()[0])
    lz0m, zll, zul, Hmax = float(lines[22].split()[0]), float(lines[22].split()[1]), float(lines[22].split()[2]), float(lines[22].split()[3])
    z0msn, z0mice = float(lines[23].split()[0]), float(lines[23].split()[1])
    lz0h = int(lines[24].split()[0])
    tcalc = int(lines[25].split()[0])
    extrapolation = int(lines[26].split()[0])
    lsnet = int(lines[27].split()[0])
    albmin, albmax = float(lines[28].split()[0]), float(lines[28].split()[1])
    emis = float(lines[29].split()[0])
    lalbedo, solzenyes, SSAfresh, radrefr = int(lines[30].split()[0]), int(lines[30].split()[1]), float(lines[30].split()[2]), float(lines[30].split()[3])
    albsnow, albice, albfirn, soot = float(lines[31].split()[0]), float(lines[31].split()[1]), float(lines[31].split()[2]), float(lines[31].split()[3])
    snowstar,tstarwet,tstardry0,tstardry10 = float(lines[32].split()[0]), float(lines[32].split()[1]), float(lines[32].split()[2]), float(lines[32].split()[3])
    penetration = int(lines[33].split()[0])
    dzrad,zradmax = float(lines[34].split()[0].replace(",","")), float(lines[34].split()[1])
    radiussn, radiusice = float(lines[35].split()[0]), float(lines[35].split()[1])
    lwcloud = int(lines[36].split()[0])
    lwmax = float(lines[37].split()[0]), float(lines[37].split()[1]), float(lines[37].split()[2])
    lwmin = float(lines[38].split()[0]), float(lines[38].split()[1]), float(lines[38].split()[2])
    depthin = float(lines[39].split()[0]), float(lines[39].split()[1]), float(lines[39].split()[2]), float(lines[39].split()[3]), float(lines[39].split()[4])
    lclimtemp,climtemp = int(lines[40].split()[0]), float(lines[40].split()[1])
    lclimprec,climprec = int(lines[41].split()[0]), float(lines[41].split()[1])
    lclimrad,climrad = int(lines[42].split()[0]), float(lines[42].split()[1])
    lclimws,climws = int(lines[43].split()[0]), float(lines[43].split()[1])
    

def read_input(ifile):
    """
    Read "awsid_HOUR-EBM.txt" forcing file in previously used FORTRAN format
    """
    df = pd.read_csv(ifile,delim_whitespace=True,header=0)

    date = pd.to_datetime(df['Date'], format='%Y-%m-%d')
    hour = pd.to_datetime(df['Hour'], format='%H:%M:%S')
     
    ilasttot =  len(df)
    iilist = np.where(date.dt.year.values == iyr)
    ilast = len(df)
    
    # Save input to array ibuf
    ibuf[iilist,0] = date.dt.year.values[iilist] # Year
    ibuf[iilist,2] = hour.dt.hour.values[iilist]+hour.dt.min.values[iilist]/60  # Hour
    for ii in iilist:
        ibuf[ii,1] = julday(ibuf[ii,0],date.dt.month.values[ii],date.dt.day.values[ii])+ibuf[ii,2]/24    # Julian day plus time
    ibuf[iilist,3] = df.T.values[iilist] # Temperature
    ibuf[iilist,4] = df.P.values[iilist] # Pressure
    ibuf[iilist,5] = df.WS.values[iilist] # Wind speed
    ibuf[iilist,6] = df.Sin.values[iilist] # Shortwave incoming radiation
    ibuf[iilist,7] = df.Sout.values[iilist] # Shortwave outgoing radiation
    ibuf[iilist,8] = df.alb.values[iilist] # Albedo
    ibuf[iilist,9] = df.q.values[iilist] # Specific humidity
    ibuf[iilist,10] = df.Lin.values[iilist] # Longwave incoming radiation
    ibuf[iilist,11] = df.Lout.values[iilist] # Longwave outgoing radiation
    ibuf[iilist,12] = df.z0m.values[iilist] # z0m = roughness length surface combi of snow and ice values.
    ibuf[iilist,13] = df.zt.values[iilist] # zt = height temperature sensor above surface
    ibuf[iilist,14] = df.zm.values[iilist] # zm = height wind speed sensor above surface
    ibuf[iilist,15] = df.Serie.values[iilist] # sonic height obs in 1 series, preferably running mean height, at t = 0 equal dsnow, dsnow must be > 0
    ibuf[iilist,16] = df.precip.values[iilist] # Precipitation
    zenith[iilist] = df.zenith.values[iilist] # Solar zenith angle

    # Check input
    # Minimum wind speed
    ibuf[:,5][ibuf[:,5] < 0.1] = 0.1
  
    # Minimum pressure
    ibuf[:,4][ibuf[:,4] < 0] = 0
  
    # Minimum reflected shortwave radiation
    ibuf[:,7][ibuf[:,7] < 0] = 0
  
    # Incoming shortwave radiation should be between Sout/albmin and Sout/albmax
    ibuf[:,6][ibuf[:,6] > ibuf[:,7]/albmax] = ibuf[:,7]/albmax
    ibuf[:,6][ibuf[:,6] < ibuf[:,7]/albmin] = ibuf[:,7]/albmin
  
    # Minimum longwave radiation
    ibuf[:,10][ibuf[:,10] < 0] = 0
    ibuf[:,11][ibuf[:,11] < 0] = 0 
    
    return ibuf

@njit
def get_errorflag(awsid, paramerror, nt):
    errorflag = np.full(nt+1, 0.0)
    for ii in range(0,nt):
        if (awsid[ii] == 0): errorflag[ii] = maxerror
        # in this case errorflag will/should be 61

        # First impact of rime = 1: negative
        if (paramerror[ii] < 0): errorflag[ii] = errorflag[ii] + 1.
        paramerror[ii] = abs(paramerror[ii])
        
        # Second impact of temperature = 16: f = temperature = 6th column
        ip = 6
        dum = int(paramerror[ii] /(10**(ip-1))) - 10*int(paramerror[ii] /(10**(ip)))
        if (round(dum) >= 1): errorflag[ii] = errorflag[ii] + 16.
        
        # Third impact of Short wave in and out, and longwave in = 8: b = Sin = 2 and c = Sout = 3, d = Lin = 4
        # was only Sin and Sout
        ip = 2
        dum = int(paramerror[ii] /(10**(ip-1))) - 10*int(paramerror[ii] /(10**(ip)))
        ip = 3
        dum2 = int(paramerror[ii] /(10**(ip-1))) - 10*int(paramerror[ii] /(10**(ip)))
        ip = 4
        dum3 = int(paramerror[ii] /(10**(ip-1))) - 10*int(paramerror[ii] /(10**(ip)))
        if ((round(dum) >= 1)| (round(dum2) >= 1) | (round(dum3) >= 1)): errorflag[ii] = errorflag[ii] + 8.
        
        # Fourth impact of wind speed = 4: a = wind speed = 1
        ip = 1
        dum = int(paramerror[ii] /(10**(ip-1))) - 10*int(paramerror[ii] /(10**(ip)))
        if (round(dum) >= 1): errorflag[ii] = errorflag[ii] + 4.
        
        # Last impact of air pressure and humidity = 2: g = RH = 7 h = P = 8
        # was Lin
        ip = 7
        dum = int(paramerror[ii] /(10**(ip-1))) - 10*int(paramerror[ii] /(10**(ip)))
        ip = 8
        dum2 = int(paramerror[ii] /(10**(ip-1))) - 10*int(paramerror[ii] /(10**(ip)))
        if ((round(dum) >= 1) | (round(dum2) >= 1)): errorflag[ii] = errorflag[ii] + 2.
        # all add up max is 30 + 1, if 30, than no data present
    errorflag[nt+1] = errorflag[nt]    
    return errorflag

@njit
def checkdata(ibuf,nt,racc_old, acclevel, climracc_old, climserie, climacclevel, climmeltlevel):
# !===============================================================================
# ! Routine that checks the input data and corrects when necessary
# ! linearly interpolates for missing values
# ! 1=year, 2=day, 3=hour, 4=T, 5=p, 6=WS, 7=Sin, 8=Sout, 9=albedo, 
# ! 10=qi, 11=Lin, 12=Lout, 13=z0m, 14=zt, 15=zm , 16=serie(running mean), 17=precip,  
# ! 18=accum, 19=melt , 20=drift determined from 16
# !===============================================================================
# USE GLOBALS_EBM , ONLY : ilast, ibuf, buf , alb_old
# USE SNOW_EBM , ONLY : hsnowmod, hsnowmod_l , dsnowh , hsnowstart , &
# &                     dsnowr , dsnowice , corrsnow
# USE INPUT_EBM , ONLY : dsnow , lcomment, ibyear, luseacc , lclimtemp, climtemp, &
# &                      lclimrad, climrad, lclimprec, climprec, lclimws, climws, &
# &                      albmin , albmax , trhoprec , rhosnprec, chstation
# USE CONSTANTS_EBM , ONLY : colmax, Tkel, StefBoltz , loutmax

    # IMPLICIT NONE

    # ! input
    # INTEGER,INTENT(IN) :: iyr
    # ! Local
    # REAL , PARAMETER :: error = -990.
    # REAL :: rh, relhum, spechum
    # REAL :: emissivity
    # REAL :: dens , densprec
    # INTEGER :: ii,j,iv
    # INTEGER :: precipconvert

    error = -990.
    # buf = -999.
    buf = ibuf

    precipconvert = 1
    if (chstation == "ant_neuma"): precipconvert = 0
    # ! 0 = nothing is done
    # ! 1 = from AWS IMAU data it was first altimeter, then converted to precip in m snow 
    # ! then with given density in convert converted to precip in m we. Here later on 
    # ! converted again to m snow. in case precip from altimeter dens conversion should equal.
    # ! trick here: convert in 'convert' with constant density of 400, than here converted 
    # ! back to m snow and then again to mm we using dens given here. This ensures consistency.

    # check all data points
    for ii in range(0,nt):

        #  check to be sure that every year each value starts with a valid value!
        if ((ii == 0 ) | (ii == nt-1)):
            for iv in range(0,colmax-3):		#16, 17-19 are snow acc, ice melt, drift
                if (ibuf[ii,iv] < error):
                    print('Start or end with error value !!! STOP',ii,iv)
                    # STOP 14


    buf[ii,3] = ibuf[ii,3] 		# temperature in K
    buf[ii,4] = ibuf[ii,4]		# pressure in Pa
    buf[ii,9] = ibuf[ii,9]	# Specific humidity in kg/kg
    
    # set lower limit of wind speed
    if (buf[ii,5] < 0.1): buf[ii,5] = 0.1

    # set upper and lower limit albedo  
    if (buf[ii,8] < albmin): buf[ii,8] = albmin
    if (buf[ii,8] > albmax): buf[ii,8] = albmax
    
    # set upper and lower limit roughness  
    if (buf[ii,12] < 1e-4): buf[ii,12] = 1e-4
    if (buf[ii,12] > 2): buf[ii,12] = 2
    
    if (precipconvert == 1):
        dens = 400.
        buf[ii,16] = buf[ii,16]/dens
        if (trhoprec == 0):	# constant value set in input
            dens = rhosnprec
        elif (trhoprec == 1): # value depends on wind speed and surface temperature
            dens = densprec(buf[ii,5],buf[ii,3])
        buf[ii,16] = buf[ii,16]*dens

    if ( (lclimtemp == 1) & (climtemp != 0.) ):
        rh = relhum(buf[ii,3],buf[ii,9],buf[ii,4]) 
        emissivity = buf[ii,10]/(StefBoltz*(buf[ii,3]**4))
        buf[ii,3] = buf[ii,3] + climtemp		# climate sensitivity
        buf[ii,9] = spechum(buf[ii,3],rh,buf[ii,4])		# change spec hum assuming rel hum remains constant
        buf[ii,10] = emissivity*StefBoltz*(buf[ii,3]**4)	# change Lin assuming constant emissivity
        
    if ((lclimrad == 1) & (climrad != 0.) ):
        buf[ii,6] = buf[ii,6] * (1.+ climrad*0.01)		# climate sensitivity
        buf[ii,7] = buf[ii,7] * (1.+ climrad*0.01)		# climate sensitivity

    # !  IF (lclimprec == 1 .and. climprec /= 0. ) THEN
    # !!    buf(ii,16) = ibuf(ii,16)*(1.+climprec*0.01)	! serie Climate sensitivity test with precip
    # !    buf(ii,17) = ibuf(ii,17)*(1.+climprec*0.01)	! precip Climate sensitivity test with precip
    # !  ENDIF
    if ((lclimws == 1) & (climws != 0.) ):
        buf[ii,5] = buf[ii,5] * (1.+ climws*0.01)		# climate sensitivity

    # in case serie from sonic height ranger is available then calculate acc, melt and drift
    if (luseacc > 0):
        buf, racc_old, acclevel = \
            getaccumandmelt(buf, nt, racc_old, acclevel)
    if ((lclimtemp == 1) | (lclimprec == 1) | (lclimrad == 1) | (lclimws == 1) ): 
        buf, climracc_old, climserie, climacclevel, climmeltlevel = \
            climaccuandmelt(buf, nt, climracc_old, climserie, climacclevel, climmeltlevel)


    # initialise other necessary parameters
    # if (iyr == ibyear):
    dsnowr = max((buf[0,15]-dsnow_glob),0.)		# initialise rest snow not yet put into a layer
    dsnow = buf[0,15]
    hsnowstart = dsnow
    hsnowmod = dsnow
    hsnowmod_l = dsnow
    dsnowh = dsnow
    corrsnow = 0.
    alb_old = buf[0,8]

    if (lcomment == 1): print('END data checking, start EBM calculation')
    
    return buf, dsnowr, dsnow, hsnowstart, hsnowmod, dsnowh, corrsnow, alb_old, racc_old, acclevel, climracc_old, climserie, climacclevel, climmeltlevel



@njit
def interp_data(buf,zenith,j,i,nstep,lid,ilast):
    """
    Linearly interpolate forcing to model timestep
    """
    sbuf = np.zeros(colmax)
    ivstart = 3
    ivend = colmax
    
    if (j == -1) :	# in case of error data continue to interpolate precip data
        ivstart = 16		# precip
        ivend = 20
        if (luseacc == 0): 
            ivend = 17		

    for iv in range(ivstart,ivend):
        if (i < ilast -1):
            grad = (buf[i+1,iv] - buf[i,iv])/nstep
            sbuf[iv] = buf[i,iv] + j*grad
            if ((iv == 16) | (iv == 19)): sbuf[iv] = buf[i+1,iv] / nstep 		# precip and drift
        elif (i == ilast - 1):
            sbuf[iv] = buf[i,iv]
            if ((iv == 16) | (iv == 19)): sbuf[iv] = buf[i,iv]  / nstep 		# precip and drift

    # Set other parameters per time step
    z0m = sbuf[12]
    zt = sbuf[13]
    zm = sbuf[14]
    # hsnow = depth snow layer, becomes negative after ice melt.
    hsnow = sbuf[15]
    drift = sbuf[19]
    if (luseacc <= 1): drift = 0.

    if ((lz0m == 1) | (lz0m == 2) | (lz0m == 3)):
        z0m = z0msn
        if (luseacc >= 2):
            # z0m implicitly > 0
            if ((lz0m <= 3) & (sbuf[15] < 0.005)):	# in this case observed surface type determines roughness
                z0m = z0mice
            elif ((lz0m == 10) & (lid == 0)):	# in this case modelled surface type determines roughness
                z0m = z0mice
        elif (lid == 0):		# of lack of more information, must be based on model surface type
            z0m = z0mice
    # else value for z0m as a function of time is given in the input file, or parameterized

    
    if (i == ilast-1): 
        grad = zenith[i]
    else: 
        grad = (zenith[i+1] - zenith[i])/nstep
    szenith = zenith[i] + j*grad
    
    sbuf[12] = z0m
    sbuf[13] = zt
    sbuf[14] = zm
    
    return sbuf, szenith, hsnow, drift

def inputradpen():
# !===============================================================================
# ! Routine that reads input information for the radiation penetration routine
# !===============================================================================
# USE RADPEN_EBM , ONLY : bandmax, lambda, SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear, &
# &                       dlambda, asymsn, qextsn, cosinglescatsn, asymice, qextice, cosinglescatice,&
# &						asymAll, cosinglescatAll, qextAll, dlambdaAll, lambdaAll
# USE INPUT_EBM , ONLY : radiussn,radiusice, lalbedo
# USE FILES , ONLY : ui

# IMPLICIT NONE

# INTEGER :: k	
# INTEGER :: ii
# REAL :: radius
    
    
    file_1 = 'IN_SolarSpectrum.txt'
    
    tmp = np.loadtxt(file_1,skiprows=1)
    Lambda = tmp[:,0]
    SolarPlateauClear = tmp[:,1]
    SolarPlateauCloudy = tmp[:,2]
    SolarSeaClear = tmp[:,3]
    
    lambdaAll = np.zeros((7,bandmax))
    dlambdaAll = np.zeros((7,bandmax))
    asymAll = np.zeros((7,bandmax))
    qextAll = np.zeros((7,bandmax))
    cosinglescatAll = np.zeros((7,bandmax))
    if ((lalbedo == 4) | (lalbedo == 5)):
        file_2 = 'IN_Mie(r=0.05mm).txt'
        tmp = np.loadtxt(file_2,skiprows=0)
        lambdaAll[0,:] = tmp[:,1]
        dlambdaAll[0,:] = tmp[:,2]
        asymAll[0,:] = tmp[:,3]
        qextAll[0,:] = tmp[:,4]
        cosinglescatAll[0,:] = tmp[:,5]
        
        file_3 = 'IN_Mie(r=0.1mm).txt'
        tmp = np.loadtxt(file_3,skiprows=0)
        lambdaAll[1,:] = tmp[:,1]
        dlambdaAll[1,:] = tmp[:,2]
        asymAll[1,:] = tmp[:,3]
        qextAll[1,:] = tmp[:,4]
        cosinglescatAll[1,:] = tmp[:,5]

        file_4 = 'IN_Mie(r=0.2mm).txt'
        tmp = np.loadtxt(file_4,skiprows=0)
        lambdaAll[2,:] = tmp[:,1]
        dlambdaAll[2,:] = tmp[:,2]
        asymAll[2,:] = tmp[:,3]
        qextAll[2,:] = tmp[:,4]
        cosinglescatAll[2,:] = tmp[:,5]
        
        file_5 = 'IN_Mie(r=0.35mm).txt'
        tmp = np.loadtxt(file_5,skiprows=0)
        lambdaAll[3,:] = tmp[:,1]
        dlambdaAll[3,:] = tmp[:,2]
        asymAll[3,:] = tmp[:,3]
        qextAll[3,:] = tmp[:,4]
        cosinglescatAll[3,:] = tmp[:,5]
        
        file_6 = 'IN_Mie(r=0.5mm).txt'
        tmp = np.loadtxt(file_6,skiprows=0)
        lambdaAll[4,:] = tmp[:,1]
        dlambdaAll[4,:] = tmp[:,2]
        asymAll[4,:] = tmp[:,3]
        qextAll[4,:] = tmp[:,4]
        cosinglescatAll[4,:] = tmp[:,5]
        
        file_7 = 'IN_Mie(r=1.0mm).txt'
        tmp = np.loadtxt(file_7,skiprows=0)
        lambdaAll[5,:] = tmp[:,1]
        dlambdaAll[5,:] = tmp[:,2]
        asymAll[5,:] = tmp[:,3]
        qextAll[5,:] = tmp[:,4]
        cosinglescatAll[5,:] = tmp[:,5]
        
        file_8 = 'IN_Mie(r=2.5mm).txt'
        tmp = np.loadtxt(file_8,skiprows=0)
        lambdaAll[6,:] = tmp[:,1]
        dlambdaAll[6,:] = tmp[:,2]
        asymAll[6,:] = tmp[:,3]
        qextAll[6,:] = tmp[:,4]
        cosinglescatAll[6,:] = tmp[:,5]
    

    else:
        for k in range(0,2):
            if (k == 0):
                radius = radiussn
            else:
                radius = radiusice

            if (radius == 0.5e-4):
                file = 'IN_Mie(r=0.05mm).txt'
            elif (radius == 1.0e-4):
                file = 'IN_Mie(r=0.1mm).txt'
            elif (radius == 2.0e-4):
                file = 'IN_Mie(r=0.2mm).txt'
            elif (radius == 3.5e-4):
                file = 'IN_Mie(r=0.35mm).txt'
            elif (radius == 5.0e-4):
                file = 'IN_Mie(r=0.5mm).txt'
            elif (radius == 1.0e-3):
                file = 'IN_Mie(r=1.0mm).txt'
            elif (radius == 2.5e-3):
                file = 'IN_Mie(r=2.5mm).txt'
            else:
                print('Grain size is ', radius)
                raise Exception("No Mie scattering file available for this grain size radius")

            # !!READ(ui,*)
            # ! DO ii = 1,bandmax
            # !  READ(ui,*) lambda(ii),dlambda(ii),asym(ii),qext(ii),cosinglescat(ii)
            # ! ENDDO

            # ! CLOSE(ui)

            if (k == 0):
                tmp = np.loadtxt(file,skiprows=0)
                Lambda = tmp[:,0]
                dlambda = tmp[:,1]
                asymsn = tmp[:,2]
                qextsn = tmp[:,3]
                cosinglescatsn = tmp[:,4]
            else:
                tmp = np.loadtxt(file,skiprows=0)
                Lambda = tmp[:,0]
                dlambda = tmp[:,1]
                asymice = tmp[:,2]
                qextice = tmp[:,3]
                cosinglescatice = tmp[:,4]


    return  (Lambda, SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear, 
             lambdaAll, dlambdaAll, asymAll, qextAll, cosinglescatAll,
             dlambda,asymsn, qextsn, cosinglescatsn,
             asymice, qextice, cosinglescatice)
