import numpy as np
import pandas as pd 
from datetime import datetime 
from numba import njit
import logging

from globals import *
from info import * 
from input import interp_data, checkdata, inputradpen, get_errorflag
from skintemperature import surft, tskin
from fluxes import turbhf, energybalance, get_z0m
from snowmodel import entemp, snowheight
from snowgrid import initgrid, initsnow
from massbalance import set_precip, set_height
from newalbedo import initgrains, newalbedo
from inittables import inittables
from radiation import cloudcover, calcalbedo, RadiationPenetration

def sebmodel_core(FORCING, indY, indX):
    """ SEB model core function, which perform the calculations on one core. 
    Params
    ======
    FORCING: xarray.Dataset
      xarray dataset which contain one grid point
    indY: 
    indX:
    Returns
    ======
    Returns all calculated variables of one grid point
    """
    # Get data from file
    #--------------------------------------------
    time = FORCING.time.values
    dhour = (((FORCING['time']-FORCING['time'].shift(time=1)).bfill(dim='time')) / np.timedelta64(1, 's')).values # seconds between forcing samples
    hour = pd.to_datetime(time).hour + pd.to_datetime(time).minute/60 # numeric hour
    day = np.asarray((FORCING['time'].values-FORCING['time'].values[0]).astype('timedelta64[D]')/np.timedelta64(1, 'D'), dtype = 'int') # days since start
    ndays = (FORCING['time'].values[-1]-FORCING['time'].values[0]).astype('timedelta64[D]')//np.timedelta64(1, 'D') # amount of days
    nmonths = (FORCING['time'].values[-1]-FORCING['time'].values[0]).astype('timedelta64[M]')//np.timedelta64(1, 'M') # amount of months
    # cmonth = pd.Timestamp(time).month # numeric cumulative motnhs
    nt = len(time)  
    if lhourlysnowout == 1:
        nt_layers = nt # Write snow layers every timestep
    else:   
        nt_layers = ndays  # Write snow layers every day
    
    # Local variables

    # surface
    Ndata = 0
    _Sin= np.full(nt, np.nan)
    _Sout = np.full(nt, np.nan)
    _Lin = np.full(nt, np.nan)
    _Loutmod = np.full(nt, np.nan)
    _Loutobs = np.full(nt, np.nan)
    _SH = np.full(nt, np.nan)
    _LE = np.full(nt, np.nan)
    _GH = np.full(nt, np.nan)
    _Restsource = np.full(nt, np.nan)
    _Source = np.full(nt, np.nan)
    _sumdivs = np.full(nt, np.nan)
    _T = np.full(nt, np.nan)
    _P = np.full(nt, np.nan)
    _WS = np.full(nt, np.nan)
    _q = np.full(nt, np.nan)
    _T0 = np.full(nt, np.nan)
    _q0 = np.full(nt, np.nan)
    _T2m = np.full(nt, np.nan)
    _q2m = np.full(nt, np.nan)
    _WS10m = np.full(nt, np.nan)
    _z0m = np.full(nt, np.nan)
    _dens_lay1 = np.full(nt, np.nan)
    _temp_lay1 = np.full(nt, np.nan)
    _dz_lay1 = np.full(nt, np.nan)

    # mass balance
    _dsnowacc = np.full(nt, np.nan)
    _icemelt = np.full(nt, np.nan)
    _hsnowmod  = np.full(nt, np.nan)
    _runoff = np.full(nt, np.nan)
    _surfwater = np.full(nt, np.nan)
    _melt = np.full(nt, np.nan)
    _surfmelt = np.full(nt, np.nan)
    _sumdrift = np.full(nt, np.nan)
    _subl = np.full(nt, np.nan)
    _icemeltmdt = np.full(nt, np.nan)
    _precip = np.full(nt, np.nan)
    _precipdt = np.full(nt, np.nan)
    _surfmeltdt = np.full(nt, np.nan)
    _meltdt = np.full(nt, np.nan)
    _runoffdt = np.full(nt, np.nan)
    _subldt = np.full(nt, np.nan)
    _sumwater = np.full(nt, np.nan)
    _topwater = np.full(nt, np.nan)
    _air_content= np.full(nt, np.nan)
    _effective_air_content = np.full(nt, np.nan)
    _errorflag = np.full(nt+1, 0.0)

    # subsurface layers
    if lwritelayers == 1:
        _z = np.full((nt_layers,nlmax), np.nan)
        _dz = np.full((nt_layers,nlmax), np.nan)
        _temp = np.full((nt_layers,nlmax), np.nan)
        _dens = np.full((nt_layers,nlmax), np.nan)
        _kice = np.full((nt_layers,nlmax), np.nan)
        _cpice = np.full((nt_layers,nlmax), np.nan)
        _rhocp = np.full((nt_layers,nlmax), np.nan)
        _energy = np.full((nt_layers,nlmax), np.nan)
        _lid = np.full((nt_layers,nlmax), np.nan)
        _mass = np.full((nt_layers,nlmax), np.nan)
        _grainsize = np.full((nt_layers,nlmax), np.nan)
        _water = np.full((nt_layers,nlmax), np.nan)
        _ice = np.full((nt_layers,nlmax), np.nan)
        _dsdz = np.full((nt_layers,nlmax), np.nan)
        _refrfrac = np.full((nt_layers,nlmax), np.nan)
    
    #--------------------------------------------
    # Prepare input array for SEB model 
    #--------------------------------------------
    ibuf = np.zeros([nt,colmax])
    ibuf[:,0] = FORCING.time.dt.year[:] # Year
    ibuf[:,1] = FORCING.time.dt.day[:] # Day
    ibuf[:,2] = FORCING.time.dt.hour[:] # Hour
    ibuf[:,3] = FORCING.T.values + 273.15 # air temperature in K
    ibuf[:,4] = FORCING.P.values * 100 # atmospheric pressure in Pa
    ibuf[:,5] = FORCING.WS.values # Wind speed in m s-1
    ibuf[:,6] = FORCING.Sin.values # Incoming shortwave radiation, definite positive in W m-2
    ibuf[:,7] = FORCING.Sout.values # Outgoing shortwave radiation, definite positive in W m-2
    ibuf[:,8] = FORCING.alb.values # Surfce broadband albedo between 0 and 1
    ibuf[:,9] = FORCING.q.values / 1000 # air specific humidity in kg kg-1
    ibuf[:,10] = FORCING.Lin.values # Incoming longwave radiation, definite positive in W m-2
    ibuf[:,11] = FORCING.Lout.values # Outgoing longwave radiation, definite positive in W m-2
    ibuf[:,16] = FORCING.precip.values # Precipitation in mm w.e
    if version == '1D':
        ibuf[:,12] = FORCING.z0m.values # Surface aerodynamic roughness length for momentum in m
        ibuf[:,13] = FORCING.zt.values # Height of T & q sensor above surface in m
        ibuf[:,14] = FORCING.zm.values # Height of WS sensor above surface in m
        ibuf[:,15] = FORCING.Serie.values # Cumulative surface level since start in m.
        zenith = FORCING.zenith.values # Solar zenith angle in degrees
        awsid = FORCING.AWSid.values # Solar zenith angle in degrees
        paramerror = FORCING.Error.values # Solar zenith angle in degrees
        Hmax = Hmax_info
    elif version == '2D':
        # ibuf[:,12] = np.full(nt, FORCING.z0m.values) #TODO implement map of albedo and z0m
        ibuf[:,12] = np.full(nt,z0msn)
        ibuf[:,13] = np.full(nt, 2.0) # zt
        ibuf[:,14] = np.full(nt, 10.) #zm
        ibuf[:,15] = np.full(nt, 0.) #Serie
        zenith = np.full(nt, 0.)
        awsid = np.full(nt, 1)
        paramerror = np.full(nt, 1000000000)
        Hmax = np.float64(1.0) #  np.float64(FORCING.Hice.values)
        albice = np.full(nt, 0.5) # np.float64(FORCING.BIA.values)
        # Hmax = np.float64(FORCING.Hice.values)
        # albice = np.float64(FORCING.BIA.values)
        # T10m = np.float64(FORCING.T10m.values)
        # dsnow = np.float64(FORCING.dsnow.values)
        # dfirn = np.float64(FORCING.dsnow.values)
        
    

    # ibuf[:,17] = FORCING.accum.values # accum determined from serie
    # ibuf[:,18] = FORCING.melt.values # melt determined from serie
    # ibuf[:,19] = FORCING.drift.values # drift determined from serie

    sumdivs = 0.
    water = np.zeros(nlmax)
    ice = np.zeros(nlmax)
    freshfrac = np.zeros(nlmax)
    dtdz = np.zeros(nlmax)
    dsdz = np.zeros(nlmax)
    refrfrac = np.zeros(nlmax)

    nlsnow = 1
    lidmax = 1
    vink = 0
    precipsum = 0.
    freshsnow = 0.
    tempprecip = 0.
    icemelt = 0.
    icemeltmin = 10.
    surfmelt = 0.
    sumdrift = 0.
    runoff = 0.
    melt = 0.
    subl = 0.
    surfwater = 0.
    Hice = Hmax #  
    icemeltmdt = 0.
    meltdt = 0.
    surfmeltdt = 0.
    runoffdt = 0.
    subldt = 0.
    sumwater = 0.
    topwater = 0.
    air_content = 0
    effective_air_content = 0
    racc_old = ibuf[0,15]
    acclevel = 0.
    climracc_old = ibuf[0,15]
    climserie =  ibuf[0,15]
    climacclevel = 0.
    climmeltlevel = 0.
    melt_l = 0.
    corrsnow = 0. # MvT TODO: implement SETTOERROR routines
    snowdays = 0.
    
    rndnr_z0 = np.random.uniform(low=0.0, high=1.0, size=None)
    rndnr_alb = np.random.uniform(low=0.0, high=1.0, size=None)
    
    if ((lalbedo == 4) | (lalbedo == 5)):
        print('Starting inittables')
        (TVals, DTDZVals, DENSVals, TAUMAT, KAPMAT, DR0MAT) = inittables()
        print('End inittables')
   
    if (penetration != 0 ):
        if lcomment == 1: print('Starting inputradpen')
        (Lambda, SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear, 
        lambdaAll, dlambdaAll, asymAll, qextAll, cosinglescatAll,
        dlambda,asymsn, qextsn, cosinglescatsn,
        asymice, qextice, cosinglescatice) = inputradpen()

    if lcomment == 1: print('Starting initgrid')
    z, dz, depthl, lid, dsnowacc, hmass, nl, nlinit  = initgrid(dsnow_glob)
    if lcomment == 1: print('Starting initsnow')
    dens, temp, mass, cpice, rhocp, kice, energy, zrad = initsnow(z,dz,lid,nl)
    if lcomment == 1: print('Starting initgrains')
    grainsize, radfresh = initgrains(lid,nl)
    
    logging.info("Starting get_errorflag")

    
    # Find errorflag
    if lcomment == 1: print('Starting get_errorflag')
    _errorflag = get_errorflag(awsid, paramerror,nt)
     
    # Check data
    if lcomment == 1: print('Starting checkdata')
    (buf, dsnowr, dsnow, hsnowstart, hsnowmod, 
    dsnowh, corrsnow, alb_old, racc_old, acclevel, 
    climracc_old, climserie, climacclevel, climmeltlevel) = \
    checkdata(
        ibuf,nt, racc_old, acclevel, climracc_old,
        climserie, climacclevel, climmeltlevel
        )

    # initialise output per timestep
    __Sin  = 0
    __Sout = 0
    __Lin  = 0
    __Loutmod = 0
    __Loutobs = 0
    __SH = 0
    __LE = 0
    __GH = 0
    __Restsource = 0
    __Source = 0
    __sumdivs = 0
    __T = 0
    __P = 0
    __WS = 0
    __q = 0
    __T0 = 0
    __q0 = 0
    __T2m = 0
    __q2m = 0
    __WS10m = 0
    __z0m = 0
    __dens_lay1 = 0
    __temp_lay1 = 0
    __dz_lay1 = 0

    # Mass balance diagnostics
    __icemelt = 0
    __dsnowacc = 0
    __hsnowmod = 0
    __runoff =0
    __surfwater = 0
    __melt = 0
    __surfmelt = 0
    __sumdrift = 0
    __subl = 0
    __precip = 0
    __precipdt = 0
    __icemeltmdt = 0
    __surfmeltdt = 0
    __meltdt = 0
    __runoffdt = 0
    __subldt = 0
    __sumwater = 0
    __topwater = 0
    __air_content = 0
    __effective_air_content = 0


                    
    #--------------------------------------------
    # TIME LOOP 
    # First, loop over nt time intervals from forcing  (tinterval)
    # Then loop over nstep amount of timesteps (tstep) within each forcing interval 
    #--------------------------------------------    
        
    for ii in np.arange(nt):      
        
        # print('ii = ',ii) 

        valerrorgap = 100
        if ((lerrorgap == 0) | ((_errorflag[ii] > 29) & (_errorflag[ii+1] > 29))): valerrorgap = 8
        if (lerrorgap >= 1): valerrorgap = 100

        if ((_errorflag[ii]< valerrorgap) & (_errorflag[ii+1] < valerrorgap)):	
        # only when enough valid observations are present in case valerrorgap = 8, 
        # very large values of valerrorgap and it will use the in input set values usually based on interpolation

            if ii > 0:
                tinterval = (time[ii] - time[ii-1])/np.timedelta64(1, "s")
                nstep = tinterval/tstep
            else:
                tinterval = (time[ii+1] - time[ii])/np.timedelta64(1, "s")
                nstep = tinterval/tstep
            for jj in np.arange(nstep):
                
                # Determine mass balance 

                # Interpolate to obtain values on higher time resolution than input data 
                sbuf, szenith, hsnow, drift = interp_data(buf,zenith,jj,ii,nstep,lid[0],nt)
                
                # Formulation calculation surface temperature
                if (tcalc == 1):	# from Lwout observations as start or definite value	
                    t0 = (sbuf[11]/StefBoltz)**(0.25)	
                elif (tcalc == 2):		# chosen equal to temperature of upper most snow/ice layer
                    t0 = temp[0]
                elif ((tcalc == 3) | (tcalc == 4)):	# extrapolated from upper 1-3 snow/ice layers, in case of 4: start value for iterations
                    t0 = surft(temp[0],temp[1],dz[0],dz[1],z[0],z[1])	# extrapolation from upper most layers snow/ice
                # Convert precipitation to m snow, calculate/set fresh snow density
                precip, water, rhosn = set_precip(t0,sbuf,water)
                    
                # In case snowheight is restricted to accumulation observations, correct for data gaps
                if (luseacc >= 2): 
                    (temp, water, ice, mass, dens, dz, z,
                    lid, nl, nlsnow, corrsnow, precip,
                    drift, icemelt, dsnowh,melt_l) = \
                    set_height(
                        sbuf,t0,dsnowh,dsnowr,rhosn,corrsnow,
                        temp, water, ice, mass, dens, dz, z, 
                        lid, nl, nlsnow, melt_l, icemelt, precip, drift
                        )
                # Add snowfall layer or remove layers
                if ((precip != 0.) | (drift != 0.)):
                    (nl, z, dz, temp, dens, kice, dtdz, cpice,
                    water, mass, ice, lid, grainsize, refrfrac,
                    sumsnow, summass, dsnowacc, hmass, vink, 
                    precipsum,freshsnow, tempprecip, dsnowr, dsnowh, sumwater,
                        topwater, air_content, effective_air_content) = \
                    snowheight(
                        t0, temp, water, dens, mass, ice, z, dz, lid, 
                        grainsize, energy, kice, dtdz, refrfrac, cpice, precip, 
                        drift,freshfrac, lidmax, rhosn, radfresh, nl, nlsnow,
                        dsnowr, vink, precipsum, freshsnow, tempprecip, dsnowh)
    
                # Compute cloud cover
                cloud = cloudcover(sbuf[3],sbuf[10]) # only necessary in aid of radiation penetration
                                                     # and for calculating albedo in new parameterization scheme (lalbedo==4)
                                                     # but calculated always for extra information 3 = temp, 10 = Lin
                
                # Compute albedo and net shortwave radiation
                if (lalbedo == 0):
                    if (lsnet == 0):
                        Snet = sbuf[6]-sbuf[7]              # Snet based on observed Sin and Sout
                    elif (lsnet == 1):
                        Snet = sbuf[7]*((1/sbuf[8])-1)      # Snet based on Sout and albedo, albedo is preferably the 24 hours running mean that removes daily cycle due to albedo!!! But also removes some of the error induced by tilt of the sensor.
                        sbuf[6] = sbuf[7]/sbuf[8]           # necessary to keep output consistent with what is used for the calculations
                    elif (lsnet == 2):
                        Snet = sbuf[6]*(1-sbuf[8])      # Snet based on Sin and albedo, albedo is preferably the 24 hours running mean that removes daily cycle due to albedo!!! But also removes some of the error induced by tilt of the sensor.
                        sbuf[7] = sbuf[6]*sbuf[8]           # necessary to keep output consistent with what is used for the calculations

                else:
                    if ((lalbedo == 4) | (lalbedo == 5)):
                        # MvT TODO check NEWALBEDO 
                        albedo = newalbedo(sbuf,szenith,cloud,z,dz,lid,nlsnow,
                                           dtdz,dens,temp,water, mass, refrfrac, 
                                           freshfrac, grainsize,radfresh,
                                           DENSVals,TVals,DTDZVals,TAUMAT,KAPMAT,DR0MAT)
                    else:
                        albedo, snowdays, alb_old = calcalbedo(szenith,lid,hmass,dsnowr,rhosn,water,precip,t0,snowdays,alb_old)
                    if (lalbedo == 5): 
                        if (lsnet == 0):
                            Snet = sbuf[6]-sbuf[7]              # Snet based on observed Sin and Sout
                        elif (lsnet == 1):
                            Snet = sbuf[7]*((1/sbuf[8])-1)      # Snet based on Sout and albedo, albedo is preferably the 24 hours running mean that removes daily cycle due to albedo!!! But also removes some of the error induced by tilt of the sensor.
                            sbuf[6] = sbuf[7]/sbuf[8]           # necessary to keep output consistent with what is used for the calculations
                        elif (lsnet == 2):
                            Snet = sbuf[6]*(1-sbuf[8])      # Snet based on Sin and albedo, albedo is preferably the 24 hours running mean that removes daily cycle due to albedo!!! But also removes some of the error induced by tilt of the sensor.
                            sbuf[7] = sbuf[6]*sbuf[8]           # necessary to keep output consistent with what is used for the calculations
                    else:
                        Snet = sbuf[7]*((1/albedo)-1)	        # Snet based on Sout and parameterised albedo
                        sbuf[6] = Snet + sbuf[7]			    # necessary to keep output consistent with what is used for the calculation
                
                # Compute penetration of shortwave radiation 
                if ((penetration == 1) & (sbuf[6] > 0.0)):
                    dsdz, sumdivs = RadiationPenetration(
                        z,dz,grainsize,dens,sbuf,Snet,lid,cloud,
                        Lambda, SolarPlateauClear, SolarPlateauCloudy, SolarSeaClear, 
                        lambdaAll, dlambdaAll, asymAll, qextAll, cosinglescatAll,
                        dlambda,asymsn, qextsn, cosinglescatsn,
                        asymice, qextice, cosinglescatice
                        )
                else:
                    dsdz = np.zeros(nlmax)
                    sumdivs = 0.

                # Compute roughness length for momentum
                (z0m, Hice, dHice) \
                = get_z0m(
                    sbuf,lid,rndnr_z0, 
                    dsnowacc, icemeltmdt, Hice, Hmax
                    )
                sbuf[12] = z0m
    
                # Formulation energy balance from skin layer temperature
                if (tcalc == 4):	# Skin layer formulation for surface temperature
                    t0, ustar, SH, LE, Ch, Cq, Chn, Cqn, t2m, q0, q2m, ws10m, densair = tskin(t0,Snet,sumdivs,kice,temp,dz,sbuf)
                    source, GH = energybalance(t0,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf)
                else: # formulation energy balance from t0 previous time step, excess heat put into first layer
                    ustar, SH, LE, Ch, Cq, Chn, Cqn, t2m, q0, q2m, ws10m, densair = turbhf(-1,t0,sbuf) # Calculate turbulent fluxes
                    source, GH = energybalance(t0,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf)
                    if (t0 >= Tkel):
                        source, GH = energybalance(Tkel,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf)		
                    if (source < 0): 
                        source=0.
                                
                # check difference from energy balance (should be =0 in case of no melt else > 0)
                restsource, GH = energybalance(t0,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf) #- source
                    
                # Calculate net longwave radiation from model results
                Lout = emis*StefBoltz*(t0)**4
                Lnet = (sbuf[10] - Lout)
                
                # # Calculate new englacial profiles
                (temp, kice, cpice, rhocp, energy, z, dz,
                water, ice, mass, dens, grainsize, dtdz,
                lid, vink, icemelt, icemeltmin,
                dsnowr, surfmelt, surfmeltdt, sumdrift, nl, nlsnow,
                runoff, runoffdt, melt, meltdt, subl, subldt, surfwater, refrfrac, 
                icemeltmdt, dsnowacc,hsnowmod_l, dsnowh, exposed_water,source) = \
                entemp(
                    t0, LE, nl, nlinit, nlsnow, lidmax, kice, cpice, rhocp, 
                    water, ice, mass, dens, lid, dz, z, temp, grainsize, dsdz,
                    source, icemelt, icemeltmin, dsnowr, precip, drift, surfmelt,
                    sumdrift, runoff, melt, subl, rhosn, surfwater, refrfrac, vink,dsnowh)
                
                # Averaged output 
                Ndata += 1
                __Sin  +=  sbuf[6]
                __Sout += sbuf[7]
                __Lin  += sbuf[10]
                __Loutmod += Lout
                __Loutobs += sbuf[11]
                __SH += SH
                __LE += LE
                __GH += GH 
                __Restsource += restsource
                __Source +=  source
                __sumdivs += sumdivs
                __T +=   sbuf[3]
                __P +=   sbuf[4]
                __WS +=   sbuf[5]
                __q +=   sbuf[9]
                __T0 +=  t0
                __q0 += q0
                __T2m += t2m
                __q2m += q2m
                __WS10m += ws10m
                __z0m += sbuf[12]
                __icemelt += icemelt
                __dsnowacc += dsnowacc
                __hsnowmod+= hsnowmod_l
                __runoff += runoff
                __surfwater += surfwater
                __melt +=  melt
                __surfmelt += surfmelt
                __sumdrift += sumdrift
                __subl +=  subl
                __precip += precipsum
                __dens_lay1 += dens[0]
                __temp_lay1 += temp[0]
                __dz_lay1 += dz[0]
                __sumwater += sumwater
                __topwater += topwater
                __air_content += air_content
                __effective_air_content += effective_air_content

                # Summed output 
                __icemeltmdt +=  icemeltmdt
                __runoffdt +=  + runoffdt
                __precipdt +=  + precip
                __surfmeltdt +=  + surfmeltdt
                __meltdt +=  + meltdt
                __subldt +=  + subldt

                # Averaged results per forcing time step 
                if (jj == (dhour[ii]-tstep)/tstep):
                    
                    # Instantaneous output 
                    _Sin[ii] = __Sin / Ndata
                    _Sout[ii] = __Sout / Ndata
                    _Lin[ii] = __Lin / Ndata
                    _Loutmod[ii] = __Loutmod / Ndata
                    _Loutobs[ii] = __Loutobs / Ndata
                    _SH[ii] = __SH / Ndata
                    _LE[ii] = __LE / Ndata
                    _GH[ii] = __GH / Ndata
                    _Restsource[ii] = __Restsource
                    _Source[ii] = __Source / Ndata
                    _sumdivs[ii] = __sumdivs / Ndata
                    _T[ii] = __T / Ndata
                    _P[ii] = __P / Ndata
                    _WS[ii] = __WS / Ndata
                    _q[ii] = __q / Ndata
                    _T0[ii] = __T0 / Ndata
                    _q0[ii] = __q0 / Ndata
                    _T2m[ii] = __T2m / Ndata
                    _q2m[ii] = __q2m / Ndata
                    _WS10m[ii] = __WS10m / Ndata
                    _z0m[ii] = __z0m / Ndata
                    _icemelt[ii] = __icemelt / Ndata
                    _dsnowacc[ii] = __dsnowacc / Ndata
                    _hsnowmod[ii] = __hsnowmod / Ndata
                    _runoff[ii] = __runoff / Ndata
                    _surfwater[ii] = __surfwater / Ndata
                    _melt[ii] = __melt / Ndata
                    _surfmelt[ii] = __surfmelt / Ndata
                    _sumdrift[ii] = __sumdrift / Ndata
                    _subl[ii] = __subl / Ndata
                    _precip[ii] = __precip / Ndata
                    _dens_lay1[ii] = __dens_lay1 / Ndata
                    _temp_lay1[ii] = __temp_lay1 / Ndata
                    _dz_lay1[ii] = __dz_lay1 / Ndata
                    _sumwater[ii] = __sumwater / Ndata
                    _topwater[ii] = __topwater / Ndata
                    _air_content[ii] = __air_content / Ndata
                    _effective_air_content[ii] = __effective_air_content / Ndata

                    # Summed output 
                    _icemeltmdt[ii] = __icemeltmdt
                    _precipdt[ii] = __precipdt
                    _surfmeltdt[ii] = __surfmeltdt
                    _meltdt[ii] = __meltdt
                    _runoffdt[ii] = __runoffdt  
                    _subldt[ii] = __subldt

                    # Subsurface layers
                    if (lhourlysnowout == 1) & (lwritelayers == 1):
                        _z[ii,:] = z 
                        _dz[ii,:] = dz 
                        _temp[ii,:] = temp 
                        _dens[ii,:] = dens 
                        _kice[ii,:] = kice 
                        _cpice[ii,:] = cpice 
                        _rhocp[ii,:] = rhocp 
                        _energy[ii,:] = energy 
                        _lid[ii,:] = lid 
                        _mass[ii,:] = mass 
                        _grainsize[ii,:] = grainsize 
                        _water[ii,:] = water 
                        _ice[ii,:] = ice 
                        _dsdz[ii,:] = dsdz  
                        _refrfrac[ii,:] = refrfrac

                    # reset output per timestep 
                    Ndata  = 0
                    __Sin  = 0
                    __Sout = 0
                    __Lin  = 0
                    __Loutmod = 0
                    __Loutobs = 0
                    __SH = 0
                    __LE = 0
                    __GH = 0
                    __Restsource = 0
                    __Source = 0
                    __sumdivs = 0
                    __T = 0
                    __P = 0
                    __WS = 0 
                    __q = 0
                    __T0 = 0
                    __q0 = 0
                    __T2m = 0
                    __q2m = 0
                    __WS10m = 0
                    __z0m =0
                    __icemelt = 0
                    __dsnowacc = 0
                    __hsnowmod = 0
                    __runoff =0
                    __surfwater = 0
                    __melt = 0
                    __surfmelt = 0
                    __sumdrift = 0
                    __subl = 0
                    __precip = 0
                    __dens_lay1 = 0
                    __temp_lay1 = 0
                    __dz_lay1 = 0
                    __icemeltmdt = 0
                    __precipdt = 0
                    __surfmeltdt = 0
                    __meltdt = 0
                    __runoffdt = 0
                    __subldt = 0
                    __sumwater = 0
                    __topwater = 0
                    __air_content= 0
                    __effective_air_content = 0          
            # end (jj) loop                   
                        
            # Results per day
            if (hour[ii] % 24 == 0):
                dd = day[ii]               
                # Subsurface layers
                if (lhourlysnowout == 0) & (lwritelayers == 1):
                    _z[dd,:] = z
                    _dz[dd,:] = dz
                    _temp[dd,:] = temp
                    _dens[dd,:] = dens
                    _kice[dd,:] = kice
                    _cpice[dd,:] = cpice
                    _rhocp[dd,:] = rhocp
                    _energy[dd,:] = energy
                    _lid[dd,:] = lid
                    _mass[dd,:] = mass
                    _grainsize[dd,:] = grainsize
                    _water[dd,:] = water
                    _ice[dd,:] = ice
                    _dsdz[dd,:] = dsdz
                    _refrfrac[dd,:] = refrfrac     
                # print progres every 100 days for 1D version only
                if (day[ii] % 100 == 0) & (version == '1D'):
                    print(round(100*day[ii]/ndays),'% done')
    # end (ii) loop
    if lwritelayers == 0:
        _z = []
        _dz = []
        _temp = []
        _dens = []
        _kice = []
        _cpice = []
        _rhocp = []
        _energy = []
        _lid = []
        _mass = []
        _grainsize = []
        _water = []
        _ice = []
        _dsdz = []
        _refrfrac = []

    return (indY,indX,_Sin,_Sout,_Lin,_Loutobs,_Loutmod,_SH,_LE,_GH,
        _Restsource, _Source, _sumdivs, _T, _P, _WS, _q, _T0, _q0, _T2m, _q2m, _WS10m, _z0m, _z, _dz,
        _temp, _dens, _kice, _cpice, _rhocp, _energy, _lid, _mass,
        _grainsize, _water, _ice, _dsdz, _refrfrac,
        _icemelt, _icemeltmdt, _dsnowacc, _hsnowmod, _runoff, _runoffdt, _surfwater,
        _melt, _meltdt, _surfmelt, _surfmeltdt, _sumdrift, _subl, _subldt, _precip, _precipdt,
        _dens_lay1, _temp_lay1, _dz_lay1, _sumwater, _topwater, _air_content, _effective_air_content,_errorflag[:-1])
