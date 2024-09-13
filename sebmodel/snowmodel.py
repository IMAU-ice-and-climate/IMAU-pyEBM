import numpy as np
from numba import njit

from globals import *
from info import *
from snowgrid import resizegrid, redefgrid
from routines import conduc, irreducible
from functions import resetdens

                
              
#===============================================================================
@njit
def entemp(t0, LE, nl, nlinit, nlsnow, lidmax, kice, cpice, rhocp, water, 
           ice, mass, dens, lid, dz, z, temp, grainsize, dsdz, source,
           icemelt, icemeltmin, dsnowr, precip, drift, surfmelt, sumdrift, runoff, melt, subl,
           rhosn, surfwater, refrfrac, vink,dsnowh):
    '''
    Routine that calculates the subsurface temperature

    1) Adapts the subsurface grid based on densification (RESIZEGRID) or if layers are too thin (REDEFGRID) 
    2) Computes the subsuface heat conductivity and heat capacity 
    3) Computes the irreducible water content 
    4) Calculates the new subsurface temperature including the radiation penetration based on the following routine:
        # ! rho*cp*dT/dt = d/dz(K*dT/dz)
        # ! dT/dt = (K/rho*cp)*dT/dz2
        # ! T(i+1) = T(i)+(d/dz(K*dT/dz))*dt/(rho*cp)
        # ! properties of the layer are defined at half the layer thickness!!
    5) Compute amount of snow melt, refreezing and percolation (SNOWCONTENT)
    6) Calculate condensation / evaporation
    7) Compute change in ice or snow thickness

    Parameters
    ----------
    t0 : float
        surface temperature
    LE : float
        Latent heat flux, positive downwards in Wm-2
    nl : int
        Number of subsurface layers
    nlinit : int
        Initial number of subsurface layers
    nlsnow : int
        Number of snow layers
    lidmax = int
        maximum subsurface layer ID (for multiyear firn)
    kice : ndarray (1D float)
        subsurface effective heat conductivity
    cpice : ndarray (1D float)
        Volumetric heat capacity of ice
    rhocp : ndarray (1D float)
        Density times volumetric heat capacity
    water : ndarray (1D float)
        subsurface liquid water content profile
    ice : ndarray (1D float)
        subsurface ice content profile
    mass : ndarray (1D float)
        subsurface mass of layer profile
    dens : ndarray (1D float)
        subsurface density profile
    lid : ndarray (1D float)
        layer ID (0 = ice, 1 = snow, 2 and larger = firn)
    z : ndarray (1D float)
        depth of layer centers, positive downwards
    dz : ndarray (1D float)
        layer thickness
    temp : ndarray (1D float)
        subsurface temperature profile
    grainsize : ndarray (1D float)
        subsurface grainsize
    dsdz : ndarray (1D float)
        Radiation penetration
    source : float
        Energy available for surface melt in Wm-2    
    icemelt : float
        Thickness ice melt. If no snowpack, all melt is ice melt and lowers surface, accumulated in this parameter
    icemeltmin : float
        Temporary variable used for icemeltout
    dsnowr : float
        Snow height that corrects for small snow layers that are added to the ice when the layer is too thin to sustain
    precip : float
        total precipitation during time step in mm w.e.
    drift : float
        snow drift in mm w.e during time step
    surfmelt : float
        Cumulative surface melt in mm w.e.
    sumdrift : float
        cumulative snow drift in mm w.e
    runoff : float
        Surface runoff in mm w.e.
    melt : float
        surface melt in mm w.e
    subl : float
        Surface sublimation in mm w.e.
    rhosn : float
        Density of snowfall set in routine SET_PRECIP
    surfwater : float
        amount of water on the surface 
    refrfrac : ndarray (1D float)
        subsurface profile of refrozen snow  
    vink : int
        switch to turn on redefining grid when layers become smaller than 0.5*dz0
    dsnowh : float
        Thickness snow layer, thickness ice melt in m ice 
 
    Returns
    -------  
    temp : ndarray (1D float)
        subsurface temperature profile
    kice : ndarray (1D float)
        subsurface heat conductivity profile
    cpice : ndarray (1D float)
        surface heat capacity profile
    rhocp : ndarray (1D float)
        Density times volumetric heat capacity
    energy : ndarray (1D float)
        energy potential layer
    z : ndarray (1D float)
        depth of layer centers, positive downwards
    dz : ndarray (1D float)
        layer thickness
    water : ndarray (1D float)
        subsurface liquid water content profile
    ice : ndarray (1D float)
        subsurface ice content profile
    mass : ndarray (1D float)
        subsurface mass of layer profile
    dens : ndarray (1D float)
        subsurface density profile
    grainsize : ndarray (1D float)
        subsurface grainsize profile
    dtdz : ndarray (1D float)
        Vertical temperature gradient
    lid : ndarray (1D float)
        layer ID (0 = ice, 1 = snow, 2 and larger = firn)
    vink : int
        switch to turn on redefining grid when layers become smaller than 0.5*dz0
    icemelt : float
        Thickness ice melt. If no snowpack, all melt is ice melt and lowers surface, accumulated in this parameter
    icemeltmin : float
        Temporary variable used for icemeltout
    dsnowr : float
        Snow height that corrects for small snow layers that are added to the ice when the layer is too thin to sustain
    surfmelt : float
        Cumulative surface melt in mm w.e.
    sumdrift : float
        cumulative snow drift in mm w.e
    nl : int
        Number of subsurface layers
    nlsnow : int
        Number of snow layers
    runoff : float
        Surface runoff in mm w.e.
    melt : float
        surface melt in mm w.e
    subl : float
        Surface sublimation in mm w.e.
    surfwater : float
        amount of water on the surface 
    refrfrac : ndarray (1D float)
        subsurface profile of refrozen snow   
    icemeltmdt : float
        ice ablation (melt + sublimation) during time step in m
    dsnowacc : float 
        snow layer in m snow
    hsnowmod_l : float
        temporary variable keeping track of snow height
    dsnowh : float
        Thickness snow layer, thickness ice melt in m ice   
        
    '''
    exposed_water = 0
    Kdtdz = np.zeros(nlmax)
    dtdz = np.zeros(nlmax)
    energy = np.zeros(nlmax)

    # if (tpdens > 0):
    z, dz, dsnowh, dsnowacc, hmass, vink = \
        resizegrid(mass, dens, lid, dz, z, nl, nlinit, lidmax,vink) #! in case of 0, constant density profile, layers remain same size, except top because of melt
    if (vink >= 1):
        water, ice, mass, dens, lid, dz, z, temp, grainsize, energy, kice, cpice, dtdz, vink, nl,dsnowr, nlsnow = \
            redefgrid(water, ice, mass, dens, lid, dz, z, nl, nlsnow, temp, grainsize, t0,dsnowr,vink) # redefine grid in case layers become thinner than 0.5*dz0
    for il in range(0,nl):
        kice[il] = conduc(dens[il],temp[il])
        cpice[il] = 152.5 + 7.122 * temp[il]
        rhocp[il] = dens[il]*cpice[il]
           
    #kice(1) = 50.0
    irrwater = irreducible(dens,nl)

    dtdz[0] = (temp[0]-t0)/dz[0]
    Kdtdz[0] = kice[0]*dtdz[0]

    if ImpExp == 1:
        # Explicit Euler method
        for il in range(1,nl):
            dzl = 0.5*(dz[il]+dz[il-1])
            dtdz[il] = (temp[il]-temp[il-1])/dzl
            Kdtdz[il] = (1./(2.*dzl))*(kice[il]*dz[il]+kice[il-1]*dz[il-1])*dtdz[il]

        # dsdz = radiation penetration
        if ((tcalc == 2.) | (tcalc == 3)): dsdz[0] = -source + dsdz[0]*dz[0]	 # (interpolation upper layer, all energy entering first layer) 
        for il in range (0,nl-1):
            if (((tcalc == 2) | (tcalc == 3)) &  (il == 0)):
                temp[il] = temp[il] + (tstep/rhocp[il])*((Kdtdz[il+1] - dsdz[il]*dz[il])/dz[il])
            else:
                temp[il] = temp[il] + (tstep/rhocp[il])*((Kdtdz[il+1] - Kdtdz[il])/dz[il] - dsdz[il])
            kice[il] = conduc(dens[il],temp[il])
            cpice[il] = 152.5 + 7.122 * temp[il]
            rhocp[il] = dens[il]*cpice[il]
            energy[il] = dens[il]*dz[il]*cpice[il]*(Tkel-temp[il])

        temp[nl-1] = temp[nl-2]			# boundary condition, no heat flux to or from lower layers
        energy[nl-1] = energy[nl-2]

    elif ImpExp == 2:
        ### Implicit method
        # Code based on IMAU-FDM: https://github.com/IMAU-ice-and-climate/IMAU-FDM 
        # see book An Introduction To Computational Fluid Dynamics by Versteeg and
        # Malalasekera chapter 8.1 to 8.3 for implicit scheme and chapter 7.2 to 7.5
        # for description Thomas algorithm for solving the set of equations
        alpha = np.zeros(nlmax)
        beta = np.zeros(nlmax)
        D = np.zeros(nlmax)
        C = np.zeros(nlmax)
        CC = np.zeros(nlmax)
        A = np.zeros(nlmax)

        th  = 0.5 # 1 for full implicit, 0.5 for Crank Nicolson (semi implicit)
        # Bottom boundary condition, no flux
        f = dz[nl-1]/( dz[nl-1]+ dz[nl-2])
        An = ((1-f)*kice[nl-1]+f*kice[nl-2])*2./(dz[nl-1]+ dz[nl-2])
        AP0 = rhocp[nl-1]*dz[nl-1]/tstep
        alpha[nl-1] = An*th
        D[nl-1]= th*An + AP0
        C[nl-1] = An*(1.-th)*temp[nl-2] + (AP0-(1.-th)*An)*temp[nl-1]
        
        for il in range (nl-2,0,-1):
            cpice[il] = 152.5 + 7.122 * temp[il]
            f = dz[il]/( dz[il]+ dz[il-1])
            As = An
            An = ((1-f)*kice[il]+f*kice[il-1])*2./(dz[il]+ dz[il-1])
            AP0 = rhocp[il]*dz[il]/tstep
            beta[il] = As*th
            alpha[il] = An*th
            D[il] = th*(As+An) + AP0
            C[il] = As*(1-th)*temp[il+1] + An*(1-th)*temp[il-1] + (AP0-(1-th)*(As+An))*temp[il]

        # Top boundary condition, fixed temperature
        As = An
        AP0 = rhocp[0]*dz[0]/tstep
        Su = 2.*kice[0]*t0/dz[0]
        Sp = -2.*kice[0]/dz[0]
        beta[0] = As*th
        D[0] = th*(As-Sp) + AP0
        C[0] = As*(1.-th)*temp[1] + (AP0-(1.-th)*(As-Sp))*temp[0]+ Su

        # forward elimination
        A[nl-1] = alpha[nl-1]/D[nl-1]
        CC[nl-1] = C[nl-1]/D[nl-1]
        for il in range (nl-2,0,-1):
            A[il] = alpha[il]/(D[il]-beta[il]*A[il+1])
            CC[il] = (beta[il]*CC[il+1]+C[il])/(D[il]-beta[il]*A[il+1])
        CC[0] = (beta[0]*CC[1]+C[0])/(D[0]-beta[0]*A[1])

        # back-substitution
        temp[0] = CC[0] 
        kice[0] = conduc(dens[0],temp[0])
        cpice[0] = 152.5 + 7.122 * temp[0]
        rhocp[0] = dens[0]*cpice[0]
        energy[0]  = dens[0]*dz[0]*cpice[0]*(Tkel-temp[0])
        for il in range (1,nl):
            temp[il] = A[il]*temp[il-1] + CC[il]
            kice[il] = conduc(dens[il],temp[il])
            cpice[il] = 152.5 + 7.122 * temp[il]
            rhocp[il] = dens[il]*cpice[il]
            energy[il]  = dens[il]*dz[il]*cpice[il]*(Tkel-temp[il])

    if (t0 < Tkel): source = 0		# source is now only melt energy, not energy used for heating
    temp, water, ice, mass, dens, kice, cpice, grainsize, z, dz, \
        dtdz, energy, lid, dsnowr, surfmelt, surfmeltdt, sumdrift, nl, nlsnow, \
            runoff, runoffdt, melt, meltdt, subl, subldt, surfwater, refrfrac, exposed_water \
        = snowcontent(LE, t0, source, temp, water, ice, mass, dens, \
            kice, cpice, grainsize, z, dz, dtdz, energy, lid, irrwater, \
                dsnowr, surfmelt, drift, sumdrift, nl, nlsnow, runoff, melt, \
                    subl, rhosn, surfwater, refrfrac)    
    # Calculate condensation / evaporation
    ll = lv
    if (LE < 0.): ll = ls
    if ((LE > 0.) & (t0 < Tkel)): ll = ls
    cond = LE*tstep/ll				# in mm w.e.
    if ((cond > 0.)  & (t0 >= Tkel)): cond = 0.

    icemeltmdt = 0
    if (dsnowh <=0.):
        icemeltmdt = (source * tstep /lm)/dens[0] #+ cond/dens(1) !in m ice per time step
        icemelt = icemelt - (source * tstep /lm)/dens[0] #+ cond/dens[0]	# in m ice
        if (cond < 0):
            icemelt = icemelt + cond/dens[0]
            icemeltmdt = icemeltmdt + cond/dens[0]
        else:
            dsnowacc = dsnowacc + cond/dens[0]  

    # icemeltout: used for output MB only, (diagnostic parameter)
    icemeltout = icemelt + max(0.,dsnowr)
    if (icemeltout > icemeltmin):
        icemeltout = icemeltmin
    else:
        icemeltmin = icemeltout

    hsnowmod = dsnowh + icemelt + dsnowr # + cond/densice (dsnowr = extra correction for adding/removing layers when all snow has melted)
    if ((precip != 0.) | (drift != 0.)): dsnowacc = dsnowacc + dsnowr
    hsnowmod_l = hsnowmod

    # dsnowh = thickness snowpack
    # icemelt = if no snowpack, all melt is ice melt and lowers surface, accumulated in this parameter
    # dsnowr =  corrects for small snow layers that are added to the ice when the layer is too thin to sustain. 

    # Check stability of the solution with criterium: D*dt/dx2 < 0.25
    #stcrit=(kice(1)/rhocp(1))*tstep/(dz(1)**2)
    #dt = 0.25*tstep/stcrit
    #IF (stcrit.ge.0.25) THEN
    ##   IF (lcomment == 1) write(*,'(4(A,f9.4))') 'stcrit= ',stcrit,' dt = ',0.25*tstep/stcrit,'dt = ',tstep*1.,' dz = ',dz(1)
    #   write(uo1,'(4(A,f9.4))') 'stcrit= ',stcrit,' dt = ',0.25*tstep/stcrit,' dt = ',tstep*1.,' dz = ',dz(1)
    # ENDIF

    # Method as described by WJvdB
    stcrit2 = (2*kice[0]/rhocp[0])*tstep/(dz[0]**2)
    stcrit1 = (kice[0]/rhocp[0])*tstep/(dz[0]**2)
    stcrit = stcrit2+(1-stcrit2)*stcrit2-stcrit2*stcrit1
    dt = 0.25*tstep/stcrit1
    if (stcrit < 0) & (ImpExp == 1):
        print('stcrit= ',stcrit,' dt = ',0.25*tstep/stcrit1,' dt = ',tstep*1.,' dz = ',dz[0])

    return (temp, kice, cpice, rhocp, energy, z, dz, water, ice, mass, 
            dens, grainsize, dtdz, lid, vink, icemelt, icemeltmin, dsnowr, surfmelt, surfmeltdt,
            sumdrift, nl, nlsnow, runoff, runoffdt, melt, meltdt, subl, subldt, surfwater, refrfrac, icemeltmdt, dsnowacc, 
            hsnowmod_l, dsnowh,exposed_water,source)



@njit
def snowheight(t0, temp, water, dens, mass, ice, z, dz, lid, 
    grainsize, energy, kice, dtdz, refrfrac, cpice, precip, 
        drift,freshfrac, lidmax, rhosn, radfresh, nl, nlsnow,
        dsnowr, vink, precipsum, freshsnow,tempprecip,dsnowh):
    '''
    Routine that calculates the changes in grid or only first grid point due to addition of layers due to snow fall
    
    Parameters
    ----------
    t0 : float
        surface temperature
    temp : ndarray (1D float)
        subsurface temperature profile
    water : ndarray (1D float)
        subsurface liquid water content profile
    dens : ndarray (1D float)
        subsurface density profile
    mass : ndarray (1D float)
        subsurface mass of layer profile
    ice : ndarray (1D float)
        subsurface ice content profile
    z : ndarray (1D float)
        depth of layer centers, positive downwards
    dz : ndarray (1D float)
        layer thickness
    lid : ndarray (1D float)
        layer ID (0 = ice, 1 = snow, 2 and larger = firn)
    grainsize : ndarray (1D float)
        subsurface grainsize
    energy : ndarray (1D float)
        energy potential layer
    kice : ndarray (1D float)
        subsurface effective heat conductivity
    dtdz : ndarray (1D float)
        Vertical temperature gradient
    cpice : ndarray (1D float)
        Volumetric heat capacity of ice
    precip : float
        total precipitation during time step in m snow
    drift : float
        snow drift in mm w.e during time step
    freshfrac : ndarray (1D float)
        subsurface fresh snow fraction
    lidmax = int
        maximum subsurface layer ID (for multiyear firn)
    rhosn : float
        Density of snowfall set in routine SET_PRECIP
    radfresh : float
        fresh snow grain radius
    nl : int
        Number of subsurface layers
    nlsnow : int
        Number of snow layers
    vink : int
        switch to turn on redefining grid when layers become smaller than 0.5*dz0
    precipsum : float
        Sum of precipitation 
    freshsnow : float
        Sum of precipitation and drift
    tempprecip : float
        temporary variable keeping track of precipitation between chnages in grid
    dsnowh : float
        Thickness snow layer, thickness ice melt in m ice 
    
    Returns
    -------
    nl : int
        Number of subsurface layers
    z : ndarray (1D float)
        depth of layer centers, positive downwards
    dz : ndarray (1D float)
        layer thickness
    temp : ndarray (1D float)
        subsurface temperature profile
    dens : ndarray (1D float)
        subsurface density profile
    kice : ndarray (1D float)
        subsurface heat conductivity profile
    dtdz : ndarray (1D float)
        Vertical temperature gradient
    cpice : ndarray (1D float)
        surface heat capacity profile
    water : ndarray (1D float)
        subsurface liquid water content profile
    mass : ndarray (1D float)
        subsurface mass of layer profile
    ice : ndarray (1D float)
        subsurface ice content profile
    grainsize : ndarray (1D float)
        subsurface grainsize profile
    refrfrac : ndarray (1D float)
        subsurface profile of refrozen snow   
    sumsnow : float
        local snow layer in m snow
    summass : float
        local snow layer in m w.e.
    dsnowacc : float 
        snow layer in m snow
    hmass : float 
        snow layer in m w.e.
    vink : int
        switch to turn on redefining grid when layers become smaller than 0.5*dz0
    precipsum : float
        Sum of precipitation 
    freshsnow : float
        Sum of precipitation and drift
    tempprecip : float
        temporary variable keeping track of precipitation between chnages in grid
    dsnowr : float
        Snow height that corrects for small snow layers that are added to the ice when the layer is too thin to sustain
    dsnowh : float
        Thickness snow layer, thickness ice melt in m ice        
        
    '''


    # It is not necessary to keep track of refrfrac when loosing/adding layers because it is calculated in SNOWCONTENT
    freshsnow = freshsnow + precip + drift
    tempprecip = tempprecip + precip*rhosn
    if(abs(freshsnow) > snowthreshold):
        precipsum = precipsum + tempprecip
        tempprecip = 0.0

        nl_old = nl
        nlsnow_old = nlsnow
        mass_old = mass
        water_old = water 
        # nlsnow = nlsnow_old + NINT((precip+dsnowr)/dz0)	! add/remove layers in steps of dz0
        nlsnow = nlsnow_old + int(round((freshsnow + dsnowr)/dz0))	# add/remove layers in steps of dz0
        
        if ( (dz[0] + (freshsnow + dsnowr) < 0.5*dz0) & (nlsnow == nlsnow_old) ): nlsnow = nlsnow_old - 1
        if (nlsnow < 0): nlsnow = 0
        nlgain = 0
        nlloss = 0
        
        if (nlsnow == nlsnow_old): # Too little snow to add a new layer
            dsnowr = freshsnow+dsnowr
            mult = 1.
            if ((lid[0] > 1) | (lid[0] == 0)): mult = 0.0
            for il in range(0,nl_old):
                fact = 0.
                if (il == 0): fact=1.0
                if (tpdens == 0):
                    dz[il] = dz[il]+mult*fact*dsnowr
                    z[il] = z[il]+(1.-0.5*fact)*mult*dsnowr
                    mass[il] = dens[il]*dz[il]
                else:
                    dens_new = rhosn	# +MAX(0.,0.5*(dens_old[il]-2.*rhosn))
                    if (dsnowr < 0): dens_new = dens[il]
                    mass[il] = mass[il]+mult*fact*dsnowr*dens_new
                    dz[il] = dz[il]+mult*fact*dsnowr
                    dens[il] = mass[il]/dz[il]
                    z[il] = z[il]+(1.-0.5*fact)*mult*dsnowr
                refrfrac[il] = refrfrac[il]*((mass_old[il]+water_old[il])/(mass[il]+water[il]))
                kice[il] = conduc(dens[il],temp[il])
                cpice[il] = 152.5 + 7.122 * temp[il]
                energy[il] = dens[il]*dz[il]*cpice[il]*(Tkel-temp[il])
                if (il == 0):
                    dtdz[il] = (temp[il]-t0)/dz[il]
                else:
                    dzl=0.5*(dz[il]+dz[il-1])
                    dtdz[il] = (temp[il]-temp[il-1])/dzl
                #  energy[il] = dens[il]*dz[il]*((152.5+7.122*Tkel)*Tkel - cpice[il]*temp[il])
            # kice(1) = 50.0

            if (dsnowr > 0.0): freshfrac[0] = (dsnowr*rhosn)/(mass[0]+water[0])
            if ((tpdens > 0) & (lid[0] > lidmax) & (dsnowr < 0.)): dsnowr = 0.
            if (mult == 1): dsnowr = 0.
        else: # Save old states in _old arrays
            temp_old = temp
            water_old = water
            dens_old = dens
            mass_old = mass
            ice_old = ice
            z_old = z
            dz_old = dz
            lid_old = lid
            grainsize_old = grainsize
            energy_old = energy
            kice_old = kice
            dtdz_old = dtdz
            refrfrac_old = refrfrac
            cpice_old = cpice

            temp = np.zeros(nlmax)
            water = np.zeros(nlmax)
            dens = np.zeros(nlmax)
            mass = np.zeros(nlmax)
            ice = np.zeros(nlmax)
            z = np.zeros(nlmax)
            dz = np.zeros(nlmax)
            lid = np.zeros(nlmax)
            grainsize = np.zeros(nlmax)
            energy = np.zeros(nlmax)
            kice = np.zeros(nlmax)
            dtdz = np.zeros(nlmax)
            refrfrac = np.zeros(nlmax)
            cpice_old = np.zeros(nlmax)
            if (nlsnow > nlsnow_old):		# snow gain
                nlgain = nlsnow - nlsnow_old
                dsnowr = (freshsnow+dsnowr) - nlgain*dz0
                temp[0:nlgain] = t0
                dens[0:nlgain] = rhosn
                water[0:nlgain] = 0.
                ice[0:nlgain] = 0.
                lid[0:nlgain] = 1
                freshfrac[0:nlgain] = 1.
                grainsize[0:nlgain] = radfresh
                refrfrac[0:nlgain] = 0.
                for il in range(0,nlgain):
                    fact = 1.
                    if (il > 0): fact = 0.
                    mass[il] = rhosn*(dz0+fact*dsnowr)
                    dz[il] = dz0 + fact*dsnowr
                    z[il] = ((il+1.)-0.5)*dz0+(1.-0.5*fact)*dsnowr
                    kice[il] = conduc(dens[il],temp[il])
                    cpice[il] = 152.5 + 7.122 * temp[il]
                    #   energy[il] = dens[il]*cpice[il]*dz[il]*(Tkel - temp[il])
                    if (il == 0):
                        dtdz[il] = (temp[il]-t0)/dz[il]
                    else:
                        dzl=0.5*(dz[il]+dz[il-1])
                        dtdz[il] = (temp[il]-temp[il-1])/dzl
                    energy[il] = dens[il]*dz[il]*cpice[il]*(Tkel-temp[il])

                # kice(1) =50.0
                for il in range(nlgain,nl_old+nlgain):
                    temp[il] = temp_old[il-nlgain]
                    dens[il] = dens_old[il-nlgain]
                    energy[il] = energy_old[il-nlgain]
                    kice[il] = kice_old[il-nlgain]
                    dtdz[il] = dtdz_old[il-nlgain]
                    cpice[il] = cpice_old[il-nlgain]
                    water[il] = water_old[il-nlgain]
                    mass[il] = mass_old[il-nlgain]
                    ice[il] = ice_old[il-nlgain]
                    dz[il] = dz_old[il-nlgain]
                    z[il] = z_old[il-nlgain]+nlgain*dz0+dsnowr
                    lid[il] = lid_old[il-nlgain]
                    grainsize[il] = grainsize_old[il-nlgain]
                    refrfrac[il] = refrfrac_old[il-nlgain]

                dzl = 0.5*(dz[nlgain]+dz[nlgain-1])
                dtdz[nlgain] = (temp_old[0]-temp[nlgain-1])/dzl # Recalculate temperature gradient of first layer below new layers
                nl = nl_old + nlgain
                dsnowr = 0.
            elif (nlsnow < nlsnow_old): # snow loss
                nlloss = nlsnow_old - nlsnow
                # dsnowr = (freshsnow+dsnowr)+nlloss*dz0
                dsnowr = (freshsnow+dsnowr) + z_old[nlloss-1] + 0.5*dz_old[nlloss-1]
                mult = 1.  # NN
                if ((lid_old[nlloss] > 1) | (lid_old[nlloss] == 0)): mult = 0.0 # NN
                for il in range(0,nl_old-nlloss):
                    fact = 0.  #NN
                    if (il == 0): fact=1.0  # NN
                    if (tpdens == 0): # NN
                        dz[il] = dz_old[il+nlloss]+mult*fact*dsnowr #NN
                        z[il] = z_old[il+nlloss]+(1.-0.5*fact)*mult*dsnowr-(z_old[nlloss-1]+0.5*dz_old[nlloss-1]) #NN
                        dens[il] = dens_old[il+nlloss]
                        mass[il] = dens[il]*dz[il] #NN
                    else: #NN
                        dens_new = dens_old[il+nlloss]
                        mass[il] = mass_old[il+nlloss]+mult*fact*dsnowr*dens_new #NN
                        dz[il] = dz_old[il+nlloss]+mult*fact*dsnowr #NN
                        dens[il] = mass[il]/dz[il] #NN
                        z[il] = z_old[il+nlloss]+(1.-0.5*fact)*mult*dsnowr -(z_old[nlloss-1]+0.5*dz_old[nlloss-1]) #NN
                    temp[il] = temp_old[il+nlloss]
                    # dens[il] = dens_old[il+nlloss]
                    energy[il] = energy_old[il+nlloss]
                    kice[il] = kice_old[il+nlloss]
                    dtdz[il] = dtdz_old[il+nlloss]
                    cpice[il] = cpice_old[il+nlloss]
                    water[il] = water_old[il+nlloss]
                    #mass[il] = mass_old[il+nlloss]
                    ice[il] = ice_old[il+nlloss]
                    #dz[il] = dz_old[il+nlloss]
                    #z[il] = z_old[il+nlloss]-(z_old[nlloss]+0.5*dz_old[nlloss])
                    lid[il] = lid_old[il+nlloss]
                    grainsize[il] = grainsize_old[il+nlloss]
                    refrfrac[il] = refrfrac_old[il+nlloss]
                if ((tpdens > 0) & (lid[0] > lidmax) & (dsnowr < 0.)): dsnowr = 0.
                if (mult == 1): dsnowr = 0.
                temp[nl-nlloss-1:nl]=temp_old[nl-1]
                dens[nl-nlloss-1:nl]=dens_old[nl-1]
                energy[nl-nlloss-1:nl]=energy_old[nl-1]
                kice[nl-nlloss-1:nl]=kice_old[nl-1]
                dtdz[nl-nlloss-1:nl]=dtdz_old[nl-1]
                cpice[nl-nlloss-1:nl]=cpice_old[nl-1]
                water[nl-nlloss-1:nl]=water_old[nl-1]
                mass[nl-nlloss-1:nl] = mass_old[nl-1]
                ice[nl-nlloss-1:nl]=ice_old[nl-1]
                dz[nl-nlloss-1:nl]=dz_old[nl-1]
                z[nl-nlloss-1:nl]=z_old[nl-1]-(z_old[nlloss-1]+0.5*dz_old[nlloss-1]) 
                lid[nl-nlloss-1:nl]=lid_old[nl-1]
                grainsize[nl-nlloss-1:nl]=grainsize_old[nl-1]
                refrfrac[nl-nlloss-1:nl]=refrfrac_old[nl-1]
                nl = nl_old - nlloss
                
        freshsnow = 0.0

    if (nlsnow < .0): nlsnow = 0

    if (tpdens == 0): dens, mass, dz, irrwater = \
        resetdens(z, lid, dens, mass, dz, nl)

    # determine information for mb parameters
    sumsnow = 0.		# local snow layer in m snow
    summass = 0.		# local snow layer in m we
    topwater = 0.      # Water content in effective FAC
    sumwater = 0.		# global liquid water content snow layer in kg m-2
    dsnowacc = 0.		# global snow layer in m snow
    hmass = 0.			# global snow layer in m we
    air_content = 0.    # firn air content FAC in kg m-2
    effective_air_content = 0.    # effective / reachable firn air content in kg m-2
    lbucket = 1

    for il in range(0,nlsnow):
      if ((lid[il] > 0) & (lid[il] <= lidmax)):
        sumsnow = z[il] + 0.5*dz[il]
        summass = summass + mass[il]
        sumwater = sumwater + water[il]
        air_content = air_content + (densice - dens[il])*dz[il]/densice
        if (dens[il] >= densclosure) & (lbucket == 1):
            effective_air_content = air_content
            topwater = sumwater
            lbucket = 0.
      if (lid[il] == 1):
        dsnowacc = sumsnow
        hmass = summass
    if (sumsnow > 0.):
      dsnowh = sumsnow
    else:
      dsnowh = 0.
      dsnowacc = 0.
      hmass = 0.

    #write(99,*) dsnowh,rhosn,dz(1),dens(1)

    if (tpdens == 0):
        coeff = 0.75*(dzdeep - dz0)/zdeep + 2.
        if ((dz[0] < 0.5*dz0) | (dz[0] > coeff*dz0)): vink = 5


    if ((nlgain != 0) | (nlloss != 0)):
        if lcomment == 2:
            print('SNOWHEIGHT: Number of layers is: ',nl_old,nl,nlsnow_old,nlsnow,nlgain,nlloss)
    
    return (nl, z, dz, temp, dens, kice, dtdz, cpice,
        water, mass, ice, lid, grainsize, refrfrac, 
        sumsnow, summass, dsnowacc, hmass, vink, precipsum, 
        freshsnow, tempprecip, dsnowr, dsnowh, sumwater, topwater, air_content, effective_air_content)
        
        
@njit
def snowcontent(LE, t0, source, temp, water, ice, mass, dens, kice, cpice, grainsize, z, dz, dtdz, energy, lid, irrwater, 
                dsnowr, surfmelt, drift, sumdrift, nl, nlsnow, runoff, melt, subl, rhosn, surfwater, refrfrac):
    '''
    Routine that calculates the amount of snow melt, refreezing and percolation
    not in here yet is changes in layer density layer thickness or depth due to melt
    
    Parameters
    ----------
    LE :  float
        Latent heat flux, positive towards surface in Wm-2
    t0 : float
        surface temperature
    source : float
        Energy available for surface melt in Wm-2
    temp : ndarray (1D float)
        subsurface temperature profile
    water : ndarray (1D float)
        subsurface liquid water content profile
    ice : ndarray (1D float)
        subsurface ice content profile
    mass : ndarray (1D float)
        subsurface mass of layer profile
    dens : ndarray (1D float)
        subsurface density profile
    kice : ndarray (1D float)
        subsurface heat conductivity profile
    cpice : ndarray (1D float)
        surface heat capacity profile
    grainsize : ndarray (1D float)
        subsurface grainsize profile
    z : ndarray (1D float)
        depth of layer centers, positive downwards
    dz : ndarray (1D float)
        layer thickness
    dtdz : ndarray (1D float)
        Vertical temperature gradient
    energy : ndarray (1D float)
        Layer energy content
    lid : ndarray (1D float)
        layer ID (0 = ice, 1 = snow, 2 and larger = firn)
    irrwater : float
        irreducible water content
    dsnowr : float
        Snow height that corrects for small snow layers that are added to the ice when the layer is too thin to sustain
    surfmelt : float
        Cumulative surface melt in mm w.e.
    drift : float
        snow drift in mm w.e during time step
    sumdrift : float
        cumulative snow drift in mm w.e
    nl : int
        Number of subsurface layers
    nlsnow : int
        Number of snow layers
    runoff : float
        Surface runoff in mm w.e.
    subl : float
        Surface sublimation in mm w.e.
    rhosn : float
        Density of snowfall set in routine SET_PRECIP
    surfwater : float
        amount of water on the surface 
    refrfrac : ndarray (1D float)
        subsurface profile of refrozen snow
    
    Returns
    -------
    
    temp : ndarray (1D float)
        subsurface temperature profile
    water : ndarray (1D float)
        subsurface liquid water content profile
    ice : ndarray (1D float)
        subsurface ice content profile
    mass : ndarray (1D float)
        subsurface mass of layer profile
    dens : ndarray (1D float)
        subsurface density profile
    kice : ndarray (1D float)
        subsurface heat conductivity profile
    cpice : ndarray (1D float)
        surface heat capacity profile
    grainsize : ndarray (1D float)
        subsurface grainsize profile
    z : ndarray (1D float)
        depth of layer centers, positive downwards
    dz : ndarray (1D float)
        layer thickness
    dtdz : ndarray (1D float)
        Vertical temperature gradient
    energy : ndarray (1D float)
        Layer energy content
    lid : ndarray (1D float)
        layer ID (0 = ice, 1 = snow, 2 and larger = firn)
    dsnowr : float
        Snow height that corrects for small snow layers that are added to the ice when the layer is too thin to sustain
    surfmelt : float
        Cumulative surface melt in mm w.e.
    sumdrift : float
        cumulative snow drift in mm w.e
    nl : int
        Number of subsurface layers
    nlsnow : int
        Number of snow layers
    runoff : float
        Surface runoff in mm w.e.
    melt : float
        surface melt in mm w.e
    subl : float
        Surface sublimation in mm w.e.
    surfwater : float
        amount of water on the surface 
    refrfrac : ndarray (1D float)
        subsurface profile of refrozen snow        
    '''

    illost = 0
    ilstlost = -1
    refr = 0
    lbucket = 1
    restwater = 0.
    lmelt = source * tstep / lm		# in mm w.e.

    runoffdt = 0
    

    if ((tcalc == 2) | (tcalc == 3)): lmelt = 0.

    sumsmelt = lmelt
    surfmeltdt = lmelt
    surfmelt = surfmelt + surfmeltdt
    if ((LE > 0.) & (t0 > Tkel)):
        ll = lv  #Condensation on ice
    else:
        ll = ls # Sublimation or deposition
    cond = LE*tstep/ll				# in mm w.e.
    water[0] = water[0] + lmelt


    il = 0
    #IF (tpdens > 0) THEN	! only when density is allowed to change

    if ((cond > 0.)  & (t0 >= Tkel)): # condensation, heat/water added to the surface
        water[0]= water[0] + cond
    else:  #IF (cond < 0., or >0 with T<Tkel) 
        if (cond > 0): # addition of ice through deposition	
            mass[0] = mass[0] + cond 
            if (tpdens > 0):
                dz[0] = dz[0] + cond/densice
                dens[0] = mass[0]/dz[0]
            else:       
                dz[0] = mass[0]/dens[0]
            # print(dz[0],' m remaining in first layer after deposition')
            cpice[0] = 152.5 + 7.122*temp[0]
            z[0] = 0.5*dz[0]
            dtdz[0] = (temp[0]-t0)/dz[0]
            energy[0] = dens[0]*dz[0]*cpice[0]*(Tkel-temp[0])
        else:   # removal of ice through sublimation
            surfen = cond 
        #    il = 1
            while (surfen < 0.):		
                mass[il] = mass[il] + surfen
                if (mass[il] < 0.):
                    surfen = mass[il]
                    mass[il] = 0.
                    dz[il] = 0.
                    water[il] = water[il] + surfen
                    energy[il] = 0.
                    z[il] = 0.
                    if (water[il] < 0.):
                        water[il] = 0.
                        surfen = water[il]
                    else:
                        water[il+1] = water[il+1] + water[il]	# not sure what to do with the water, maybe also remove with cond?
                        water[il] = 0.
                    il = il+1
                else: 
                    dz[il] = mass[il]/dens[il]
                    # print(dz[0],'m remaining in first layer after sublimation')
                    cpice[il] = 152.5 + 7.122*temp[il]
                    z[il] = 0.5*dz[il]
                    dtdz[il] = (temp[il]-t0)/dz[il]
                    if (il == nl-1):
                        energy[il] = energy[il-1]
                    else:
                        energy[il] = dens[il]*dz[il]*cpice[il]*(Tkel-temp[il])
                    surfen = 0.
        if ((dz[0] <= 0) & (lcomment == 1)):
            print('SNOWCONTENT: dz[0] negative!!! ',dz[0],cond,mass[0],drift,il)

        #    dz[0] = dz[0] + cond/dens[0]
        # IF ((dz[0] <= 0) .and. (lcomment == 1)) THEN
            # WRITE(*,'(/,a,4f15.8,i3,/)') 'SNOWCONTENT: dz[0] negative!!! ',dz[0],cond,mass[0],drift,il
    sumdrift = sumdrift + drift		# in case luseacc = 0, drift = 0
    # sumdrift = sumdrift + sbuf(20)		! in case luseacc = 0, drift = 0
    # surfen = -lmelt - drift*rhosn
    surfen = -lmelt 

    if ((surfen < 0.) & (dsnowr > 0.)):
        dsnowr = dsnowr + surfen/rhosn
        if (dsnowr < 0.): dsnowr = 0.

    # il = 1
    while (surfen < 0.):		# removal of mass through melt
        mass[il] = mass[il] + surfen
        if (mass[il] < 0.):
            surfen = mass[il]
            mass[il] = 0.
            dz[il] = 0.
            energy[il] = 0.
            z[il] = 0.
            water[il+1] = water[il+1] + water[il]
            water[il] = 0.
            il = il+1
        else:
            dz[il] = mass[il]/dens[il]
            cpice[il] = 152.5 + 7.122*temp[il]
            z[il] = 0.5*dz[il]
            dtdz[il] = (temp[il]-t0)/dz[il]
            # print(dz[0],'m remaining in first layer after melt')
            if (il == nl-1):
                energy[il] = energy[il-1]
            else:
                energy[il] = dens[il]*dz[il]*cpice[il]*(Tkel-temp[il])
            surfen = 0.

    if (il > 0):
        illost = il
        if lcomment == 1:
            print('SNOWCONTENT: ',il,' Layer(s) completely melted away or removed by drift')
        for il in range(0,nl-illost):
            temp[il] = temp[il+illost]
            water[il] = water[il+illost]
            ice[il] = ice[il+illost]
            mass[il] = mass[il+illost]
            dens[il] = dens[il+illost]
            kice[il] = kice[il+illost]
            cpice[il] = cpice[il+illost]
            grainsize[il] = grainsize[il+illost]
            dz[il] = dz[il+illost]
            if (il == 0):
                z[il] = 0.5*dz[il]
                dtdz[il] = (temp[il]-t0)/dz[il]
            else:
                z[il] = 0.5*dz[il] + z[il-1]
                dzl = 0.5*(dz[il]+dz[il-1])
                dtdz[il] = (temp[il]-temp[il-1])/dzl
            energy[il] = dens[il]*dz[il]*cpice[il]*(Tkel-temp[il])
            lid[il] = lid[il+illost]
        nl = nl - illost
        nlsnow = nlsnow - illost
        if (nlsnow < 0): nlsnow = 0
        illost = 0
    
    if (lrefr == 0):
        for il in range(0,nl):
            if (il == 0): water[il] = water[il] + surfwater 
            runoffdt = runoffdt + water[il]
            water[il] = 0.
    else:       
        refrfrac = np.zeros(nlmax)
        for il in range(0,nl):
            # first calculate the amount of melt or refreezing depending on layer temperature
            # energy >0 in case T<Tkel, refreezing
            # else energy < 0, melt
            if (il == 0):
                water[il] = water[il] + surfwater
            else:
                if (dens[il]>=densclosure): lbucket = 0.
                if (lbucket == 1):
                    water[il] = water[il] + restwater #Add water that percolates from layer(s) above
                    restwater = 0.
                else:
                    water[il] = water[il] # No percolation allowed though ice layer, all remaining water stored in "restwater"
            if (temp[il] >= Tkel):	# melt
                water[il] = water[il] - energy[il]/lm
                sumsmelt = sumsmelt - energy[il]/lm
                if (tpdens > 0): 
                    mass[il] = mass[il] + energy[il]/lm
                    if (mass[il] < 0.):
                        mass[il+1] = mass[il+1]+mass[il]
                        water[il+1]=water[il+1]+water[il]
                        mass[il] = 0.
                        water[il] = 0.
                        illost = illost + 1  
                        if (ilstlost < 0): ilstlost = il        
                        if (lcomment == 1):
                            print('SNOWCONTENT: ',il,illost,' Layer(s) completely melted away')
                temp[il] = Tkel
            elif ((temp[il] < Tkel) & (water[il] > 0.) & (dens[il] < densice)):	# refreezing
                maxwater = denswater * dz[il]* (densice-dens[il])/densice
                refrwater = water[il]
                refr = 1
                if (maxwater < water[il]): refrwater = maxwater
                if (energy[il] > lm*refrwater):	# refreeze all water
                    ice[il] = ice[il] + refrwater
                    water[il] = water[il] - refrwater
                    temp[il] = temp[il] + lm*refrwater/(dens[il]*cpice[il]*dz[il])
                #      energy[il] = energy[il]-(lm*refrwater)
                #      temp_old = temp[il]
                    if (tpdens > 0): 
                        mass[il] = dens[il]*dz[il] + refrwater
                        dens[il] = mass[il]/dz[il]
                #        temp[il] = solvetemp(temp_old,energy[il],mass[il])
                        energy[il] = dens[il]*dz[il]*cpice[il]*(Tkel-temp[il])
                        cpice[il] = 152.5+7.122*temp[il]
                        if (dens[il] > densice):
                            dens[il] = densice
                            restwater = restwater + water[il] 
                            water[il] = 0.
                        kice[il] = conduc(dens[il],temp[il])
                        if (il == 0):
                            dtdz[il] = (temp[il]-t0)/dz[il]
                        else:
                            dzl=0.5*(dz[il]+dz[il-1])
                            dtdz[il]=(temp[il]-temp[il-1])/dzl
                    else:
                        energy[il] = dens[il]*dz[il]*cpice[il]*(Tkel-temp[il])
                        kice[il] = conduc(dens[il],temp[il])
                        cpice[il] = 152.5+7.122*temp[il]
                        if (il == 0):
                            dtdz[il] = (temp[il]-t0)/dz[il]
                        else:
                            dzl=0.5*(dz[il]+dz[il-1])
                            dtdz[il]=(temp[il]-temp[il-1])/dzl
                else:									# refreeze till energy content is gone
                    refrwater = energy[il]/lm
                    energy[il] = 0.
                    water[il] = water[il] - refrwater
                    ice[il] = ice[il] + refrwater
                    if (tpdens > 0): 
                        mass[il] = dens[il]*dz[il] + refrwater
                        dens[il] = mass[il]/dz[il]
                    temp[il] = Tkel
                    
            elif (dens[il] >= densice):
                restwater = restwater + water[il]
                if (ice[il] > mass[il]): ice[il] = mass[il]
                if (lid[il] == 0): ice[il] = 0.
                water[il] = 0.
            
            # Then calculate the amount of percolating water
            maxwater = denswater * irrwater[il] * dz[il]* (densice-dens[il])/densice
            if (dens[il] >= (densice - 0.001*densice)): maxwater = 0.
            if (water[il] > maxwater):			# too much water in a layer
                restwater = restwater + water[il] - maxwater
                water[il] = maxwater

            if ((tpdens > 0) & (dens[il] < densice)):
                dens[il] = densification(il,nl,dens,temp,z,dz,t0)
            
            if (refr == 1): # Refreezing has taken place
                if (lslush == 0):
                    refrfrac[il] = refrwater/(mass[il]+water[il])
                    refr = 0
                else:
                    refrfrac[il] = refrwater/mass[il] #Water will be changed by slush routine, do it this way to store data already, and multiply by mass[il]/(mass[il]+water[il]) in slush routine to get it done

    slushdepth = 0.
    if ((lslush == 1) & (restwater > 0.)): 
        slushdepth, surfwater, restwater, refrfrac, water, refrbool = slush(restwater,refr,refrfrac,lid, water, mass, dens, dz, nl)

    runoffdt =  runoffdt + restwater
    runoff = runoff + runoffdt

    # MvT start melt lake routine at this point if runoff > 0 but slope = 0
    if (runoffdt > 0) & (surfangle == 0 ) & (llake == 1): 
        exposed_water = 1
        print('exposed water')
    meltdt = sumsmelt
    melt = melt + meltdt
    if ((cond > 0.) & (t0 >= Tkel)):
        subldt = 0.
    else:
        subldt = cond
    subl = subl + subldt

    if (illost > 0):
        if (lcomment == 1):
            print('SNOWCONTENT: adapt grid when layers are lost',ilstlost,nl,illost)

        for il in range(ilstlost,nl-illost):
            temp[il] = temp[il+illost]
            water[il] = water[il+illost]
            ice[il] = ice[il+illost]
            mass[il] = mass[il+illost]
            dens[il] = dens[il+illost]
            kice[il] = kice[il+illost]
            cpice[il] = cpice[il+illost]
            grainsize[il] = grainsize[il+illost]
            dz[il] = dz[il+illost]
            if (il == 0):
                z[il] = 0.5*dz[il]
                dtdz[il] = (temp[il]-t0)/dz[il]
            else:
                z[il] = 0.5*dz[il] + z[il-1]
                dzl=0.5*(dz[il]+dz[il-1])
                dtdz[il] = (temp[il]-temp[il-1])/dzl
            energy[il] = dens[il]*dz[il]*cpice[il]*(Tkel-temp[il])
            lid[il] = lid[il+illost]
        nl = nl-illost
        if (ilstlost <= nlsnow): nlsnow = nlsnow - illost
        if (nlsnow < 0): nlsnow = 0

        
    # MvT TODO if ( water[0] < -99.) STOP 20

    return (temp, water, ice, mass, dens, kice, cpice, grainsize, z, dz, dtdz, 
            energy, lid, dsnowr, surfmelt, surfmeltdt, sumdrift, nl, nlsnow, runoff, runoffdt,
            melt, meltdt, subl, subldt, surfwater, refrfrac,exposed_water)


@njit
def slush(restwater,refrbool,refrfrac,lid, water, mass, dens, dz, nl):
    '''
    Routine that calculates the time scale of runoff
    based on Zuo and Oerlemans, 1996
    
    Parameters
    ----------
    restwater : float
        amount of remaining liquid water after percolation
        layer index
    refrbool : int
        flag = 1 if refreezing has taken place in SNOWCONTENT routine
    refrfrac : ndarray (1D float)
        subsurface profile of refrozen snow 
    lid : ndarray (1D float)
        layer ID (0 = ice, 1 = snow, 2 and larger = firn)   
    water : ndarray (1D float)
        subsurface liquid water content profile
    mass : ndarray (1D float)
        subsurface mass of layer profile
    dens : ndarray (1D float)
        subsurface density profile
    dz : ndarray (1D float)
        layer thickness
    Returns
    -------
    slushdepth : float
        depth of the water saturated snow layer
    surfwater : float
        amount of water on the surface 
    restwater : float
        amount of remaining liquid water after percolation
    refrfrac : ndarray (1D float)
        subsurface profile of refrozen snow
    water : ndarray (1D float)
        subsurface liquid water content profile
    refrbool : int
        flag = 1 if refreezing has taken place in SNOWCONTENT routine
        
    '''
    
    cdouble = 1.
    c1 = tausteep
    c2 = tauhor - tausteep
    c3 = -np.log( (tau1 - tausteep)/(tauhor - tausteep))/np.tan(cdouble*np.pi/180)

    cdouble = surfangle

    timestar = c1+c2*np.exp(c3*np.tan(cdouble*np.pi/180))		# time scale in days
    timefactsurf = (tstep/(24.*3600))/timestar		# timescale per timestep for surface runoff, fraction to runoff per time step
    if (timefactsurf > 1.): timefactsurf = 1.				# 1 = immediate runoff everything
    timefactint = timefactsurf/slfact							# timescale per timestep for runoff in snow, slower than on surface, fraction to runoff per time step

    surplus = restwater*(1.-timefactint)
    restwater = restwater * timefactint
    slushdepth = 0.
    surfwater = 0.
    for il in range(nl,1,-1):
        if ((lid[il] > 0) & (dens[il] < densice-0.001*densice)): # in snow/firn, water cannot be retained in ice 
            maxwater = denswater * dz[il]* (densice-dens[il])/densice
            if (water[il] < maxwater):
                if (surplus > (maxwater-water[il])):
                    surplus = surplus - maxwater + water[il]
                    water[il] = maxwater
                    if (surplus < 0.): surplus = 0.
                    slushdepth = slushdepth + dz[il]
                else:
                    water[il] = water[il] + surplus
                    surplus = 0.
                    slushdepth = slushdepth + (water[il]/maxwater)*dz[il]
        if (refrbool == 1): # Refreezing has occurred in snowcontent subroutine
            refrfrac[il] = refrfrac[il]*(mass[il]/(mass[il]+water[il]))
            refrbool = 0
    if (surplus > 0.):	# water on ice
        surfwater = surplus*(1.-timefactsurf)
        restwater = restwater + surplus*timefactsurf

    return slushdepth, surfwater, restwater, refrfrac, water, refrbool

@njit
def densification(il,nl,dens,temp,z,dz,t0):
    '''
    Routine that calculates the densification of the snowpack due to gravity and
    vapour transport see Herron and Langway (1980), Li and Zwally (2004), Colbeck (1993),
    Helsen et al. (2008), Arthern et al (2010), Ligtenberg et al (2011).
    
    Parameters
    ----------
    il : int
        layer index
    nl : int
        number of subsurface layers
    dens : ndarray (1D float)
        subsurface density profile    
    temp : ndarray (1D float)
        subsurface temperature profile    
    z : ndarray (1D float)
        depth of layer centers, positive downwards
    dz : ndarray (1D float)
        layer thickness
    t0 : float
        surface temperature
        
    Returns
    -------
    dens: float
        density of layer at index il
        
    '''


    dens_o = dens[il]

    if ((il > 0) & (il < nl-1)):
        temp1 = temp[il-1]
        temp2 = temp[il]
        temp3 = temp[il+1]
        z1 = z[il-1]
        z2 = z[il]
        z3 = z[il+1]
    elif (il == 0):
        temp1 = t0
        temp2 = temp[il]
        temp3 = temp[il+1]
        z1 = 0.
        z2 = z[il]
        z3 = z[il+1]
    elif (il == nl-1):
        temp1 = temp[il-1]
        temp2 = temp[il]
        temp3 = temp[il]
        z1 = z[il-1]
        z2 = z[il]
        z3 = zdeep
        
    # input accyear is in m w.e.
    coeftime = 1./(365.25*24.*3600.)		# factor converting accumulation rate per year to per second

    if (tpdens == 2): # Herron and Langway 1980
        if (dens_o <= 550.):
            Ec = 10160. 
            K0 = 11.
            alpha = 1.
        else:
            Ec = 21400. 
            K0 = 575.
            alpha = 0.5
        coef1 = -Ec/(rstar*temp2) 
        coef2 = accyear**alpha
        coef3 = (densice - dens_o)		# /densice
        grav = K0*np.exp(coef1)*coef2*coef3*coeftime
        drhodt =  grav 

    elif ((tpdens == 3) | (tpdens == 4)):	# Li and Zwally 2004 (as described in Arthern et al. 2004)
        templim = temp[il]
        if (templim > tempcut): templim = tempcut
        Ec = 883.8*((np.abs(templim - Tkel))**(-0.885))			#KJ/mol *1000. problem when you multiply with 1000
        K0G = 8.36*((np.abs(templim - Tkel))**(-2.061))		#mm2/yr according to article but I think mm2/kg
        beta = 139.21-0.542*T10m
        if (beta < betadens): beta = betadens
        K0 = beta*K0G									
        coef1 = -Ec/(rstar*temp2) 
        coef2 = accyear
        coef3 = (densice - dens_o)		#/densice
        grav = K0*np.exp(coef1)*coef2*coef3*coeftime
        drhodt =  grav 
        if (tpdens == 4):							# adition of effects of vapor transport Li and Zwally 2004
            coef1 = es0/(rv*rv)
            coef2 = (ls-rv*temp1)/(temp1**3)
            coef3 = np.exp((ls/rv)*((1./Tkel)-(1./temp1)))
            coef4 = (temp2-temp1)/(z2-z1)
            J1 = -Deff*coef1*coef2*coef3*coef4
            coef2 = (ls-rv*temp2)/(temp2**3)
            coef3 = np.exp((ls/rv)*((1./Tkel)-(1./temp2)))
            coef4 = (temp3-temp2)/(z3-z2)
            J2 = -Deff*coef1*coef2*coef3*coef4
            dJdz = -(J2-J1)/dz[il]
            drhodt =  drhodt + dJdz 
            if (drhodt < 0. ): drhodt = 0.

    elif (tpdens == 5): # Helsen 2008
        templim = temp[il]
        if (templim > tempcut): templim = tempcut
        Ec = 883.8*((np.abs(templim - Tkel))**(-0.885)) # KJ/mol *1000. problem when you multiply with 1000
        K0G = 8.36*((np.abs(templim - Tkel))**(-2.061)) # mm2/yr according to article but I think mm2/kg
        beta = 76.138-0.28965*T10m
        if (beta < betadens): beta = betadens
        K0 = beta*K0G									
        coef1 = -Ec/(rstar*temp2) 
        coef2 = accyear
        coef3 = (densice - dens_o)		# /densice
        coef4 = 0.0
        grav = K0*np.exp(coef1)*coef2*coef3*coeftime
        drhodt =  grav 

    elif ((tpdens == 6) | (tpdens == 7)): # Arthern et al. 2010, Ligtenberg 2011
        Ec =  60.*1000.				# 60 kJ/mol
        Eg = 42.4*1000.				# 42.4 kJ/mol
        if (dens_o <= 550.):
            K0 = 0.07
            coef4 = 1.415-0.147*np.log(accyear*1000.)
        else:
            K0 = 0.03
            coef4 = 2.335-0.288*np.log(accyear*1000.)
        if (coef4 < 0.25): coef4 = 0.25
        if (tpdens == 6): coef4 = 1. # Arthern
        coef1 = -(Ec/(rstar*temp2)) + (Eg/(rstar*T10m))
        coef2 = g*accyear*1000.	 # accyear must be in mm per year
        coef3 = (densice - dens_o)	# /densice
        grav = K0*np.exp(coef1)*coef2*coef3*coef4*coeftime
        drhodt =  grav

    dens[il] = drhodt*tstep + dens_o

    if (dens[il] > densice-accur): dens[il] = densice
    if ((dens[il] < densnow) | (dens[il] > densice)):
        print("DENSIFICATION error",il,dens[il],dens_o,grav,densnow)
        # raise ValueError("DENSIFICATION error",il,dens[il],dens_o,coef1,coef2,coef3,coef4,grav,densnow)

    return dens[il]