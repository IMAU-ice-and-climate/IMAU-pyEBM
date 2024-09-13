#===============================================================================
#
# All routines related to the calculation of the turbuent heat fluxes
#
#===============================================================================
from numba import njit
import numpy as np

from globals import *
from info import *
from routines import * 

@njit
def energybalance(t0,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf):
    #===============================================================================
    # This function calculates the energy balance as a function of the surface temperature (K)

    # ! input
    # REAL :: t0
    # ! local
    # REAL :: l, term1, term2, term3, term4, term4a, term4b
    # REAL :: energybalance, spechum
    
    zt = sbuf[13]

    if (t0 < Tkel): l = ls
    if (t0 >= Tkel): l = lv

    q0 = spechum(t0, 1., sbuf[4])

    # net radiation
    term1 = Snet - sumdivs + (sbuf[10]-emis*StefBoltz*(np.power(t0,4.)))
    # sensible heat flux
    term2 = Ch*densair*cp*sbuf[5]*(sbuf[3]+(zt*g/cp)-t0)
    # latent heat flux
    term3 = Cq*densair*l*sbuf[5]*(sbuf[9]-q0)

    # subsurface energy flux
    #kice(1) = 50.0 ! ECMWF trick for glaciers and ice caps
    if (extrapolation == 1):
        term4 = -(kice[0]*(t0-temp[0]))/dz[0]
        #  term4 = -(kice(1)*(t0-temp(1)))/(0.5*dz(1))
    elif (extrapolation == 2):
        term4a = kice[0]*(t0-temp[0])/dz[0]
        #  term4a = kice(1)*(t0-temp(1))/(0.5*dz(1))
        term4b = ((kice[0]*dz[0]+kice[1]*dz[1])/(dz[1]+dz[1]))*((temp[0]-temp[1])/(0.5*(dz[0]+dz[1])))
        term4 = -(dz[0]*(2.0*term4a - term4b) + dz[1]*term4a)/(dz[0]+dz[1])

    if (((tcalc == 2) | (tcalc == 3)) & (t0 < Tkel)):
        GH = -(term1 + term2 + term3)
    else:
        GH = term4

    if ((tcalc == 2) | (tcalc == 3)): term4 = 0.0
    #IF (t0>=Tkel) term4 = 0.

    return (term1 + term2 + term3 + term4, GH)

@njit
def turbhf(skiter,t0,sbuf):
    # USE GLOBALS_EBM , ONLY : sbuf, t0, q0, SH, LE, densair, Ch, Cq, zt, &
    # &                        zm, t2m, q2m, ws10m , ustar,thstar,qstar, psim, psih, psiq , &
    # &                        psim0, psih0, psiq0, z0q,z0h,z0m, Hice
    # USE CONSTANTS_EBM , ONLY : lv, ls, karman, g, cp, rd, eps, Lmocrit , Tkel
    # USE INPUT_EBM , ONLY : lcomment, lz0m, Hmax, lz0h
    # !USE SNOW_EBM , ONLY : dsnowacc, melt
    # USE FILES , ONLY : uo1

    # IMPLICIT NONE
    # !input
    # INTEGER, INTENT(IN) :: skiter
    # !Local
    # INTEGER :: iter, itermax
    # INTEGER :: crit
    # INTEGER :: landreas
    # REAL :: t,q,ws,p,theta
    # REAL :: l
    # REAL :: Lmo,Lmo_old
    # REAL :: psimin,z0hmin,z0qmin
    # REAL :: Chn,Cqn

    if (t0 < Tkel): l = ls
    if (t0 >= Tkel): l = lv

    t = sbuf[3]
    q = sbuf[9]
    ws = sbuf[5]
    p = sbuf[4]
    zm = sbuf[14]
    zt = sbuf[13]
    z0m = sbuf[12]
    


    densair = p / (rd*t)
    q0 = spechum(t0, 1.0, sbuf[4])

    # First step, calculate u*, th*, q* under neutral conditions
    ustar = (karman * (ws - 0.0 )) / np.log(zm/z0m)

    if (lz0h == 0):
        z0h = z0m*0.1
        z0q = z0m*0.1
    else:
        z0h, z0q = z0h_model(ustar, densair, t, z0m)

    theta = t+zt*g/cp

    thstar = (karman * (t - t0)) / np.log(zt/z0h)
    qstar = (karman * (q - q0)) / np.log(zt/z0q)

    Lmo = 10.E4
    Lmo_old = 1.E4
    psimin = -2.0
    z0hmin = 1.0E-10
    z0qmin = 1.0E-10
    itermax = 40
    iter = 0

    psim = 0.0
    psih = 0.0
    psiq = 0.0
    psim0 = 0.0
    psih0 = 0.0
    psiq0 = 0.0

    if (theta > t0): crit = 1 # stable conditions
    if (theta <= t0): crit = -1 # unstable conditions

    while ((np.abs((Lmo_old-Lmo)/Lmo_old) > Lmocrit) & (iter < itermax)):
        iter = iter + 1
        # Now add stability and iterate
        psim, psih, psiq = stabilityfunc(crit, Lmo, zm, zt, zt)
        psim0, psih0, psiq0 = stabilityfunc(crit, Lmo, z0m, z0h, z0q)
        if (crit == 1):	# Limited stability correction
            if (Lmo == 0.): iter = iter+itermax
            if (psim < psimin): psim = psimin		
            if (psih < psimin): psih = psimin
            if (psiq < psimin): psiq = psimin
            if (psim0 < psimin): psim0 = psimin		
            if (psih0 < psimin): psih0 = psimin
            if (psiq0 < psimin): psiq0 = psimin

        # Recalculate the u*, th* and q*
        ustar = (karman * (ws - 0.0 )) / ( np.log(zm/z0m) - psim + psim0 )
        if ((crit == 1) & (ustar < 0.)): ustar = (karman * (ws - 0.0 )) / np.log(zm/z0m)

        if (lz0h == 0):
            z0h = z0m*0.1
            z0q = z0m*0.1
        else:
            z0h, z0q = z0h_model(ustar, densair, t, z0m)

        if (z0h < z0hmin): z0h = z0hmin
        if (z0q < z0qmin): z0q = z0qmin

        thstar = (karman * (theta - t0)) / ( np.log(zt/z0h) - psih + psih0 )
        qstar = (karman * (q - q0)) / ( np.log(zt/z0q) - psiq + psiq0 )

        Lmo_old = Lmo
        Lmo = (ustar**2) / ( (karman*g/t) * (thstar + t * ((1.-eps)/eps) * qstar ) )
        crit = 1 # Stable
        if (Lmo < 0.): crit = -1 # Unstable

    if (iter >= itermax): # no solution found
        if (lcomment == 1):
            print('TURBHF more than itermax iterations necessary')
        #    WRITE(*,'(/,A,i3,A,/,2I5,16f16.8,/)') 'TURBHF more than',itermax,' iterations necessary',&
        #       skiter,iter,Lmo,Lmocrit,ABS((Lmo_old-Lmo)/Lmo_old),psim,psih,psim0,psih0,t,t0,ws,z0m,z0h,z0q,&
        #       ustar,thstar,qstar
        #   WRITE(*,'(/,A,i3,A,/,4f16.8,/)') 'TURBHF more than',itermax,' iterations necessary',&
        #&       Lmo,Lmocrit,ABS((Lmo_old-Lmo)/Lmo_old),crit*1.0
        #  STOP
        #WRITE(uo1,'(/,A,i3,A,/,2I5,16f16.8,/)') 'TURBHF  more than',itermax,' iterations necessary',&
        #&       skiter,iter,Lmo,Lmocrit,ABS((Lmo_old-Lmo)/Lmo_old),crit*1.0!,psim,psih,psim0,psih0,t,t0,ws,z0m,z0h,z0q,&
        #!&       ustar,thstar,qstar
        SH = 0.
        LE = 0.
        Ch = 0.
        Cq = 0.
        Chn = 0.
        Cqn = 0.
        psim = 0.
        psih = 0.
        psiq = 0.
        psim0 = 0.
        psih0 = 0.
        psiq0 = 0.
    elif (Lmo == 0.):
        SH = 0.
        LE = 0.
        Ch = 0.
        Cq = 0.
        Chn = 0.
        Cqn = 0.
        psim = 0.
        psih = 0.
        psiq = 0.
        psim0 = 0.
        psih0 = 0.
        psiq0 = 0.
    else:
        SH = densair*cp*thstar*ustar
        LE = densair*l*qstar*ustar
        Ch = (karman**2.)/((np.log(zm/z0m)-psim+psim0)*(np.log(zt/z0h)-psih+psih0))
        Cq = (karman**2.)/((np.log(zm/z0m)-psim+psim0)*(np.log(zt/z0q)-psiq+psiq0))
        Chn = (karman**2.)/(np.log(zm/z0m)*np.log(zt/z0h))
        Cqn = (karman**2.)/(np.log(zm/z0m)*np.log(zt/z0q))
        
    t2m = t0+thstar/karman*(np.log(2./z0h)-psih+psih0)-2*g/cp
    q2m = q0+qstar/karman*(np.log(2./z0q)-psiq+psiq0)
    ws10m = ustar/karman*(np.log(10./z0m)-psim+psim0)
    if (ws10m<0.1): ws10m = 0.1		# can become negative for very low wind speeds

    return ustar, SH, LE, Ch, Cq, Chn, Cqn, t2m, q0, q2m, ws10m, densair

@njit
def z0h_model(ustar, densair, temp, z0m):
    # !===============================================================================
    # !Calculate roughness length for heat and moisture, necessary to calculate th* and q* according to Andreas
    # !===============================================================================
    # USE GLOBALS_EBM , ONLY : z0m
    # USE INPUT_EBM , ONLY : chstation, lz0h

    # IMPLICIT NONE
    # !Input
    # REAL, INTENT(IN) :: ustar, densair, temp
    # !Output
    # REAL, INTENT(OUT) :: z0h, z0q
    # !Local
    # REAL :: Restar , visc, mu
    # REAL :: z0_min

    z0_min = 1.0E-10

    mu = (-0.039225 + 0.0079067 * temp - 5.1515E-6 * (np.power(temp,2))) / 1.0E5
    # visc = the kinematic fluid viscosity of air (m2/s) for a given pressure (Pa) and temperature (K)
    visc = mu/densair
    Restar = ustar*z0m/visc

    if (Restar <= 0.135): # Smooth regime
        z0h = z0m * np.exp( 1.250 )
        z0q = z0m * np.exp( 1.610 )
    elif ((Restar > 0.135) & (Restar < 2.5)): # Transition regime
        z0h = z0m * np.exp( 0.149 - 0.550*np.log(Restar) )
        z0q = z0m * np.exp( 0.351 - 0.628*np.log(Restar) )
    elif (Restar >= 2.5): # Rough Regime
        z0h = z0m * np.exp( 0.317 - 0.565*np.log(Restar) - 0.183*np.power(np.log(Restar),2) )
        z0q = z0m * np.exp( 0.396 - 0.512*np.log(Restar) - 0.180*np.power(np.log(Restar),2) )
    if ((z0m >= 0.001) & (lz0h == 2)): # Extension to very rough terrain by Smeets and VdBroeke 2008.
        z0h = z0m * np.exp( 1.5 - 0.2*np.log(Restar) - 0.11*np.power(np.log(Restar),2) ) 
        z0q = z0h
    elif ((z0m >= 0.001) & (lz0h == 3)): # Updated parameterization from Van Tiggelen et al (2021) 
        z0h = z0m * np.exp( 1.5 - 0.15*np.log(Restar) - 0.16*np.power(np.log(Restar),2) ) 
        z0q = z0h

    if (z0h < z0_min): z0h = z0_min
    if (z0q < z0_min): z0q = z0_min
    
    return z0h, z0q

@njit
def stabilityfunc(crit, Lmo, hm, hh, hq):
    # USE CONSTANTS_EBM , ONLY : pi
    # IMPLICIT NONE
    # !Input
    # INTEGER :: crit
    # REAL :: Lmo, hm, hh, hq
    # !Output
    # REAL :: psim, psih, psiq
    # !Local
    # REAL :: fim, fih, fiq
    # REAL :: psim_min

    psim_min = -2.0
    
    # Check the used functions!!!!
    if (crit == -1): # unstable conditions Dyer 1974
    # fim = (1.0 - 20.0 * hm / Lmo )**0.25		!momentum
    # fih = (1.0 - 15.0 * hh / Lmo )**0.5		!temperature
    # fiq = (1.0 - 15.0 * hq / Lmo )**0.5		!humidity
    #  Functions of Dyer 1974
        fim = np.power((1.0 - 16.0 * hm / Lmo ),(0.25)) # momentum
        fih = np.power((1.0 - 16.0 * hh / Lmo ),(0.5)) # temperature
        fiq = np.power((1.0 - 16.0 * hq / Lmo ),(0.5)) # humidity
        psim = 2. * np.log( (1.+fim)/2. ) + np.log( (1.+np.power(fim,2))/2. ) - 2.0*np.arctan(fim) + np.pi/2.0
        psih = 2. * np.log( (1.+fih)/2. )
        psiq = 2. * np.log( (1.+fiq)/2. )
    elif (crit == 1): # stable conditions Holtslag en de Bruin 1988
        # psim = -( 1 + 6.25* hm / Lmo )**0.8		!momentum
        # psih = -( 1 + 9.375* hh / Lmo )**0.8		!temperature
        # psiq = -( 1 + 9.375* hq / Lmo )**0.8		!humidity
        # Functions of Holtslag and de Bruijn 1988
        psim = -( (0.7 * hm / Lmo) + 0.75 * ((hm / Lmo) - (5. / 0.35)) * np.exp(-0.35 * (hm / Lmo)) + 0.75 * 5. / 0.35 ) # momentum
        if (psim < psim_min): psim = psim_min
        psih = psim	 # temperature
        psiq = psim	# humidity

    return psim, psih, psiq

@njit
def get_z0m(sbuf,lid,rndnr, dsnowacc, ablation, Hice, Hmax):

    z0m = sbuf[12] # by default, z0m is given in forcing file
    dHice = 0.
    if (lz0m == 0): # prescribed 
        return z0m, Hice, dHice
    elif ((lz0m == 1) | (lz0m == 2) | (lz0m == 3)): # observed surface type determines roughness
        if (lz0m == 2): # z0m snow is random
            z0msn_tmp = np.power(10,(rndnr[0]*(np.log10(zul)-np.log10(zll)))+np.log10(zll))
            z0mice_tmp = z0mice
        elif (lz0m == 3): # z0m ice is random
            z0mice_tmp = np.power(10,(rndnr[0]*(np.log10(zul)-np.log10(zll)))+np.log10(zll))
            z0msn_tmp = z0msn
        else:
            z0mice_tmp = z0mice
            z0msn_tmp = z0msn
        z0m = z0msn_tmp
        if (luseacc >= 2): # if we can use height ranger data
            if (sbuf[17] < 0.005): # ice
                z0m = z0mice_tmp
        elif (lid[0] == 0): # or lack of more information, must be based on model surface type
            z0m = z0mice_tmp
    elif (lz0m == 10): # modelled surface type determines roughness
        if (lid[0] == 0): z0m = z0mice # ice
        else: z0m = z0msn # snow
    elif (lz0m == 4): # parameterized roughness in time
        z0m, Hice, dHice = z0m_model(dsnowacc, ablation, Hmax, Hice)
    else:
        return z0m, Hice, dHice
        
    return z0m, Hice, dHice
 
@njit
def z0m_model(dsnowacc, ablation, Hmax, Hice):
    """
    Calculate roughness length for momentum based on Raupach (1994) and Van Tiggelen et al. (2021) 

    Parameters
    ----------
    dsnowacc : float
        snow height at current timestep
    ablation : float
        ice ablation at previous timestep
    Hmax : float
        maximum height of ice obstacles
    Hice : float
        height of ice obstacles at previous timestep
        

    Returns
    -------
    z0m: float 
        roughness length for momentum at current timestep
    Hice: float 
        height of ice obstacles at current timestep
    dHice: float
        change of height of ice obstacles at current timestep
    """

    # Parameters
    Hmin = Hmax / 2    # minimum height of ice obstacles [m]
    dHsub = 1e-3       # obstacle height reduction due to sublimation [m/day]
    dHabl = 0.1        # fraction of ice ablation contributing to increase of ice obstacle height [-]
    f = 8              # Obstacles per 100 m [#]
    Cs10 = 1.2071e-3   # 10m drag coefficient for skin friction
    c1 = 7.5           # empirical constant for displacement height
    PsiH = 0.193       # roughness layer correction for wind profile
    
    # Model for obstacle height
    Hsnow = dsnowacc 

    if (ablation > .0):
        dHice = dHabl * ablation 
    else:
        dHice = - dHsub * (float(tstep)/(3600*24))

    if (Hsnow > Hice):
        Hice = Hice
        dHice = 0.
    else:
        Hice = Hice + dHice 

    if (Hice > Hmax):
        Hice = Hmax 
        dHice = 0.
    elif (Hice < Hmin):
        Hice = Hmin 
        dHice = 0.
        
    # Obstacle height 
    Ho = Hice - Hsnow 
    if (Ho < 0.01): Ho = 0.01
    if (Ho > Hmax): Ho = Hmax

    # Obstacle frontal area ratio
    Lo = (f/100)*Ho 

    # drag coefficient for form drag, adapted from Garbrecht at al. (2002)
    if (Ho <= 2.5527):
        Cr = 0.5 * (0.185 + (0.147*Ho))
    else:
        Cr = 0.11 * np.log(Ho/0.2)

    #IF (Hsnow.ge.0.01) Cr = Cr / 2 # drag reduction when snow is present, not used

    # Displacement height
    d = Ho * (1 - (1 - np.exp(-(c1*Lo)**0.5)) / (c1*Lo)**0.5)

    # Drag coefficient for skin friction at z = Ho
    Cs = (Cs10**(-0.5) - (1/karman) * (np.log((10-d) / (Ho-d))-PsiH))**(-2)

    # Model for u/ustar
    gamma = (Cs + Cr*Lo)**(-0.5)
    
    # Roughness length for momentum
    z0m = (Ho-d) * np.exp(-karman * gamma) * np.exp(PsiH)
    
    return z0m, Hice, dHice