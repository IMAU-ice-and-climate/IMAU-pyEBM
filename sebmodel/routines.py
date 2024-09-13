#===============================================================================
#
# File with general usefull routines and functions
#
#===============================================================================
from numba import njit
import numpy as np

from globals import *
from info import *

@njit
def julday(year,month,day):
    #===============================================================================
    # Function which calculates the Julian day (or day number) based on input 
    # year, month and day
    #===============================================================================
    dmth =  [31,28,31,30,31,30,31,31,30,31,30,31]

    # Determine whether it is a leap year
    lyr = 0 
    if (int(year/4) == year/4): lyr = 1
    if (int(year/100) == year/100): lyr = 0
    if (int(year/400) == year/400): lyr = 1

    if (lyr == 0): dmth[1] = 28
    if (lyr == 1): dmth[1] = 29

    jday = day
    for k in range(1,month-1):
        jday = jday + dmth[k]
    return jday

@njit
def spechum(t,rh,p):
	#===============================================================================
	# Function that calculates the specific humidity based on given 
	# temperature relative humidity and air pressure
	#===============================================================================

	# saturation vapor pressure
	# k=0 with respect to water, k=1 with respect to ice
	if (t < Tkel): k=1
	if (t >= Tkel): k=0
	esat = es(t,k)
	# specific humidity
	q = (eps*esat*rh) / (p-(1-eps)*esat)

	return q

@njit
def relhum(t,q,p):
	#===============================================================================
	# Function that calculates the specific humidity based on given 
	# temperature relative humidity and air pressure
	#===============================================================================
	
    # saturation vapor pressure
	# k=0 with respect to water, k=1 with respect to ice
	if (t < Tkel): k = 1
	if (t >= Tkel): k = 0
	esat = es(t,k)
	# relative humidity
	rh = (q*(p-(1.-eps)*esat)) / (eps*esat)

	return rh

@njit
def es(t,k):
	#===============================================================================
	# Function that calculates the vapor pressure of water over ice or 
	# a water surface given temperature
	#===============================================================================
 
    fact1 = 1./rv
    fact2a = (lv + beta * Tkel)
    fact2b = (ls + beta * Tkel)
    fact3 = ( (1/Tkel) - (1/t) )
    fact4 = beta * np.log(t/Tkel)

    if (k == 0):  # with respect to water (lv)
        es = es0 * np.exp(fact1 * ((fact2a*fact3)-fact4) )
    elif (k ==1):	# with respect to ice (ls)
        es = es0 * np.exp(fact1 * ((fact2b*fact3)-fact4) )
    return es

@njit
def conduc(dens,temp):
    # !===============================================================================
	# ! Function that calculates the effective conductivity of snow/ice
	# ! Density in g/kg=kg/m3
	# ! most of these are defined between 100-700 kg/m3 density
	# !===============================================================================
	# USE INPUT_EBM , ONLY : tpcond , densice
	# !USE CONSTANTS_EBM , ONLY : densice
	# IMPLICIT NONE
	# REAL :: conduc, dens, temp
	# REAL :: Kice

	if (tpcond == 1): #(4) Von Dussens equation T range? , 
		conduc = (0.21e-01+0.42e-03*dens+0.22e-08*dens**3)
	elif (tpcond == 2):		#(5) Sturm equation 156 <= rho <= 600
		conduc = (0.138-1.01e-3*dens+3.233e-6*dens**2)
	elif (tpcond == 3):	#(1) Douville equation
		conduc = 2.2*((dens/densice)**1.88)
	elif (tpcond == 4):	#(2) Jansson
		conduc = 0.02093+0.7953e-3*dens+2.512e-12*dens**4
	elif (tpcond == 5):	#(3) Anderson
		conduc = 0.021+2.5*(dens*0.001)**2
	elif (tpcond == 6):	#(6) Ostin&Anderson (6, lowest)
		conduc = -0.00871 + 0.439e-3 * dens + 1.05e-6 * dens*dens
	elif (tpcond == 7):	#(7) Arthern 1998
		conduc = 2.1*((dens/densice)**2.0)
	elif (tpcond == 8):
		Kice = 7.8042*np.exp(-0.0046507*temp)	# data from internet, fit in kaleidagraph
		conduc = ((2.*dens)/(3.*(densice-dens)))*Kice	# Cox et al., 2014 TCD
		if (dens >= densice-0.001): conduc = Kice
	elif (tpcond == 9):		# (9) tuned value
		conduc = 50.0
	elif (tpcond==10):	#!(10) Calonne et al., 2019 GRL (Eq.5 for T=Tref=-3C)
		theta =  1 / (1+np.exp(-2*0.02*(dens-450.)))
		kreff = 2.107 + 0.003618*(dens-917.)
		krefs = 0.024-1.23e-4*dens + 2.5e-6*dens**2
		conduc = (1 - theta)*krefs + theta*kreff
	# write(*,*) conduc,temp
	
	return conduc

@njit
def densprec(ws,Tskin):
# !===============================================================================
# ! Function that calculates the density of freshly fallen snow as function of
# ! wind speed and surface temperature output in g/kg=kg/m3
# ! gives value between 260 and 350
# ! based on Lenaerts et al. 2012, JGR
# !===============================================================================
# USE INPUT_EBM , ONLY : rhosnprec
# USE CONSTANTS_EBM , ONLY : Tkel
# IMPLICIT NONE
# REAL :: ws, Tskin
# REAL :: wslocal, tslocal
# REAL :: densprec
	wslocal = ws
	if (wslocal > 10.): wslocal = 10.
	tslocal = Tskin
	if (tslocal > Tkel): tslocal = Tkel
	densprec = 97.49 + 4.49 * wslocal + 0.769 * tslocal
	if (densprec < rhosnprec): densprec = rhosnprec

	return densprec

@njit
def irreducible(dens,nl):
    # !===============================================================================
    # ! Routine that calculates the potential for irreducible water content of a snow layer
    # ! based on Schneider and Jansson, 2004 J. Glaciol., 50(168), 25-34 or
    # ! Coleou and Lesaffre, 1998
    # !===============================================================================
    # USE SNOW_EBM , ONLY : nl, irrwater, dens
    # USE INPUT_EBM , ONLY : tpirre , cirre , densice
    # USE CONSTANTS_EBM , ONLY : denswater !, densice 

    # IMPLICIT NONE
    # !local
    # INTEGER :: il
    # REAL :: porosity
    # REAL :: percirrwater		! irreducible water content in % of mass according to schneider (mass of water devided by sum of masses of water and snow)
    # REAL :: dencap				! density of capillary water (kg/m3)
    # REAL :: denpor				! density of water when all pores are filled completely (kg/m3)
    # REAL :: dencol				! density of water when maximum amount according to Coleou is filled (kg/m3)
    irrwater = np.zeros(nl)
    
    if (tpirre > 1):
        for il in range(0,nl-1):
            if (dens[il] > (densice-0.001*densice)):
                irrwater[il] = 0.
            else:
                porosity = (densice - dens[il]) / densice
                if (tpirre == 2):
                    percirrwater = ((0.057*porosity)/(1.-porosity)) + 0.017
                elif (tpirre == 3):
                    percirrwater = 0.0143*np.exp(3.3*porosity)
                dencol = percirrwater/(1. - percirrwater) * dens[il]
                denpor = porosity*denswater
                dencap = dencol
                if (dencap > denpor): dencap = denpor
                irrwater[il] = dencap/(porosity*denswater)
    else:
        for il in range(0,nl-1):
            irrwater[il] = cirre
    
    return irrwater

@njit
def tempprec(t,q,p):
	# !===============================================================================
	# FUNCTION tempprec(t,q,p)
	# !===============================================================================
	# ! Function that calculates the temperature of freshly fallen snow or rain as function of
	# ! near surface air temperature, specific humidity and air pressure based on the psychrometric energy balance
	# ! from Harder P, Pomeroy J (2013) Hydrol. Process https://doi.org/10.1002/hyp.9799
	# !===============================================================================
	# USE CONSTANTS_EBM , ONLY : Tkel, rd, ls, lv
	# IMPLICIT NONE
	# REAL,INTENT(IN) :: t, q, p
	# REAL :: D, lambdat, lambdaf, esat
	# REAL :: densair, ll, qsat
	# REAL :: tempprec
	# REAL :: spechum
	if(t>Tkel):
		ll = lv # Evaporation
	else:
		ll = ls # Sublimation 

	densair = p / (rd*t)
	D = 2.06e-5 * (t/Tkel)**(1.75)
	lambdat = 6.3e-5*t+6.73e-3

	qsat = spechum(t,1.0,p)
	tempprec = t + (D/lambdat)*ll*(densair*(q-qsat))

	return tempprec
'''
def settoerror(icount):
	# !===============================================================================
	# USE INPUT_EBM , ONLY : dsnow
	# USE GLOBALS_EBM , ONLY : errorflag , buf , sbuf , t0, q0, t2m, q2m, ws10m, cloud , &
	# &                        Snet , albedo , Lnet, SH, LE, source, restsource, GH , &
	# &                        Ch,Cq,ustar,thstar,qstar,psim,psih,psiq,psim0, psih0, &
	# &                        psiq0, z0h,z0q,z0m,zt,zm , paramerror
	# USE SNOW_EBM , ONLY : hsnow, hsnowmod , hmass, precipsum, &
	# &                     dsnowacc , icemeltout, runoff, melt , subl , slushdepth , surfwater , &
	# &                     corrsnow , hsnowmod_l
	# USE RADPEN_EBM , ONLY : sumdivs
	# USE CONSTANTS_EBM , ONLY : errorval , colmax 

	# IMPLICIT NONE

	# INTEGER :: icount , iv , ip
	# INTEGER :: idT,idS,idWS,idL,idall
	# REAL :: valerror
	# REAL :: dum, dum2

	idT = 0 	! temperature = 
	idS = 0 	! short wave radiation
	idWS = 0	! wind speed
	idL = 0		! long wave radiation
	idall = 0	! no data at all

	valerror = errorval - 90.

	ip = 6	!temperature
	dum = INT(paramerror(icount) /(10**(ip-1))) - 10*INT(paramerror(icount) /(10**(ip)))
	IF (dum == 9) THEN
	idT = 1
	ENDIF
	ip = 2	!Sin
	dum = INT(paramerror(icount) /(10**(ip-1))) - 10*INT(paramerror(icount) /(10**(ip)))
	ip = 3	!Sout
	dum2 = INT(paramerror(icount) /(10**(ip-1))) - 10*INT(paramerror(icount) /(10**(ip)))
	IF (NINT(dum) == 9 .or. NINT(dum2) == 9) THEN
	idS = 1
	ENDIF
	ip = 1	!wind speed
	dum = INT(paramerror(icount) /(10**(ip-1))) - 10*INT(paramerror(icount) /(10**(ip)))
	IF (dum == 9) THEN
	idWS = 1
	ENDIF
	ip = 4	!Lin
	dum = INT(paramerror(icount) /(10**(ip-1))) - 10*INT(paramerror(icount) /(10**(ip)))
	IF (dum == 9) THEN
	idL = 1
	ENDIF

	DO iv = 4,colmax
	sbuf(iv) = buf(icount,iv)
	ENDDO

	IF (idT == 1) THEN
	sbuf(4) = valerror		! T
	sbuf(11) = valerror/1000.	! q
	ENDIF

	IF (idS == 1) THEN
	sbuf(7) = valerror		! Sin	
	sbuf(8) = valerror		! Sout
	sbuf(9) = valerror		! albedo
	ENDIF

	IF (idWS == 1) THEN
	sbuf(6) = valerror		! ws
	ENDIF

	IF (idL == 1) THEN
	sbuf(11) = valerror		! lin
	ENDIF

	t2m = valerror
	t0 = valerror
	q2m = valerror/1000.
	q0 = valerror/1000.
	ws10m = valerror
	cloud = valerror
	albedo = valerror

	sumdivs = valerror	 
	Snet = valerror	
	restsource = valerror
	source = valerror
	Lnet = valerror
	SH = valerror
	LE = valerror
	GH = valerror

	runoff = valerror
	melt = valerror
	subl = valerror
	slushdepth = valerror
	surfwater = valerror

	Ch = valerror
	Cq = valerror
	ustar = valerror
	thstar = valerror
	qstar = valerror/1000.
	psim0 = valerror
	psim = valerror
	psih0 = valerror
	psih = valerror
	psiq0 = valerror
	psiq = valerror
	z0m = valerror
	z0h = valerror
	z0q = valerror
	zt = valerror
	zm = valerror

	! = start snow, start sbuf(16), current sbuf(16), 16 = serie mean
	hsnow = sbuf(16)		!current value serie
	corrsnow = hsnow - hsnowmod_l

	hsnowmod = valerror
	hmass = valerror
	precipsum = valerror
	dsnowacc = valerror
	icemeltout = valerror

END SUBROUTINE SETTOERROR


!===============================================================================
SUBROUTINE SET_RANDOM_SEED
!===============================================================================
!Subroutine that seeds the random number generator depending on system time
!===============================================================================
	IMPLICIT NONE
	
	INTEGER				:: i, clock
	INTEGER, PARAMETER	:: n=33
	INTEGER				:: seed(n)
	
	CALL SYSTEM_CLOCK(COUNT=clock)
	seed = clock + 37 * (/ (i - 1, i = 1, n) /)
	CALL RANDOM_SEED(PUT = seed)
END SUBROUTINE SET_RANDOM_SEED
!===============================================================================
SUBROUTINE RANDHUM(q,temp,press,tempold,pressold)
!===============================================================================
!Subroutine that randomly disturbs the relative humidity measurement
!===============================================================================
	USE CONSTANTS_EBM, ONLY: Tkel
	USE INPUT_EBM, ONLY: chstation
	
	IMPLICIT NONE
	
	REAL, INTENT(IN)	:: temp, press, tempold, pressold
	REAL, INTENT(INOUT)	:: q
	!Already defined functions
	REAL	:: relhum, spechum
	!Local
	REAL	:: nudge, rh
	
	rh = relhum((tempold+Tkel),(q/1000.),(pressold*100.))
	
	CALL RANDOM_NUMBER(nudge)
	nudge = (2.*nudge)-1.
	
	IF(chstation.eq."ant_neuma") THEN
	 rh = rh*(1.+0.05*nudge) !5% error
	ELSE
	 IF(rh .gt. 0.9) THEN
	  rh = rh*(1.+0.03*nudge) !If rel hum >90%: 3% error
	 ELSE
	  rh = rh*(1.+0.02*nudge) !Else 2% error
	 ENDIF
	ENDIF
	IF(rh .gt. 1.) rh = 1.
	
	q = (spechum((temp+Tkel),rh,(press*100.)))*1000.
	
END SUBROUTINE RANDHUM
!===============================================================================
FUNCTION solvetemp(temp_old,energy,mass)
!===============================================================================
! Function that finds new layer temperature from its energy and mass
!===============================================================================
    USE CONSTANTS_EBM, ONLY: Tkel
    
    IMPLICIT NONE
    
    !Input
    REAL, INTENT(IN)    :: temp_old, energy, mass
    !Output
    REAL    :: solvetemp
    
    !Local
    REAL,PARAMETER  :: a=7.122, b=152.5
    REAL            :: c, discriminant
    REAL            :: sol1, sol2
    
    c = (energy/mass) - (152.5+7.122*Tkel)*Tkel
    
    discriminant = b**2. - (4.*a*c)
	IF(discriminant.lt.0) THEN
	 WRITE(*,*) "Can't solve second order polynomial!"
	 STOP 60
	END IF
	
	!Usually it's solution with positive sign
    sol1 = ((-b)+SQRT(discriminant))/(2.*a)
    IF(ABS(sol1-temp_old) .lt. 10.) THEN
     solvetemp = sol1
    ELSE
     sol2 = ((-b)-SQRT(discriminant))/(2.*a)
    
     !Now choose solution that is closest to previous value
     IF(ABS(sol1 - temp_old) < ABS(sol2 - temp_old)) THEN
      solvetemp = sol1
     ELSE
      solvetemp = sol2
     ENDIF
    ENDIF
END FUNCTION
!===============================================================================
FUNCTION solvequadpos(a,b,c)
!===============================================================================
! Function that solves a simple quadratic function using the abc-formula
!===============================================================================
	IMPLICIT NONE
	
	!Input
	REAL,INTENT(IN) 			:: a, b, c
	!Output
	REAL	:: solvequadpos
	
	!Local
	REAL	:: discriminant
	
	discriminant = b**2.-(4.*a*c)
	IF(discriminant.lt.0) THEN
	 WRITE(*,*) "Can't solve second order polynomial!"
	 STOP 60
	END IF
	
	solvequadpos = ((-b)+SQRT(discriminant))/(2.*a)
	
END FUNCTION solvequadpos
!===============================================================================
FUNCTION solvequadneg(a,b,c)
!===============================================================================
! Function that solves a simple quadratic function using the abc-formula
!===============================================================================
	IMPLICIT NONE
	
	!Input
	REAL,INTENT(IN) 			:: a, b, c
	!Output
	REAL	:: solvequadneg
	
	!Local
	REAL	:: discriminant
	
	discriminant = b**2.-(4.*a*c)
	IF(discriminant.lt.0) THEN
	 WRITE(*,*) "Can't solve second order polynomial!"
	 STOP 61
	END IF
	
	solvequadneg = ((-b)-SQRT(discriminant))/(2*a)
	
END FUNCTION solvequadneg
'''