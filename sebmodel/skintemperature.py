from numba import njit

from globals import *
from info import *
from fluxes import turbhf, energybalance
from routines import spechum

@njit
def surft(t1,t2,dz1,dz2,z1,z2):
    #===============================================================================
    # T0 calculation by means of linear extrapolation of two upper layers
    # temperature
    #===============================================================================
    # USE CONSTANTS_EBM , ONLY : tinterv, taccur, Tkel

    # REAL :: surft
    # REAL,INTENT(IN) :: t1, t2, dz1, dz2, z1, z2
    # REAL :: tgrad,tsurf

    tgrad = (t2-t1)/(z2-z1)
    tsurf = t1 - tgrad*z1

    # Surface temperature of melting snow cannot exceed 0.0 C or 273.16 K
    if (tsurf > Tkel): tsurf = Tkel

    return tsurf

@njit
def tskin(t0,Snet,sumdivs,kice,temp,dz,sbuf):
#===============================================================================
## 1=year, 2=day, 3=hour, 4=T, 5=p, 6=WS, 7=Sin, 8=Sout, 9=albedo, 10=running mean albedo, 
# 11=qi, 12=Lin, 13=Lout, 14=z0m, 15=zt, 16=zm , 17=serie(running mean), 18=precip,  
#===============================================================================
    # USE GLOBALS_EBM , ONLY : t0 , q0 , sbuf, source, densair 
    # USE CONSTANTS_EBM , ONLY : tinterv, taccur, Tkel, rd
    # USE INPUT_EBM , ONLY : lcomment
    # USE FILES , ONLY : uo1

  
    # INTEGER :: iter,itermax
    # REAL :: t01,t02,dt0,t0old
    # REAL :: t0old1,dt01
    # REAL :: bisection, energybalance, spechum, falsepos
    
    zm = sbuf[14]
    zt = sbuf[13]
    z0m = sbuf[12]

    ustar, SH, LE, Ch, Cq, Chn, Cqn, t2m, q0, q2m, ws10m, densair = turbhf(0,t0,sbuf)

    t01 = t0 - tinterv
    t02 = t0 + tinterv
    dt0 = 2*taccur
    dt01 = 2*taccur
    t0old = t01

    iter = 0
    itermax = 40

    source = 0.

    while ((dt0 > taccur) & (dt01 > 0)):
        iter = iter + 1

        t0old1 = t0old
        t0old = t0
        
        if (energybalance(Tkel,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf)[0] >= 0.0):
            t0 = Tkel
        else:
            if (t02 > 273.2): t02 = 273.2 # Resulting skin temperature will certainly not exceed this value
            #  t0 = bisection(t01,t02,taccur,yr,ii)
            t0 = falsepos(t01,t02,taccur,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf)
        
        if (t0 >= Tkel): t0 = Tkel
        q0 = spechum(t0, 1.0, sbuf[4])
        
        ustar, SH, LE, Ch, Cq, Chn, Cqn, t2m, q0, q2m, ws10m, densair = turbhf(iter,t0,sbuf)

        # source = energybalance(t0)
        if (t0 >= Tkel): source, _  = energybalance(Tkel,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf)
        if (source < 0.): source = 0.				
                        
        dt0 = abs(t0 - t0old)			# TEST TEST TEST TEST
        dt01 = abs(t0 - t0old1)

        if (iter >= itermax):		# no solution found
            if (lcomment == 1):
                print('TSKIN more than itermax necessary')
                # WRITE(*,'(/,A,i3,A,/,I5,6f16.8,/)') 'TSKIN more than',itermax,' iterations necessary',iter,taccur,dt0,t0,t0old,source
            #  STOP
            # WRITE(uo1,'(/,A,i3,A,/,I5,5f16.8,/)') 'TSKIN more than',itermax,' iterations necessary',iter,taccur,dt0,t0,t0old,source
            t0 = 0.5*(t0+t0old)
            dt0 = 0.

    if ((dt0 > taccur) & (dt01 == 0.) & (iter > itermax*0.5)): # no solution found
    #  IF (lcomment == 1) THEN
    #    WRITE(*,'(/,A,/,I5,5f16.8,/)') 'TSKIN no solution found, varies between ',&
    #&       iter,taccur,dt0,t0,t0old
    #  ENDIF
        # WRITE(uo1,'(/,A,/,I5,5f16.8,/)') 'TSKIN no solution found, varies between ',&&       iter,taccur,dt0,t0,t0old
        t0 = 0.5*(t0+t0old)
        
    return t0, ustar, SH, LE, Ch, Cq, Chn, Cqn, t2m, q0, q2m, ws10m, densair

@njit
def falsepos(a,b,acc,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf):
# ===============================================================================
#  This function is similar to the one above but it uses the "false position method"
# ===============================================================================
# 	IMPLICIT NONE
# 	REAL				:: a, b, c, acc
# 	REAL				:: falsepos, energybalance
# 	REAL				:: fa, fb, fc, dx, rtb
# 	INTEGER				:: j
# 	INTEGER, PARAMETER	:: jmax = 40
# 	REAL, PARAMETER		:: ebacc = 0.001 !Accuracy of energy balance rest term
# 	!	b>a, so f(b)<f(a) (because f(x) is monotonically descending)
    fa, _ = energybalance(a,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf)	
    fb, _ = energybalance(b,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf)	
    dx = abs(a-b)
    jmax = 40 
    ebacc = 0.001
    while ((fa*fb) >= 0.0): # If we don't capture the root
        if(fb >= 0): # If we're left of the root
            a = b
            b = b + dx
        else: # Right
            b = a
            a = a - dx
        fa, _  = energybalance(a,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf)
        fb, _  = energybalance(b,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf)

    j = 1
    c = b - fb*((b-a)/(fb-fa))
    fc, _  =  energybalance(c,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf)

    while((j < jmax) & (dx > acc) & (abs(fc) > ebacc)):
        if (fa*fc >= 0.0):
            a = c
            fa = fc
        else:
            b = c
            fb = fc
        dx = abs(a-b)

        j = j + 1
        c = b - fb*((b-a)/(fb-fa))
        fc, _  = energybalance(c,Snet,sumdivs,Ch,Cq,densair,kice,temp,dz,sbuf)

    if (j == jmax):
        print('Check iterations in routine falsepos ')
        # WRITE(*,'(A,I2,A)') "Maximum (>", jmax, ") number of false position iterations!"
        # STOP 70
    return c
