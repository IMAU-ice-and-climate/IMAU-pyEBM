#!/usr/bin/env python3
"""
===============================================================================
All routines related to the initialisation of arrays and parameters.
===============================================================================
"""
""" Standard library imports """
import pandas as pd
import numpy as np1 
""" Third-party import"""

"""Local imports"""
from info import * 
# from input import * 
from globals import *

def setconstants():
    """
    Decalare constants
    """
    
def settozero():
    """
    Declare arrays and set them to zero
    """
    
    ibuf = np.zeros([rowmax,colmax])
    zenith = np.zeros([rowmax,1])

'''
!===============================================================================
SUBROUTINE SETTOZERO
!===============================================================================
USE GLOBALS_EBM, ONLY:	ilasttot , ibuf, buf, sbuf, dhour, chdtg , & 
&                       zenith, dbuf , dalbedo , dalbedo_obs , &
&						mbuf , mcount , ssbuf , scount ,clbuf, clcount, clmth , ybuf , ycount , perbuf , percount , &
&						t0, q0, Snet, Lnet, SH, LE, GH, source, restsource, Ch, Cq , snowdays , &
&						buflast , errorflag , merror , serror , &
&						yerror , clerror , derror , awsid , dcount , &
&						dt0buf , dt0sum, dt0stdsum , dt0sumcount, dt02sum , &
&                       t0sum, t0obssum, t0t0obssum, t02sum, t0obs2sum, &
&                       pearsonR , daccsum, dacc2sum , dmelsum, dmel2sum , &
&						dsmbsum, dsmbsumcount, dsmb2sum , icemeltcorrmod , icemeltcorrobs ,  &
&						cumicemeltobs , snowcorrobs , accobs , paramerror, Hice
USE SNOW_EBM, ONLY:		dz, z, temp, dens, mass, cpice, rhocp, kice, &
&                       water, energy, icemeltmdt, &
&						irrwater, ice, lid, vink, freshfrac, refrfrac, &
&						hsnowmod , dsnowr , dsnowh , dsnowacc ,  hmass ,&
&						runoff , surfmelt, melt , subl , precipsum , icemelt , &
&                       surfwater , slushdepth, &
&						sumrunoff, sumsurfmelt, summelt , icemeltmin , &
&						startmelt , endmelt, stmeltout , ndmeltout , &
&						winterbal, summerbal, annualbal, winterbalmod, summerbalmod, annualbalmod , &
&                       wintermelt, summermelt, annualmelt , winterbalmsn, summerbalmsn, annualbalmsn,&
&						mbout , mpdd , mtsumpdd , mnrmeltdys, cpdd , ctsumpdd , cnrmeltdys , &
&						spdd , stsumpdd , snrmeltdys, ypdd , ytsumpdd , ynrmeltdys, &
&                       dtdz, grainsize
USE CONSTANTS_EBM, ONLY: rowmax , colmax , nlmax , mmax
USE INPUT_EBM, ONLY:	lcomment , ibyear, ilyear, lz0m, Hmax
USE RADPEN_EBM, ONLY:	dsdz, sumdivs , zrad , nlradmax
USE PKM_EBM, ONLY:		drdry

IMPLICIT NONE

ALLOCATE (dhour(rowmax))
ALLOCATE (ibuf(rowmax,colmax))
ALLOCATE (buf(rowmax,colmax))
ALLOCATE (sbuf(colmax))
ALLOCATE (dbuf(mmax))
ALLOCATE (dcount(mmax))
ALLOCATE (mbuf(mmax))
ALLOCATE (mcount(mmax))
ALLOCATE (ssbuf(mmax))
ALLOCATE (scount(mmax))
ALLOCATE (clbuf(12,mmax))
ALLOCATE (clcount(12,mmax))
ALLOCATE (ybuf(mmax))
ALLOCATE (ycount(mmax))
ALLOCATE (perbuf(mmax))
ALLOCATE (percount(mmax))
ALLOCATE (chdtg(rowmax))
ALLOCATE (paramerror(rowmax))
ALLOCATE (zenith(rowmax))
ALLOCATE (errorflag(rowmax))
ALLOCATE (buflast(colmax))
ALLOCATE (awsid(rowmax+1))
ALLOCATE (dt0buf((ilyear - ibyear + 1)*rowmax))

ALLOCATE (dz(nlmax))
ALLOCATE (z(nlmax))
ALLOCATE (temp(nlmax))
ALLOCATE (dtdz(nlmax))
ALLOCATE (dens(nlmax))
ALLOCATE (mass(nlmax))
ALLOCATE (cpice(nlmax))
ALLOCATE (rhocp(nlmax))
ALLOCATE (kice(nlmax))
ALLOCATE (dsdz(nlmax))
ALLOCATE (water(nlmax))
ALLOCATE (irrwater(nlmax))
ALLOCATE (energy(nlmax))
ALLOCATE (ice(nlmax))
ALLOCATE (grainsize(nlmax))
ALLOCATE (freshfrac(nlmax))
ALLOCATE (refrfrac(nlmax))
ALLOCATE (lid(nlmax))
ALLOCATE (zrad(nlradmax))

ALLOCATE (drdry(nlmax))

IF (lcomment == 1) WRITE(*,'(/,A)') "Initialise all arrays to zero"
ilasttot = 0

dhour = 0
ibuf = 0.0
buf = 0.0
sbuf = 0.0
dbuf = 0.0
mbuf = 0.0
ssbuf = 0.0
 clbuf = 0.0
ybuf = 0.0
perbuf = 0.0
dcount = 0
mcount = 0
scount = 0
 clcount = 0
 clmth = 0
ycount = 0
percount = 0
zenith = 0.
dalbedo = 0.
dalbedo_obs = 0.
errorflag = 0.
derror = 0.
merror = 0.
serror = 0.
yerror = 0.
 clerror = 0.
buflast = -999.9
awsid = 0.0
paramerror = 0

lid = 0		!ice
dz = 0.0
z = 0.0
temp = 0.0
dtdz = 0.0
dens = 0.0
mass = 0.0
 cpice = 0.0
rhocp = 0.0
kice = 0.0
dsdz = 0.0
water = 0.0
irrwater = 0.0
energy = 0.0
ice = 0.0
freshfrac = 0.0
refrfrac = 0.0

zrad = 0.

grainsize = 0.0
drdry = 0.0

t0 = 0.0
q0 = 0.0
Snet = 0.0
sumdivs = 0.0
Lnet = 0.0
SH = 0.0
LE = 0.0
GH = 0.0
source = 0.0
restsource = 0.0
 Ch = 0.0
 Cq = 0.0

vink = 0
hsnowmod = 0.0
dsnowr= 0.0

runoff = 0.
surfmelt = 0.
melt = 0.
subl = 0.
surfwater = 0.
slushdepth = 0.

sumrunoff = 0.
sumsurfmelt = 0.
summelt = 0.

precipsum = 0.
icemelt = 0.
icemeltmdt = 0.
icemeltmin = 10.
dsnowh = 0.
dsnowacc = 0.
hmass = 0.

snowdays = 0.

dt0buf = -999.0
dt0sum = 0.
dt0sumcount = 0
dt02sum = 0.
dt0stdsum = 0.
t0sum = 0.
t0obssum = 0.
t0t0obssum = 0.
t02sum = 0.
t0obs2sum = 0.
pearsonR = 0.
daccsum = 0.
dacc2sum = 0.
dmelsum = 0.
dmel2sum = 0.
dsmbsum = 0.
dsmbsumcount = 0
dsmb2sum = 0.
icemeltcorrmod = 0.
icemeltcorrobs = 0.
snowcorrobs = 0.
 cumicemeltobs = 0.
accobs = 0.

startmelt = 999. 
endmelt = -999.
stmeltout = -999.
ndmeltout = -999.

winterbal = 0.	!-999.
summerbal = 0.	!-999.
annualbal = 0.	!-999.
winterbalmsn = 0.	!-999.
summerbalmsn = 0.	!-999.
annualbalmsn = 0.	!-999.
winterbalmod = 0.	!-999.
summerbalmod = 0.	!-999.
annualbalmod = 0.	!-999.
wintermelt = 0.	!-999.
summermelt = 0.	!-999.
annualmelt = 0.	!-999.
mbout = -999.

mpdd = 0.
mtsumpdd = 0.
mnrmeltdys = 0.
spdd = 0.
stsumpdd = 0.
snrmeltdys = 0.
 cpdd = 0.
 ctsumpdd = 0.
 cnrmeltdys = 0.
ypdd = 0.
ytsumpdd = 0.
ynrmeltdys = 0.

if (lz0m == 4) Hice = Hmax

END SUBROUTINE SETTOZERO

!===============================================================================
SUBROUTINE FREEARRAYS
!===============================================================================

USE GLOBALS_EBM, ONLY:	ibuf, buf, sbuf, dhour, chdtg , zenith , buflast, dbuf, &
&						mbuf, ssbuf , clbuf , ybuf , perbuf , errorflag , awsid , dcount , mcount, &
&						scount , ycount , percount , paramerror , dt0buf
USE SNOW_EBM, ONLY:		dz, z, temp, dens, mass, cpice, rhocp, kice, water, energy, &
&						irrwater, ice , lid, dtdz, freshfrac, refrfrac, grainsize
USE RADPEN_EBM, ONLY:	dsdz , zrad
USE PKM_EBM, ONLY:		drdry

IMPLICIT NONE

DEALLOCATE (ibuf)
DEALLOCATE (buf)
DEALLOCATE (sbuf)
DEALLOCATE (dhour)
DEALLOCATE (chdtg)
DEALLOCATE (paramerror)
DEALLOCATE (zenith)
DEALLOCATE (errorflag)
DEALLOCATE (buflast)
DEALLOCATE (dbuf)
DEALLOCATE (dcount)
DEALLOCATE (mbuf)
DEALLOCATE (mcount)
DEALLOCATE (ssbuf)
DEALLOCATE (scount)
DEALLOCATE (clbuf)
DEALLOCATE (ybuf)
DEALLOCATE (ycount)
DEALLOCATE (perbuf)
DEALLOCATE (percount)
DEALLOCATE (awsid)
DEALLOCATE (dt0buf)

DEALLOCATE (dz)
DEALLOCATE (z)
DEALLOCATE (temp)
DEALLOCATE (dtdz)
DEALLOCATE (dens)
DEALLOCATE (mass)
DEALLOCATE (cpice)
DEALLOCATE (rhocp)
DEALLOCATE (kice)
DEALLOCATE (dsdz)
DEALLOCATE (water)
DEALLOCATE (irrwater)
DEALLOCATE (energy)
DEALLOCATE (ice)
DEALLOCATE (freshfrac)
DEALLOCATE (refrfrac)
DEALLOCATE (lid)

DEALLOCATE (grainsize)
DEALLOCATE (drdry)

DEALLOCATE (zrad)

END SUBROUTINE FREEARRAYS

!===============================================================================
'''