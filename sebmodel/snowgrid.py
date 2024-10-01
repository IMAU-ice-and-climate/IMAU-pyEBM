
import numpy as np
from numba import njit

from globals import *
from info import *
from routines import conduc, irreducible
from functions import snowdens, snowtemp
# from variables import *
#===============================================================================
#
# Routines related to the sub-surface snow model excluding the calculation
# of the surface temperature
#
#!===============================================================================
@njit
def initsnow(z,dz,lid,nl):
    # !===============================================================================
    # ! Initialises the sub-surface profiles of density, temperature and water content
    # !===============================================================================
    # USE SNOW_EBM , ONLY : nl, z, dz, temp, dens, mass, rhocp, cpice, kice, lid, energy
    # USE INPUT_EBM , ONLY : chstation, lcomment ,dfirn , penetration, tcalc, dzrad
    # USE CONSTANTS_EBM , ONLY : Tkel 
    # USE RADPEN_EBM , ONLY : zrad , nlradmax
    # USE FILES , ONLY : uo1

    # IMPLICIT NONE
    # !Local
    # INTEGER :: il
    # REAL :: conduc
    # REAL :: snowtemp, snowdens
    # REAL :: snowtemp_STGL, snowdens_STGL
    # !Assuming dry snow

    # IF (lcomment == 1) WRITE(*,'(/,A)') 'Initialise the subsurface, temp and dens'
    # WRITE(uo1,'(/,A)') 'Initialise the subsurface, temp and dens'
    nlradmax = int(zradmax/dzrad)
    
    dens = np.zeros(nlmax)
    temp = np.zeros(nlmax)
    mass = np.zeros(nlmax)
    cpice = np.zeros(nlmax)
    rhocp = np.zeros(nlmax)
    kice = np.zeros(nlmax)
    energy = np.zeros(nlmax)
    zrad = np.zeros(nlradmax)
    
    for il in range(0,nl):
    # !!!determine general profile based on either measurements or theory!
        dens[il] = snowdens(z[il],lid[il])
        temp[il] = snowtemp(z[il])
        mass[il] = dens[il]*dz[il]
        cpice[il] = 152.5+7.122*temp[il]			# from Patterson1994
        rhocp[il] = dens[il]*cpice[il]
        kice[il] = conduc(dens[il],temp[il])		# effective conductivity
        if (il == nl-1):
            energy[il] = energy[il-1]
        else:
            energy[il] = dens[il]*dz[il]*cpice[il]*(Tkel-temp[il])
        #  energy(il) = dens(il)*dz(il)*((152.5+7.122*Tkel)*Tkel - cpice(il)*temp(il))
    #kice(1) = 50.0
    irrwater = irreducible(dens,nl)

    if (penetration == 1):
        for il in range(0,nlradmax):
            zrad[il] = ((il+1)-0.5)*dzrad

    # ! insert temperature profile only changes results significantly for first
    # ! two days. Put something in based on day of year and given amplitude.
    
    return dens, temp, mass, cpice, rhocp, kice, energy, zrad

@njit
def initgrid(dsnow):
# !===============================================================================
# ! Initialises the sub-surface grid
# !===============================================================================
# USE SNOW_EBM , ONLY : nl, nlinit, nlsnow, dz, z , lid , lidmax , dsnowacc , &
# &                     sumdrift , cumdrift , cumdriftobs , hmass , mass
# USE INPUT_EBM , ONLY : lcomment, dz0, dzdeep, zdeep, dsnow, dfirn
# USE CONSTANTS_EBM , ONLY : nlmax
# USE FILES , ONLY : uo1


    # No annual layers
    small = 0.0000001

    depth = 0.0
    coeff = (dzdeep - dz0)/zdeep
    id = 1		# snow layer
    nlsnow = 1
    if ((dsnow < dz0) & (dfirn < dz0)):
        id = 0			# ice layer
        nlsnow = 0
    elif (dsnow < dz0):
        id = 2			# firn layer, counted as nlsnow
    zfirn = dsnow + dfirn
    if (zfirn > zdeep): zfirn = zdeep

    depthl = np.zeros(nlmax)
    dz = np.zeros(nlmax)
    z = np.zeros(nlmax)
    lid = np.zeros(nlmax)
    mass = np.zeros(nlmax)
    
    depthl[0] = dz0 
    dz[0] = dz0							# thickness of layer
    z[0] = 0.5*dz0						# depth of layer = half way the layer, z positive downwards
    depth = depth + depthl[0]
    lid[0] = id

    il = 1
    while (depthl[il-1] < zdeep):
        depthl[il] = depthl[il-1] + coeff*depth + dz0
        dz[il] = depthl[il] - depthl[il-1]
        z[il] = depth + 0.5*dz[il]
        lid[il] = id
        if ((depthl[il] > dsnow) & (depthl[il] < zfirn) & (id > 0) & (id < 2)):
            depthl[il] = dsnow
            dz[il] = depthl[il] - depthl[il-1]
            z[il] = depth + 0.5*dz[il]
            depth = depthl[il]
            id = 2
            if (zfirn == 0.): id = 0
            if (dz[il] < dz0):
                depthl[il] = dsnow
                dz[il-1] = depthl[il] - depthl[il-1] + dz[il-1]
                z[il-1] = z[il-1] + 0.5*dz[il]
                il = il-1
                depthl[il] = dsnow
                depth = depthl[il]
        elif ((depthl[il] > zfirn) & (depthl[il] < zdeep) & (id > 0)):
            depthl[il] = zfirn
            dz[il] = depthl[il] - depthl[il-1]
            z[il] = depth + 0.5*dz[il]
            depth = depthl[il] 
            id = 0
            if (dz[il] < dz0):
                depthl[il] = zfirn
                dz[il-1] = depthl[il] - depthl[il-1] + dz[il-1]
                z[il-1] = z[il-1] + 0.5*dz[il]
                il = il - 1
                depthl[il] = zfirn
                depth = depthl[il]
        elif (depthl[il] >= zdeep - small): 
            depthl[il] = zdeep
            dz[il] = depthl[il] - depthl[il-1]
            if (dz[il] < 0.5*dz0):
                dz[il] = dz[il-1]
                depthl[il] = depth + dz[il]
            z[il] = depth + 0.5*dz[il]
            break
        else:
            depth = depthl[il]
            
        if (lid[il] > 0): nlsnow = nlsnow + 1
        il = il + 1
        if (il == nlmax):
            raise Exception("'Number of layers exceeds maximum possible layers")
        # WRITE(*,*) 'Number of layers exceeds maximum possible layers'

    # if (lid[il] > 0): nlsnow = nlsnow + 1
    nl = il
    nlinit = nl

    sumdrift = 0.
    cumdrift = 0.
    cumdriftobs = 0.
    dsnowacc = 0.
    hmass = 0.
    
    lidmax = 1
    for il in range(0,nlsnow):
        if (lid[il] == 1):	#this mb year snow layer
            dsnowacc = z[il] + 0.5*dz[il]	# in m snow
            hmass = hmass + mass[il]		#in m we
    
    if (lcomment == 1):
        print('initial number of layers is: ',nl,nlsnow,dz[nl-1],z[nl-1],lidmax)
    return z, dz, depthl, lid, dsnowacc, hmass, nl, nlinit

@njit
def resizegrid(mass, dens, lid, dz, z, nl, nlinit, lidmax, vink):
# !===============================================================================
# ! resizes the sub-surface grid after densification
# ! mass = rho*dz, dz = mass/rho
# !===============================================================================
# USE SNOW_EBM , ONLY : nl, z, dz, mass, dens, vink , nlinit , lid , hmass , &
# &                     dsnowh , dsnowacc , lidmax
# USE INPUT_EBM , ONLY : dz0, dzdeep, zdeep  , densice , lcomment
# USE CONSTANTS_EBM , ONLY : nlmax , accur 

# IMPLICIT NONE

# ! local
# REAL :: sumsnow , summass , coeff
# INTEGER :: il

    coeff = 0.75*(dzdeep - dz0)/zdeep + 2.

    sumsnow = 0.	# local snow layer in m snow
    summass = 0.	# local snow layer in m we
    dsnowacc = 0.	# snow layer in m snow
    hmass = 0.		# snow layer in m we

    dz[0] = mass[0]/dens[0]
    z[0] = 0.5*dz[0]
    if ((dz[0] < 0.5*dz0) | (dz[0] > coeff*dz0)): 
        vink = 1
        # print('vink = 1',nl,dz[0],z[0],mass[0],dens[0])
    if (lid[0] == 1):
        dsnowacc = dz[0]
        hmass = mass[0]

    if ((lid[0] > 0) & (lid[0] <= lidmax)):
        sumsnow = dz[0]
        summass = mass[0]


    if ((z[0] > 0.5*dz[0]) | (dz[0] <= 0.)):
        il = 0
        if (lcomment == 1) : 
            print('1 resizegrid',nl,il,dz[il],z[il],mass[il],dens[il])

    for il in range(1,nl):
        dz[il] = mass[il]/dens[il]
        z[il] = 0.5*(dz[il]+ dz[il-1]) + z[il-1]
        if (dz[il] < 0.5*dz0): 
            vink = 2
            # print(nl,'resizegrid',nl,il,dz[il],z[il],mass[il],dens[il])
        if (dz[il] <= 0): 
            print(nl,'resizegrid',nl,il,dz[il],z[il],mass[il],dens[il])
        # if ((dz[il] <= 0.) & (lcomment == 1)):
            # WRITE(*,'(a,2i6,4f12.5)') 'nl RESIZEGRID ',nl,il,dz(il),z(il),mass(il),dens(il)
        if ((lid[il] > 0) & (lid[il] <= lidmax)):
            sumsnow = z[il] + 0.5*dz[il]
            summass = summass + mass[il]
        if (lid[il] == 1):
            dsnowacc = sumsnow
            hmass = summass

    if (sumsnow > 0.):
        dsnowh = sumsnow
    else:
        dsnowh = 0.
        dsnowacc = 0.
        hmass = 0.

    # limit number of layers to speed up calculations
    if ((nl > 3.0*nlinit) | (nl >= nlmax-1)): vink = 3
    
    return z, dz, dsnowh, dsnowacc, hmass, vink
@njit
def redefgrid(water, ice, mass, dens, lid, dz, z, nl, nlsnow, temp, grainsize, t0, dsnowr,vink):
    # !===============================================================================
    # ! redefines the sub-surface grid
    # ! now only works in case only snow is present, no ice layer below
    # !===============================================================================

    # No annual layers or ice below the snow

    # It is not necessary to keep track of refrfrac when splitting/merging layers because it is calculated in SNOWCONTENT

    small = 0.0000001

    nl_old = nl
    nlsnow_old = nlsnow

    nl = 0
    nlsnow = 0

    temp_old = temp
    water_old = water
    dens_old = dens
    mass_old = mass
    ice_old = ice
    z_old = z
    dz_old = dz
    lid_old = lid
    grainsize_old = grainsize

    temp_tmp = 0.
    water_tmp = 0.
    dens_tmp = 0.
    mass_tmp = 0.
    ice_tmp = 0.
    z_tmp = 0.
    dz_tmp = 0.
    grainsize_tmp = 0.
    

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
    cpice = np.zeros(nlmax)
    dtdz = np.zeros(nlmax)

    coeff = 0.75*(dzdeep - dz0)/zdeep

    depth = 0.0				# depth of upper layer boundary, depthl = depth of lower boundary
    depthl = np.zeros(nlmax)
    ike = -1
    ik = -1
    il = 0
    check = 0
    while ((il <= nlmax-1) & (ike < (nl_old - 1)) & (depth < z_old[nl_old-1]+0.5*dz_old[nl_old-1]) & (nl < nlmax)):
        depthl[il] = depth*(1.+coeff) + dz0
        iks = ike + 1
        ik = iks
        ilp = 0
        while ((ik < nl_old - 1) & (depthl[il] > (z_old[ik]+0.5*dz_old[ik]+small)) & (lid_old[ik+1] == lid_old[iks])):
            ik = ik + 1
        ike = ik
        if ((ike == 0) & (depthl[il] > 0.5*dz0) & (lid_old[ik] != lid_old[ik+1]) & (z_old[ik] <= 0.5*dz0)):
            ike = 1
            dsnowr = mass_old[0]/dens_old[0] #dz_old(1)
        ilp = 0
        depthnext = depthl[il]*(1.+coeff) + dz0
        while (depthnext <= (z_old[ike]+0.5*dz_old[ike])):
            depthnext = depthnext*(1.+coeff) + dz0
            ilp = ilp + 1
        if (iks == ike) :
            check = 1 # nothing to do with the layer
            if (ilp > 0): 
                check = 3	# split layers
        else:
            check = 2 # merge layers
            if ((ike == 1) & ((z_old[ike]+0.5*dz_old[ike]) > 1.5*dz0)):
                check = 4		# first merge and then split in 2 layers
            elif (ilp > 1):
                check = 5		# first merge and then split in multiple layers
            if ((dens_old[iks] > 700.) &  (dz_old[iks] > dz0)):
                check = 1		# nothing to do with the layer
                ike = iks
        # print('REDEFGRID, il = ',il,', ike = ',ike,', iks = ',iks,', ilp = ',ilp, 'ik = ',ik, ', nl_old = ',nl_old, ', nl = ',nl, ', check = ',check)

        if ((check <= 2) | (check == 4)):		# keep layers as is or merge them
            for ik in range(iks,ike+1):
                temp[il] = temp[il] + temp_old[ik]*dz_old[ik]
                grainsize[il] = grainsize[il] + grainsize_old[ik]*dz_old[ik]
                water[il] = water[il] + water_old[ik]
                ice[il] = ice[il] + ice_old[ik]
                mass[il] = mass[il] + mass_old[ik]
                dz[il] = dz[il] + dz_old[ik]
                dens[il] = mass[il]/dz[il]

            temp[il] = temp[il]/dz[il]
            grainsize[il] = grainsize[il]/dz[il]
            dens[il] = mass[il]/dz[il]
            if (dens[il] > densice): dens[il] = densice
            z[il] = z_old[ike]+0.5*dz_old[ike]-0.5*dz[il]
            lid[il] = lid_old[ike]
            depth = z[il] + 0.5*dz[il]
            if (lid[il] > 0): nlsnow = nlsnow + 1
            nl = nl + 1
            il = il + 1
            if (check == 4):
                water[il-1] = 0.5*water[il-1]
                ice[il-1] = 0.5*ice[il-1]
                mass[il-1] = 0.5*mass[il-1]
                dz[il-1] = 0.5*dz[il-1]
                z[il-1] = z[il-1]-0.5*dz[il-1]
                temp[il] = temp[il-1]
                grainsize[il] = grainsize[il-1]
                water[il] = water[il-1]
                ice[il] = ice[il-1]
                mass[il] = mass[il-1]
                dens[il] = dens[il-1]
                if (dens[il] > densice): dens[il] = densice
                dz[il] = dz[il-1]
                z[il] = z[il-1] + dz[il]
                lid[il] = lid[il-1]
                if (lid[il] > 0): nlsnow = nlsnow + 1
                nl = nl + 1
                il = il + 1
        else:
            if (check == 5): # first merge and then split in multiple layers
                for ik in range(iks,ike+1):
                    temp_tmp = temp_tmp + temp_old[ik]*dz_old[ik]
                    grainsize_tmp = grainsize_tmp + grainsize_old[ik]*dz_old[ik]
                    water_tmp = water_tmp + water_old[ik]
                    ice_tmp = ice_tmp + ice_old[ik]
                    mass_tmp = mass_tmp + mass_old[ik]
                    dz_tmp = dz_tmp + dz_old[ik]

                temp_tmp = temp_tmp/dz_tmp
                grainsize_tmp = grainsize_tmp/dz_tmp
                dens_tmp = mass_tmp/dz_tmp
                if (dens_tmp > densice): dens_tmp = densice
                z_tmp = z_old[ike]+0.5*dz_old[ike]-0.5*dz_tmp
                lid_tmp = lid_old[ike]
            else: # split layers
                temp_tmp = temp_old[ike]
                grainsize_tmp = grainsize_old[ike]
                water_tmp = water_old[ike]
                ice_tmp = ice_old[ike]
                mass_tmp = mass_old[ike]
                dens_tmp = dens_old[ike]
                dz_tmp = dz_old[ike]
                z_tmp = z_old[ike]
                lid_tmp = lid_old[ike]

            depths = depth
            depthe = depth*(1.+coeff) + dz0
            for ill in range(0,ilp+2):
                ddepth = depthe-depths
                temp[il] = temp_tmp
                grainsize[il] = grainsize_tmp
                water[il] = water_tmp*ddepth/dz_tmp
                ice[il] = ice_tmp*ddepth/dz_tmp
                mass[il] = mass_tmp*ddepth/dz_tmp
                dens[il] = dens_tmp
                dz[il] = ddepth
                lid[il] = lid_old[ike]
                if (il == 0):
                    z[il] = 0.5*ddepth
                else:
                    z[il] = z[il-1] + 0.5*ddepth + 0.5*dz[il-1]
                depth = depthe
                if (lid[il] > 0):
                    nlsnow = nlsnow + 1
                depths = depthe
                depthe = depths*(1.+coeff) + dz0
                if (ill == ilp): depthe = (z_tmp+0.5*dz_tmp)
                il = il + 1
                nl = nl + 1

            temp_tmp = 0.
            grainsize_tmp = 0.
            water_tmp = 0.
            ice_tmp = 0.
            mass_tmp = 0.
            dz_tmp = 0.

    
    if (z[nl-1]+0.5*dz[nl-1] < (zdeep - dzdeep)):		#add layer to maintain sufficient deep model
        nl = nl + 1
        dz[nl-1] = zdeep - z[nl-2]+0.5*dz[nl-2]
        z[nl-1] = z[nl-2]+0.5*dz[nl-2]+0.5*dz[nl-1]
        temp[nl-1] = temp[nl-2]
        grainsize[nl-1] = grainsize[nl-2]
        water[nl-1] = 0.
        ice[nl-1] = 0.
        lid[nl-1] = lid[nl-2]
        dens[nl-1] = dens[nl-2]
        if (dens[nl-2] >= densice):
            dens[nl-1] = densice
            lid[nl-1] = 0
        mass[nl-1] = dens[nl-1]*dz[nl-1]
        if (lid[nl-2] == 1): lid[nl-1] = 2
        if (lcomment == 1):
            print('REDEFGRID: added layer at bottom',nl,lid[nl-1],z[nl-1],dz[nl-1],z[nl-2],dz[nl-2],dzdeep,dens[nl-1])
        # IF (lcomment == 1) WRITE(*,'(a,2i6,6f12.5)') 'REDEFGRID: added layer at bottom',nl,lid[nl-1],&
        # &                    z[nl-1],dz[nl-1],z(nl-1),dz(nl-1),dzdeep,dens[nl-1]
        
    for il in range(0,nl):
        kice[il] = conduc(dens[il],temp[il])
        cpice[il] = 152.5 + 7.122 * temp[il]
        if (il==0):
            dtdz[il] = (temp[il]-t0)/dz[il]
        else: 
            dzl=0.5*(dz[il]+dz[il-1])
            dtdz[il] = (temp[il]-temp[il-1])/dzl
        energy[il] = dens[il]*dz[il]*cpice[il]*(Tkel-temp[il])
    #kice(1) = 50.0
    # if lcomment == 1:
        # print('REDEFGRID: Redefine grid ', nl_old,nl,nlsnow_old,nlsnow,dz[0],dz_old[0],z_old[0],vink,check)
    vink = 0
    
    return water, ice, mass, dens, lid, dz, z, temp, grainsize, energy, kice, cpice, dtdz, vink, nl, dsnowr, nlsnow
