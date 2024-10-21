chstation =  "ant_aws14"
version =  "1D"
lcomment = 1                         # Logical [1] = yes or [0] = no comments to screen
lerrorgap = 0                        # jump over data gaps = 0 or interpolate = 1 or interpolate and use all data for averaging = 2 (errorflag)
lhourlysnowout = 0                   # write snow characteristics every [0] day or [1] hour, default is 0 because it slows down enormously, use only when hourly values are needed for evaluation
lwritelayers = 1                     # write subsurface layers to output [1] or not [0]. Use [0] for faster 2D simulations 
ImpExp = 2                           # Euler explicit [1] or Implicit [2] method to solve for heat diffusion. When [2] , full implicit or Crank-Nicolson may be adapted in the code
tstep = 1800                         # Time step, choose small enough to keep subsurface calculations stable in case ImpExp = 1. In case ImpExp = 2, be mindful of the accuracy of implicit solver
dz0 = 0.01                           # Size of uppermost grid cell (m)
dzdeep = 2                           # Size of lowest  grid cell (m)
zdeep = 10                           # Depth of grid (m)
llake = 0                            # Switch to keep activate lake routine [0]
densice = 916.                       # ice density (kg/m3)
densfirn = 700.                      # firn density (kg/m3)
densclosure = 920.                   # firn pore closure density after which no liquid water is allowed to percolate (kg/m3)
rhosninit = 400.0                    # Density of initial snow pack (kg/m3)
rhosnprec = 400.0                    # Density of fresh snowfall (kg/m3)
trhoprec = 0                         # type of determination of snow fall density 0 = constant, 1 is f(ws,t) (in case 1, rhosnprec is used as lower limit)
tpcond = 10                          # Type of effective conductivity equation 1 = Von Dussen, 2 = Sturm, 3=Douville, 4=Jansson, 5=Anderson, 6=Ostin&Anderson, 7=Arthern, 8=Cox, 9=fixed value, 10= Calonne2019 
T10m = -14.1908                      # Measured 10 m snow temperature (C) 
dsnow_glob = 10.                     # Thickness of snow layer (1 year), at start or on 1 Jan. from acc in m over ice/firn surface
dfirn = 10.                          # Thickness of firn layer, (if ice deeper than zdeep, set to zdeep), at start or on 1 Jan. from acc in m over ice/firn surface
luseacc = 3                          # Use sonic data for: 0= not at all, 1=only for evaluation of ice melt and accumulation, 2= 1 + restrict accum based on precip with measured sonic altim.
tpdens = 0                           # Type of densification routine, (tpdens: 0 = no densification, 1 = Herron and Langway 1980, 2 = Li and Zwally 2004, 3 = Li and Zwally plus vapor transport, 4 = Helsen 2008, 5 = Arthern 2010, 6 = Ligtenberg 2011, -1 = no dry densification, density profile reset after each time step, snow height folows observations)
lrefr = 1                            #  lrefr:1/0 yes/no refreezing
tpirre = 3                           # Type of irreducible water content routine 0 = no water percolation, 1 = constant value cirre (fraction), 2 = Coleou and Lesaffre 1998, 3 = Schneider and Jansson 2004
cirre = 0.02                         # constant value of  irreducible water content if tpirre = 1
lslush = 0                           # yes (1) or no (0) slush formation possible based on Zuo and Oerlemans 1996.
surfangle = 0                        # parameter in slush formation routine : surface angle
tausteep = 0.05                      # parameter in slush formation routine :  runoff time scale of steep slopes > 5
tauhor = 20                          # parameter in slush formation routine : runoff time scale od horizontal slopes
tau1 = 2                             # parameter in slush formation routine : runoff time scales on 1deg slopes
slfact = 10                          # parameter in slush formation routine : factor difference between surface runoff and within snowpack
accyear = 0.5                        # Annual average accumulation for densification calculation (m w.e.)
lz0m = 1                             # Switch roughness length for momentum from file (0), set below (1) or a random value between zll and zul (2) or ice (3), or parameterized using Hmax based on Van Tiggelen et al (2021) (4) 
zll = 0.0005                         #	
zul = 0.05                           #
z0msn = 0.0001                       # Roughness length for momentum for snow in case lz0m = 1 (m)
z0mice = 0.001                       # Roughness length for momentum for ice in case lz0m = 1 (m)
Hmax_info  = 0.2                     # Maximum height of ice obstacles, used in case lz0m = 4
lz0h = 3                             # Switch roughness length for heat and moisture, 0.1*z0m (0), using Andreas 1987 parameterization (1), using Smeets and van den Broeke 2008 parameterization (2) or using Van Tiggelen et al (2021) parameterization (3). For smooth surfaces (z0m < 1mm) and lz0m > 0 Andreas 1987 is used.
tcalc = 4                            # Formulation calculation surface temperature 1 = from Lout observations, 2 = equal to temp uppermost layer, 3 = extrapolated from upper most layers, 4 = skin layer formulation
extrapolation = 1                    # extrapolation to surface temp 1 = based on only upper layer, 2 = based on upper 2 layers
lsnet = 1                            # switch calculate Snet from Sin and Sout (0), or calculate Snet and Sin based on Sout and running mean albedo (1) 
albmin = 0.20                        # minimum albedo used to set error values to and prevent unrealistic values
albmax = 0.96                        # maximum albedo used to set error values to and prevent unrealistic values
emis = 1                             # Emissivity of the snow/ice surface 
lalbedo = 0                          # albedo from 0 = observations,  or parameterised: 1 = Oerlemans and Knap (1998), 2 = Bougamont et al. (2005), 3 = Douville et al. (1995), 4 = Kuipers Munneke et al. (2011),
solzenyes = 0                        #  yes(1) or no(0) correction on albedo for zenith angle, fresh snow specific surface area and refrozen snow grain radius (for lalbedo option 4)
SSAfresh = 60                        # fresh snow specific surface area (for lalbedo option 4)
radrefr = 0.001                      # refrozen snow grain radius (for lalbedo option 4)
albsnow = 0.8                        # albedo snow (used by lalbedo options 1, 2, 3)
albice = 0.5                         # albedo ice(used by lalbedo options 1, 2, 3)
albfirn = 0.6                        #  albedo firn/old snow (used by lalbedo options 1, 2, 3)
soot = 0.0                           # soot concentration in ppmw (used by lalbedo option 4)
snowstar = 1.2                       # characteristic depth scale snow (mm w.e.)(used by lalbedo options 1, 2, 3)
tstarwet = 10                        # time scales (days) of decay for wet snow (used by lalbedo options 1, 2)
tstardry0 = 40                       # dry surface at 0C(used by lalbedo option 2)
tstardry10 = 2000                    # dry surface at -10C (used by lalbedo option 2)
penetration = 0                      # radiation penetration 0 = off, 1 = on
dzrad = 0.001                        # skin layer thickness for radiation penetration calculation (m)
zradmax = 7                          # max depth of radiation penetration (m)
radiussn = 1e-4                      # radius of snow particles for radiation penetration routine
radiusice = 2.5e-3                   # radius of ice particles for radiation penetration routine
lwcloud = 1                          # cloud cover from lwin (1) or not (0), if not then set to 0.
lwmax = 314.0252,  3.7240,  0.0046   # polynomial to describe upper limit of lwin as function of temperature (C), corresponding to max cloud cover
lwmin = 219.3003,  2.2600,  0.0046   # polynomial to describe lower limit of lwin as function of temperature (C), corresponding to min cloud 
lclimtemp = 0                        # switch climate sensitivity test temperature
climtemp = 1.0                       # climate sensitivity test temperature change in K
lclimrad = 0                         # switch climate sensitivity test shortwave radiation down
climrad = 5.0                        # climate sensitivity test shortwave radiation down in %
lclimws = 0                          # switch climate sensitivity test wind speed change
climws = 5.0                         # climate sensitivity test shortwave wind speed changein %
lclimprec = 0                        # switch climate sensitivity test precipitation change
climprec = 10.0                      # climate sensitivity test precipitation change in %