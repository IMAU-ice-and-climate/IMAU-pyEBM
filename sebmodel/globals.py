"""
    Definition of global variables
"""

northing = 'rlat'
easting = 'rlon'
full_field = False
lat = 68.0
lon = -50.0
compression_level = 2                                       # Choose value between 1 and 9 (highest compression)
#-----------------------------------
# PARALLELIZATION 
#-----------------------------------
name = 'Worker'
slurm_use = False                                         # use SLURM
port = 8786
cores = 32                                   # One grid point per core, do not change
nodes = 1                               # grid points submitted in one sbatch script
memory_per_process = 3
memory = memory=str(cores * memory_per_process) + 'GB'                          # memory per processes in GB
shebang = '#!/bin/bash'
slurm_parameters = ['--error=slurm.err','--output=slurm.out','--time=1-00:00:00',  ]
slurm_local_directory = '/perm/rumv/sebmodel/seb-model_python/dask-worker-space'


rowmax = 100000        # max amount of lines, 1 year of half hourly data
colmax = 20            # max amount of input parameters  = 9 + year + jday + hour + some extra for accumulation and running means
nlmax = 200           # max number of layers in the snow 
mmax = 35              # max number of daily, montly and annual parameters into output
secday = 86400         # seconds in 1 day
Tkel = 273.15          # 0 degC in degK
rd = 287.05            # gas constant of dry air (J/Kkg)
rv = 461.51            # gas constant of moist air (J/Kkg)
rstar = 8.314          # universal gas constant (JK-1mol-1)
eps = 0.622            # Rd/Rv
lm = 3.34e5            # latent heat of melting (Jkg-1)
ls = 2.834e6           # latent heat of sublimation (Jkg-1)
lv = 2.501e6           # latent heat of vaporisation (Jkg-1)
beta = 2317            # constant for calculation es (J K-1 kg-1)
es0 = 610.78           # water vapour pressure at melting point Pa (Triple point pressure)
densnow = 150.0        # density of snow, lower limit density snow layer (g/kg=kg/m3)
denswater = 1000.0     # density of water (g/kg=kg/m3)
betadens = 	.0         # 3 constant for calculation densification of snow (tuning), used now as minimum value
tempcut = 272.5        # cut off temperature for gravitational densification
Deff = 1.1e-4          # Vapor diffusion coefficient (m2s-1)
StefBoltz = 5.67e-8    # Stephan Boltzman constant
karman = 0.4           # Von Karman constant
g = 9.81               # gravity accelaration
cp = 1005.0            # specific heat air at constant pressure (J K-1 kg-1)                               
loutmax = 315.0        # 315 W/m2 to limit the outgoing longwave radiation to melting level
tinterv = 5.0          # search interval for t0 (K)
taccur = 0.005         # accuracy of final t0 (K)
Lmocrit = 0.01         # accuracy of final MO-length (relative)
accur = 1.e-10         # accuracy in general
errorval = -908.0      # errorvalue set when no EB is calculated
maxerror = 32.0                                       
snowthreshold = 0.001  # Threshold for cumulated snow to be added as snowfall event
bandmax = 118		#amount of lines to make up solar spectrum


