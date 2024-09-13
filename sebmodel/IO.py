"""
 Read forcing and write output to netcdf file
"""

''' Standard library imports '''
import os
import sys
import xarray as xr
import pandas as pd
import numpy as np
import time
from datetime import timedelta
import configparser

''' Local imports '''
from globals import *
# from variables import *
from info import * 








class IOClass:

    def __init__(self, FORCING=None):
        """ Init IO Class"""
        
        # Read variable list from file
        config = configparser.ConfigParser()
        config.read('./variables')
        self.output = config['vars']['output']
        self.full = config['vars']['full']
        
        # Initialize FORCING
        self.FORCING = FORCING
        self.OUTPUT = None
      
        # If local IO class is initialized we need to get the dimensions of the dataset
        if FORCING is not None:
            self.time = self.FORCING.dims['time']
            self.date = self.FORCING.dims['date']


    #==============================================================================
    # Creates the input FORCING and reads the restart file, if necessary. The function
    # returns the FORCING xarray dataset which contains all input variables.
    #==============================================================================
    def create_forcing_file(self):
        """ Returns the FORCING xarray dataset"""
        
        self.init_forcing_dataset()  #  read FORCING according to the times defined in info.py
        self.time = self.FORCING.dims['time']
        self.date = self.FORCING.dims['date']
        if version == '1D':
            self.ny = self.FORCING.dims['lat']
            self.nx = self.FORCING.dims['lon']
        elif version == '2D':
            self.ny = self.FORCING.dims['y']
            self.nx = self.FORCING.dims['x']

        return self.FORCING
    #==============================================================================
    # The functions create_output_file creates and initializes the OUTPUT xarray FORCIdatasetNGset
    #==============================================================================
    #----------------------------------------------
    # Creates the output xarray dataset
    #----------------------------------------------
    def create_output_file(self):
        """ Returns the OUTPUT xarray """
        self.init_output_dataset()
        return self.OUTPUT
    
    #==============================================================================
    # The init_forcing_dataset read the input NetCDF dataset and stores the FORCING in
    # an xarray. 
    #==============================================================================
    #----------------------------------------------
    # Reads the input FORCING into a xarray dataset 
    #----------------------------------------------
    def init_forcing_dataset(self):
        """     
        WS        ::   Wind speed (at height zm over ground) [m s-1]
        Sin       ::   Incoming shortwave radiation [W m-2]
        Sout      ::   Outgoing shortwave radiation [W m-2]
        Lin       ::   Incoming longwave radiation [W m-2] 
        Lout      ::   Outgoing longwave radiation [W m-2]
        T         ::   Air temperature (at height zt over ground) [C]
        q         ::   Specific humidity (at height zt over ground) [g kg-1]
        P         ::   Air pressure (at height zt over ground) [hPa]
        alb       ::   Surface broadband albedo [-]
        zenith    ::   solar zenith angle [rad]
        Serie     ::   relative change in surface height [m]
        precip    ::   precipitation [kg m-2 s-1]
        zt        ::   height of temperature/humidity sensor above surface [m]
        zm        ::   height of wind speed sensor above surface [m]
        z0m       ::   aerodynamic roughness length [m]
        """
    
        # Open input dataset
        # self.FORCING = xr.open_dataset(os.path.join(FORCING_path,'input',input_netcdf))
        def read_forcing(file,version):
            if version == '1D':
                df = pd.read_csv(file,delim_whitespace=True,header=0)
                # date = pd.to_datetime(df['Date'], format='%Y-%m-%d')
                # hour = pd.to_datetime(df['Hour'], format='%H:%M:%S')
                df['time'] = pd.to_datetime(df['Date'] + ' ' + df['Hour'])
                ds = df.set_index('time').to_xarray()
                ds.coords['lon'] = np.array([lon])
                ds.coords['lat'] = np.array([lat])
            elif version == '2D':
                ds = xr.open_dataset(file)
                # ds = ds.rename_vars({northing: 'y',easting: 'x'}, inplace=True)
                ds = ds.rename({northing: 'y',easting: 'x'})

            return ds
        
        if version == '1D':
            ifile = '../input/' + chstation + '/' + chstation +'_HOUR-EBM.txt'
        elif version == '2D':
            ifile = '../input/' + chstation + '/' + chstation +'_HOUR-EBM.nc'
        self.FORCING = read_forcing(ifile,version)
        self.FORCING['time'] = np.sort(self.FORCING['time'].values)
        self.FORCING['date'] = np.arange(np.datetime64(self.FORCING['time'].values[0]), np.datetime64(self.FORCING['time'].values[-1]), np.timedelta64(1, "D"))
        start_interval = str(self.FORCING.time.values[0])[0:16]
        end_interval = str(self.FORCING.time.values[-1])[0:16]
        time_steps = str(self.FORCING.dims['time'])
        
        print('\n Maximum available time interval from %s until %s. Time steps: %s \n\n' % (start_interval, end_interval, time_steps))

        print('--------------------------------------------------------------')
        print('Checking input FORCING .... \n')

        # Define a auxiliary function to check the validity of the FORCING
        def check(field, max, min):
            '''Check & correct the validity of the input FORCING '''
            if np.nanmax(field) > max or np.nanmin(field) < min:
                print('WARNING %s in input out of range MAX: %.2f MIN: %.2f \n' % (str.capitalize(field.name), np.nanmax(field), np.nanmin(field)))
            xr.where(field < min,field,min)
            xr.where(field > max,field,max)
            
        # Check if FORCING is within valid bounds
        if ('T' in self.FORCING):
            print('Temperature FORCING (T) present')
            check(self.FORCING.T, 50.0, -90.0)
        if ('q' in self.FORCING):
            print('Specific humidity FORCING (q) present')
            check(self.FORCING.q, 100.0, 0.0)
        if ('Sin' in self.FORCING):
            print('Incoming hortwave FORCING (Sin) present')
            check(self.FORCING.Sin, 1600.0, 0.0)
        if ('Sout' in self.FORCING):
            print('Outgoing shortwave FORCING (Sout) present')
            check(self.FORCING.Sin, 1600.0, 0.0)
        if ('WS' in self.FORCING):
            print('Wind velocity FORCING (WS) present')
            check(self.FORCING.WS, 50.0, 0.0)
        if ('precip' in self.FORCING):
            print('Precipitation FORCING (precip) present')
            check(self.FORCING.precip, 40.0, 0.0)
        if ('P' in self.FORCING):
            print('Pressure FORCING (P) present')
            check(self.FORCING.P, 1500.0, 500.0)
        if ('alb' in self.FORCING):
            print('Albedo FORCING (alb) present')
            check(self.FORCING.alb, 1.0, 0.0)
        if ('Lin' in self.FORCING):
            print('Incoming longwave FORCING (Lin) present')
            check(self.FORCING.Lin, 400.0, 100.0)
        if ('Lout' in self.FORCING):
            print('Outgoing longwave FORCING (Lout) present')
            check(self.FORCING.Lout, 400.0, 100.0)
        if ('Serie' in self.FORCING):
            print('Serie FORCING (Serie) present')
            check(self.FORCING.Serie, 100.0, -100.0)
        if ('z0m' in self.FORCING):
            print('Aerodynamic roughness FORCING (z0m) present')
            check(self.FORCING.z0m, 10., 0.0000001)
        if ('zt' in self.FORCING):
            print('Height of T-q sensor FORCING (zt) present')
            check(self.FORCING.zt, 20., 0.1)
        if ('zm' in self.FORCING):
            print('Height of WS sensor FORCING (zm) present')
            check(self.FORCING.zm, 20.0, 0.1)
            
    #==============================================================================
    # The init_output_dataset sets creates the output dataset which stores the results
    # from the individual runs. After the dataset has been filled with the
    # runs from all workers the dataset is written to disc.
    #==============================================================================
    #----------------------------------------------
    # Initializes the output xarray dataset
    #----------------------------------------------
    def init_output_dataset(self):
        """ This function creates the result file 
        Args:
            
            self.FORCING    ::  self.FORCING structure 
            
        Returns:
            
            self.OUTPUT  ::  one-dimensional self.RESULT structure"""
        
        # Coordinates
        self.OUTPUT = xr.Dataset()
        self.OUTPUT.coords['time'] = self.FORCING.coords['time']
        self.OUTPUT.coords['date'] = self.FORCING.coords['date']
        #self.OUTPUT.coords['lat'] = self.FORCING.coords['lat']
        #self.OUTPUT.coords['lon'] = self.FORCING.coords['lon']

        # Global attributes from info.py

        # Global attributes from constants.py
        self.OUTPUT.attrs['tstep'] = tstep
        self.OUTPUT.attrs['dz0'] = dz0
        self.OUTPUT.attrs['zdeep'] = zdeep
        """TO-DO: ADD MORE HERE"""

        # Variables given by the input dataset
        # if 1D version
        # self.add_variable_along_time(self.OUTPUT, self.FORCING.T, 'T', 'C', 'Air temperature')
        # self.add_variable_along_time(self.OUTPUT, self.FORCING.Lin, 'Lin', 'W m-2', 'Incoming longwave radiation')
        
        """ 
        if 2D version
        #self.add_variable_along_latlon(self.OUTPUT, self.FORCING.HGT, 'HGT', 'm', 'Elevation')
        #self.add_variable_along_latlon(self.OUTPUT, self.FORCING.MASK, 'MASK', 'boolean', 'Glacier mask')
        if ('SLOPE' in self.FORCING):
            self.add_variable_along_latlon(self.OUTPUT, self.FORCING.SLOPE, 'SLOPE', 'degrees', 'Terrain slope')
        if ('ASPECT' in self.FORCING):
            self.add_variable_along_latlon(self.OUTPUT, self.FORCING.ASPECT, 'ASPECT', 'degrees', 'Aspect of slope')
        self.add_variable_along_latlontime(self.OUTPUT, self.FORCING.T, 'T', 'C', 'Air temperature')
        self.add_variable_along_latlontime(self.OUTPUT, self.FORCING.q, 'q', '%', 'Specific humidity')
        self.add_variable_along_latlontime(self.OUTPUT, self.FORCING.WS, 'WS', 'm s-1', 'Wind velocity')
        self.add_variable_along_latlontime(self.OUTPUT, self.FORCING.P, 'P', 'hPa', 'Air pressure')
        self.add_variable_along_latlontime(self.OUTPUT, self.FORCING.Sin, 'Sin', 'W m-1', 'Incoming shortwave radiation')
        self.add_variable_along_latlontime(self.OUTPUT, self.FORCING.Sout, 'Sout', 'W m-1', 'Outgoing shortwave radiation')
        self.add_variable_along_latlontime(self.OUTPUT, self.FORCING.Lin, 'Lin', 'W m-1', 'Incoming longwave radiation')
        self.add_variable_along_latlontime(self.OUTPUT, self.FORCING.Lout, 'Lout', 'W m-1', 'Outgoing longwave radiation')
        
        if ('precip' in self.FORCING):
            self.add_variable_along_latlontime(self.OUTPUT, self.FORCING.precip, 'precip', 'kg m-2 s-1','precipitation')
        else:
            self.add_variable_along_latlontime(self.OUTPUT, np.full_like(self.FORCING.T, np.nan), 'precip', 'kg m-2 s-1','Total precipitation')
        
        if ('Serie' in self.FORCING):
            self.add_variable_along_latlontime(self.OUTPUT, self.FORCING.Serie, 'Serie', 'm', 'Serie')
        """
        print('\n') 
        print('Output dataset ... ok')
        return self.OUTPUT
  

    #==============================================================================
    # This function creates the global numpy arrays which store the variables.
    # The global array is filled with the local results from the workers. Finally,
    # the arrays are assigned to the OUTPUT dataset and is stored to disc
    #==============================================================================
    def create_global_output_arrays(self):
        if ('Sin' in self.output):
            self.Sin = np.full((self.time,self.ny,self.nx), np.nan)
        if ('Sout' in self.output):
            self.Sout = np.full((self.time,self.ny,self.nx), np.nan)
        if ('Lin' in self.output):
            self.Lin = np.full((self.time,self.ny,self.nx), np.nan)
        if ('Loutobs' in self.output):
            self.Loutobs = np.full((self.time,self.ny,self.nx), np.nan)
        if ('Loutmod' in self.output):
            self.Loutmod = np.full((self.time,self.ny,self.nx), np.nan)   
        if ('SH' in self.output):
            self.SH = np.full((self.time,self.ny,self.nx), np.nan)
        if ('LE' in self.output):
            self.LE = np.full((self.time,self.ny,self.nx), np.nan)
        if ('GH' in self.output):
            self.GH = np.full((self.time,self.ny,self.nx), np.nan)
        if ('Restsource' in self.output):
            self.Restsource = np.full((self.time,self.ny,self.nx), np.nan)
        if ('Source' in self.output):
            self.Source = np.full((self.time,self.ny,self.nx), np.nan)
        if ('sumdivs' in self.output):
            self.sumdivs = np.full((self.time,self.ny,self.nx), np.nan)
        if ('T0' in self.output):
            self.T0 = np.full((self.time,self.ny,self.nx), np.nan)
        if ('q0' in self.output):
            self.q0 = np.full((self.time,self.ny,self.nx), np.nan)   
        if ('T2m' in self.output):
            self.T2m = np.full((self.time,self.ny,self.nx), np.nan)
        if ('q2m' in self.output):
            self.q2m = np.full((self.time,self.ny,self.nx), np.nan)
        if ('WS10m' in self.output):
            self.WS10m = np.full((self.time,self.ny,self.nx), np.nan)
        if ('z0m' in self.output):
            self.z0m = np.full((self.time,self.ny,self.nx), np.nan)
        if ('icemelt' in self.output):
            self.icemelt = np.full((self.time,self.ny,self.nx), np.nan)
        if ('dsnowacc' in self.output):
            self.dsnowacc = np.full((self.time,self.ny,self.nx), np.nan)
        if ('hsnowmod' in self.output):
            self.hsnowmod = np.full((self.time,self.ny,self.nx), np.nan)
        if ('runoff' in self.output):
            self.runoff = np.full((self.time,self.ny,self.nx), np.nan)
        if ('runoffdt' in self.output):
            self.runoffdt = np.full((self.time,self.ny,self.nx), np.nan)
        if ('surfwater' in self.output):
            self.surfwater = np.full((self.time,self.ny,self.nx), np.nan)
        if ('melt' in self.output):
            self.melt = np.full((self.time,self.ny,self.nx), np.nan)
        if ('meltdt' in self.output):
            self.meltdt = np.full((self.time,self.ny,self.nx), np.nan)
        if ('surfmelt' in self.output):
            self.surfmelt = np.full((self.time,self.ny,self.nx), np.nan)
        if ('surfmeltdt' in self.output):
            self.surfmeltdt = np.full((self.time,self.ny,self.nx), np.nan)
        if ('sumdrift' in self.output):
            self.sumdrift = np.full((self.time,self.ny,self.nx), np.nan)
        if ('subl' in self.output):
            self.subl = np.full((self.time,self.ny,self.nx), np.nan)
        if ('subldt' in self.output):
            self.subldt = np.full((self.time,self.ny,self.nx), np.nan)
        if ('icemeltmdt' in self.output):
            self.icemeltmdt = np.full((self.time,self.ny,self.nx), np.nan)
        if ('precip' in self.output):
            self.precip = np.full((self.time,self.ny,self.nx), np.nan)
        if ('precipdt' in self.output):
            self.precipdt = np.full((self.time,self.ny,self.nx), np.nan)
        if ('dens_lay1' in self.output):
            self.dens_lay1 = np.full((self.time,self.ny,self.nx), np.nan)
        if ('temp_lay1' in self.output):
            self.temp_lay1 = np.full((self.time,self.ny,self.nx), np.nan)
        if ('dz_lay1' in self.output):
            self.dz_lay1 = np.full((self.time,self.ny,self.nx), np.nan)
        if ('sumwater' in self.output):
            self.sumwater = np.full((self.time,self.ny,self.nx), np.nan)
        if ('topwater' in self.output):
            self.topwater = np.full((self.time,self.ny,self.nx), np.nan)
        if ('air_content' in self.output):
            self.air_content = np.full((self.time,self.ny,self.nx), np.nan)
        if ('effective_air_content' in self.output):
            self.effective_air_content = np.full((self.time,self.ny,self.nx), np.nan)

        if ('errorflag' in self.output):
            self.errorflag = np.full((self.time,self.ny,self.nx), np.nan)
            
        if lhourlysnowout == 1:
            if ('temp' in self.full):
                self.temp = np.full((self.time,self.ny,self.nx,nlmax), np.nan)
            if ('dens' in self.full):
                self.dens = np.full((self.time,self.ny,self.nx,nlmax), np.nan)  
            if ('kice' in self.full):
                self.kice = np.full((self.time,self.ny,self.nx,nlmax), np.nan)
            if ('cpice' in self.full):
                self.cpice = np.full((self.time,self.ny,self.nx,nlmax), np.nan)  
            if ('rhocp' in self.full):
                self.rhocp = np.full((self.time,self.ny,self.nx,nlmax), np.nan)
            if ('energy' in self.full):
                self.energy = np.full((self.time,self.ny,self.nx,nlmax), np.nan)  
            if ('z' in self.full):
                self.z = np.full((self.time,self.ny,self.nx,nlmax), np.nan)
            if ('dz' in self.full):
                self.dz = np.full((self.time,self.ny,self.nx,nlmax), np.nan)  
            if ('lid' in self.full):
                self.lid = np.full((self.time,self.ny,self.nx,nlmax), np.nan)  
            if ('mass' in self.full):
                self.mass = np.full((self.time,self.ny,self.nx,nlmax), np.nan)  
            if ('grainsize' in self.full):
                self.grainsize = np.full((self.time,self.ny,self.nx,nlmax), np.nan)  
            if ('water' in self.full):
                self.water = np.full((self.time,self.ny,self.nx,nlmax), np.nan)  
            if ('ice' in self.full):
                self.ice = np.full((self.time,self.ny,self.nx,nlmax), np.nan)  
            if ('dsdz' in self.full):
                self.dsdz = np.full((self.time,self.ny,self.nx,nlmax), np.nan)           
            if ('refrfrac' in self.full):
                self.refrfrac = np.full((self.time,self.ny,self.nx,nlmax), np.nan) 
        elif lhourlysnowout == 0:    
            if ('temp' in self.full):
                self.temp = np.full((self.date,self.ny,self.nx,nlmax), np.nan)
            if ('dens' in self.full):
                self.dens = np.full((self.date,self.ny,self.nx,nlmax), np.nan)  
            if ('kice' in self.full):
                self.kice = np.full((self.date,self.ny,self.nx,nlmax), np.nan)
            if ('cpice' in self.full):
                self.cpice = np.full((self.date,self.ny,self.nx,nlmax), np.nan)  
            if ('rhocp' in self.full):
                self.rhocp = np.full((self.date,self.ny,self.nx,nlmax), np.nan)
            if ('energy' in self.full):
                self.energy = np.full((self.date,self.ny,self.nx,nlmax), np.nan)  
            if ('z' in self.full):
                self.z = np.full((self.date,self.ny,self.nx,nlmax), np.nan)
            if ('dz' in self.full):
                self.dz = np.full((self.date,self.ny,self.nx,nlmax), np.nan)  
            if ('lid' in self.full):
                self.lid = np.full((self.date,self.ny,self.nx,nlmax), np.nan)  
            if ('mass' in self.full):
                self.mass = np.full((self.date,self.ny,self.nx,nlmax), np.nan)  
            if ('grainsize' in self.full):
                self.grainsize = np.full((self.date,self.ny,self.nx,nlmax), np.nan)  
            if ('water' in self.full):
                self.water = np.full((self.date,self.ny,self.nx,nlmax), np.nan)  
            if ('ice' in self.full):
                self.ice = np.full((self.date,self.ny,self.nx,nlmax), np.nan)  
            if ('dsdz' in self.full):
                self.dsdz = np.full((self.date,self.ny,self.nx,nlmax), np.nan)           
            if ('refrfrac' in self.full):
                self.refrfrac = np.full((self.date,self.ny,self.nx,nlmax), np.nan) 
    #==============================================================================
    # This function assigns the local results from the workers to the global
    # numpy arrays. The y and x values are the lat/lon indices.
    #==============================================================================
    def copy_local_to_global(self,y,x,local_Sin,local_Sout,local_Lin,local_Loutobs,local_Loutmod,local_SH,local_LE,
        local_GH, local_Restsource, local_Source, local_sumdivs,local_T0, local_q0, local_T2m, local_q2m, local_WS10m, local_z0m,
            local_z, local_dz, local_temp, local_dens, local_kice, local_cpice, local_rhocp, local_energy,local_lid,local_mass,
            local_grainsize, local_water, local_ice,local_dsdz,local_refrfrac,
             local_icemelt, local_icemeltmdt, local_dsnowacc,local_hsnowmod,local_runoff,local_runoffdt, 
             local_surfwater, local_melt,local_meltdt, local_surfmelt, local_surfmeltdt,local_sumdrift, 
             local_subl, local_subldt, local_precip,local_precipdt,local_dens_lay1,
             local_temp_lay1,local_dz_lay1,local_sumwater,local_topwater,local_air_content,local_effective_air_content,local_errorflag):


        if ('Sin' in self.output):
            self.Sin[:,y,x] = local_Sin
        if ('Sout' in self.output):
            self.Sout[:,y,x] = local_Sout
        if ('Lin' in self.output):
            self.Lin[:,y,x] = local_Lin
        if ('Loutobs' in self.output):
            self.Loutobs[:,y,x] = local_Loutobs
        if ('Loutmod' in self.output):
            self.Loutmod[:,y,x] = local_Loutmod
        if ('SH' in self.output):
            self.SH[:,y,x] = local_SH
        if ('LE' in self.output):
            self.LE[:,y,x] = local_LE
        if ('GH' in self.output):
            self.GH[:,y,x] = local_GH
        if ('Restsource' in self.output):
            self.Restsource[:,y,x] = local_Restsource
        if ('Source' in self.output):
            self.Source[:,y,x] = local_Source
        if ('sumdivs' in self.output):
            self.sumdivs[:,y,x] = local_sumdivs
        if ('T0' in self.output):
            self.T0[:,y,x] = local_T0
        if ('q0' in self.output):
            self.q0[:,y,x] = local_q0
        if ('T2m' in self.output):
            self.T2m[:,y,x] = local_T2m
        if ('q2m' in self.output):
            self.q2m[:,y,x] = local_q2m
        if ('WS10m' in self.output):
            self.WS10m[:,y,x] = local_WS10m
        if ('z0m' in self.output):
            self.z0m[:,y,x] = local_z0m
        if ('icemelt' in self.output):
            self.icemelt[:,y,x] = local_icemelt
        if ('dsnowacc' in self.output):
            self.dsnowacc[:,y,x] = local_dsnowacc
        if ('hsnowmod' in self.output):
            self.hsnowmod[:,y,x] = local_hsnowmod
        if ('runoff' in self.output):
            self.runoff[:,y,x] = local_runoff
        if ('runoffdt' in self.output):
            self.runoffdt[:,y,x] = local_runoffdt
        if ('surfwater' in self.output):
            self.surfwater[:,y,x] = local_surfwater
        if ('melt' in self.output):
            self.melt[:,y,x] = local_melt
        if ('meltdt' in self.output):
            self.meltdt[:,y,x] = local_meltdt
        if ('surfmelt' in self.output):
            self.surfmelt[:,y,x] = local_surfmelt
        if ('surfmeltdt' in self.output):
            self.surfmeltdt[:,y,x] = local_surfmeltdt
        if ('sumdrift' in self.output):
            self.sumdrift[:,y,x] = local_sumdrift
        if ('subl' in self.output):
            self.subl[:,y,x] = local_subl
        if ('subldt' in self.output):
            self.subldt[:,y,x] = local_subldt
        if ('icemeltmdt' in self.output):
            self.icemeltmdt[:,y,x] = local_icemeltmdt
        if ('precip' in self.output):
            self.precip[:,y,x] = local_precip
        if ('precipdt' in self.output):
            self.precipdt[:,y,x] = local_precipdt
        if ('dens_lay1' in self.output):
            self.dens_lay1[:,y,x] = local_dens_lay1
        if ('temp_lay1' in self.output):
            self.temp_lay1[:,y,x] = local_temp_lay1
        if ('dz_lay1' in self.output):
            self.dz_lay1[:,y,x] = local_dz_lay1
        if ('sumwater' in self.output):
            self.sumwater[:,y,x] = local_sumwater
        if ('topwater' in self.output):
            self.topwater[:,y,x] = local_topwater
        if ('air_content' in self.output):
            self.air_content[:,y,x] = local_air_content
        if ('effective_air_content' in self.output):
            self.effective_air_content[:,y,x] = local_effective_air_content
        if ('errorflag' in self.output):
            self.errorflag[:,y,x] = local_errorflag        
            
        if lwritelayers == 1:
            if ('temp' in self.full):
                self.temp[:,y,x,:] = local_temp 
            if ('dens' in self.full):
                self.dens[:,y,x,:] = local_dens
            if ('kice' in self.full):
                self.kice[:,y,x,:] = local_kice
            if ('cpice' in self.full):
                self.cpice[:,y,x,:] = local_cpice
            if ('rhocp' in self.full):
                self.rhocp[:,y,x,:] = local_rhocp
            if ('energy' in self.full):
                self.energy[:,y,x,:] = local_energy
            if ('z' in self.full):
                self.z[:,y,x,:] = local_z
            if ('dz' in self.full):
                self.dz[:,y,x,:] = local_dz
            if ('lid' in self.full):
                self.lid[:,y,x,:] = local_lid
            if ('mass' in self.full):
                self.mass[:,y,x,:] = local_mass
            if ('grainsize' in self.full):
                self.grainsize[:,y,x,:] = local_grainsize
            if ('water' in self.full):
                self.water[:,y,x,:] = local_water
            if ('ice' in self.full):
                self.ice[:,y,x,:] = local_ice
            if ('dsdz' in self.full):
                self.dsdz[:,y,x,:] = local_dsdz
            if ('refrfrac' in self.full):
                self.refrfrac[:,y,x,:] = local_refrfrac
    #==============================================================================
    # This function adds the global numpy arrays to the OUTPUT dataset which will
    # be written to disc.
    #==============================================================================
    def write_outputs_to_file(self):
        # self.add_variable_along_time(self.OUTPUT, self.Lin, 'Lin', 'W m-2', 'Incoming longwave radiation')
        if ('Sin' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.Sin, 'Sin', 'W m-2', 'Incoming shortwave radiation') 
        if ('Sout' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.Sout, 'Sout', 'W m-2', 'Outgoing shortwave radiation') 
        if ('Lin' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.Lin, 'Lin', 'W m-2', 'Incoming longwave radiation') 
        if ('Loutobs' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.Loutobs, 'Loutobs', 'W m-2', 'Outgoing longwave radiation observed') 
        if ('Loutmod' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.Loutmod, 'Loutmod', 'W m-2', 'Outgoing longwave radiation modelled') 
        if ('SH' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.SH, 'SH', 'W m-2', 'Sensible heat flux positive downwards') 
        if ('LE' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.LE, 'LE', 'W m-2', 'Latent heat flux positive downwards') 
        if ('GH' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.GH, 'GH', 'W m-2', 'Ground heat flux') 
        if ('Restsource' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.Restsource, 'Restsource', 'W m-2', 'Residual energy from SEB') 
        if ('Source' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.Source, 'Source', 'W m-2', 'Energy available for surface melt') 
        if ('sumdivs' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.sumdivs, 'sumdivs', 'W m-2', 'Total of shortwave radiation penetrated in the snow') 
        if ('T0' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.T0, 'T0', 'K', 'Surface temperature modelled') 
        if ('q0' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.q0, 'q0', 'kg kg-1', 'Surface specific humidity') 
        if ('T2m' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.T2m, 'T2m', 'K', '2m air temperature') 
        if ('q2m' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.q2m, 'q2m', 'kg kg-1', '2m specific humidity') 
        if ('WS10m' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.WS10m, 'WS10m', 'm s-1', '10m wind speed') 
        if ('z0m' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.z0m, 'z0m', 'm', 'Aerodynamic roughness length for momentum') 
        if ('icemelt' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.icemelt, 'icemelt', 'm', 'cumulative ice melt in m ice') 
        if ('dsnowacc' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.dsnowacc, 'dsnowacc', 'm', 'Snow height') 
        if ('hsnowmod' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.hsnowmod, 'hsnowmod', 'm', 'Cumulative modelled surface height') 
        if ('runoff' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.runoff, 'runoff', 'mm w.e.', 'cumulative runoff') 
        if ('runoffdt' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.runoffdt, 'runoffdt', 'mm w.e. dt^-1', 'runoff per timestep') 
        if ('surfwater' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.surfwater, 'surfwater', 'mm w.e.', 'Liquid water on surface') 
        if ('melt' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.melt, 'melt', 'mm w.e.', 'cumulative melt') 
        if ('meltdt' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.meltdt, 'meltdt', 'mm w.e. dt^-1', 'melt per timestep') 
        if ('surfmelt' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.surfmelt, 'surfmelt', 'mm w.e.', 'cumulative surface melt') 
        if ('surfmeltdt' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.surfmeltdt, 'surfmeltdt', 'mm w.e. dt^-1', 'surface melt per timestep') 
        if ('sumdrift' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.sumdrift, 'sumdrift', 'mm w.e.', 'sumdrift') 
        if ('subl' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.subl, 'subl', 'mm w.e.', 'cumulative sublimation') 
        if ('subldt' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.subldt, 'subldt', 'mm w.e. dt^-1', 'sublimation per time step') 
        if ('icemeltmdt' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.icemeltmdt, 'icemeltmdt', 'm', 'ice melt in m ice per time step') 
        if ('precip' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.precip, 'precip', 'mm w.e.', 'cumulative precipitation') 
        if ('precipdt' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.precipdt, 'precipdt', 'mm w.e.', 'precipitation per time step') 
        if ('dens_lay1' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.dens_lay1, 'dens_lay1', 'kg m-3', 'density first model level') 
        if ('temp_lay1' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.temp_lay1, 'temp_lay1', 'K', 'temperature first model level') 
        if ('dz_lay1' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.dz_lay1, 'dz_lay1', 'm', 'thickness first model level')
        if ('sumwater' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.sumwater, 'sumwater', 'kg m-2', 'liquid water content')
        if ('topwater' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.topwater, 'topwater', 'kg m-2', 'liquid water content above first ice layer')
        if ('air_content' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.air_content, 'air_content', 'kg m-2', 'firn air content')
        if ('effective_air_content' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.effective_air_content, 'effective_air_content', 'kg m-2', 'effective/reachable firn air content')
                    
        if ('errorflag' in self.output):
            self.add_variable_along_latlontime(self.OUTPUT, self.errorflag, 'errorflag', ' ', 'Flag for error in data. Bad data in forcing when not zero.') 
        if lhourlysnowout == 1:
            if ('temp' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.temp, 'temp', 'K', 'subsurface temperature') 
            if ('dens' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.dens, 'dens', 'kg m^-3', 'subsurface density') 
            if ('kice' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.kice, 'kice', 'K', 'subsurface thermal conductivity') 
            if ('cpice' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.cpice, 'cpice', ' ', 'subsurface heat capacity') 
            if ('rhocp' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.rhocp, 'rhocp', ' ', ' ') 
            if ('energy' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.energy, 'energy', ' ', 'subsurface energy content') 
            if ('z' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.z, 'z', 'm', 'subsurface layer depth') 
            if ('dz' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.dz, 'dz', 'm', 'subsurface layer thickness') 
            if ('lid' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.lid, 'lid', '-', 'layer id 0 = ice (for glacier purposes), 1 = snow 2 = firn + older firn') 
            if ('mass' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.mass, 'mass', 'kg', 'mass layer (dens * dz)') 
            if ('grainsize' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.grainsize, 'grainsize', ' ', 'snow grain size') 
            if ('water' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.water, 'water', ' ', 'water content layer') 
            if ('ice' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.ice, 'ice', ' ', 'refrozen water content layer') 
            if ('dsdz' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.dsdz, 'dsdz', ' ', 'radiation penetration') 
            if ('refrfrac' in self.full):
                self.add_variable_along_latlonlayertime(self.OUTPUT, self.refrfrac, 'refrfrac', ' ', 'refrfrac') 
        elif lhourlysnowout == 0:
            if ('temp' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.temp, 'temp', 'K', 'subsurface temperature') 
            if ('dens' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.dens, 'dens', 'kg m^-3', 'subsurface density') 
            if ('kice' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.kice, 'kice', 'K', 'subsurface thermal conductivity') 
            if ('cpice' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.cpice, 'cpice', ' ', 'subsurface heat capacity') 
            if ('rhocp' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.rhocp, 'rhocp', ' ', ' ') 
            if ('energy' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.energy, 'energy', ' ', 'subsurface energy content') 
            if ('z' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.z, 'z', 'm', 'subsurface layer depth') 
            if ('dz' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.dz, 'dz', 'm', 'subsurface layer thickness') 
            if ('lid' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.lid, 'lid', '-', 'layer id 0 = ice (for glacier purposes), 1 = snow 2 = firn + older firn') 
            if ('mass' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.mass, 'mass', 'kg', 'mass layer (dens * dz)') 
            if ('grainsize' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.grainsize, 'grainsize', ' ', 'snow grain size') 
            if ('water' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.water, 'water', ' ', 'water content layer') 
            if ('ice' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.ice, 'ice', ' ', 'refrozen water content layer') 
            if ('dsdz' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.dsdz, 'dsdz', ' ', 'radiation penetration') 
            if ('refrfrac' in self.full):
                self.add_variable_along_latlonlayerdate(self.OUTPUT, self.refrfrac, 'refrfrac', ' ', 'refrfrac') 
            
    #----------------------------------------------
    # Time averaging of output 
    #----------------------------------------------
    def average_output(self,OUTPUT):
        """ This function averages the output NetCDF in time  (Daily)

        Returns:
             self.RESTART  ::  xarray structure
        """
        
        if lhourlysnowout == 1:
            self.OUTPUT_D = OUTPUT
            self.OUTPUT_M = OUTPUT
            self.OUTPUT_Y = OUTPUT
            self.OUTPUT_SEAS = OUTPUT
        else:
            self.OUTPUT_D = OUTPUT.drop_dims('date')
            self.OUTPUT_M = OUTPUT.drop_dims('date')
            self.OUTPUT_Y = OUTPUT.drop_dims('date')
            self.OUTPUT_SEAS = OUTPUT.drop_dims('date')
        
        self.OUTPUT_D = self.OUTPUT_D.resample(time='1D').mean()
        self.OUTPUT_M = self.OUTPUT_M.resample(time='1M').mean()
        self.OUTPUT_Y = self.OUTPUT_Y.resample(time='1Y').mean()
        self.OUTPUT_SEAS = self.OUTPUT_SEAS.resample(time='1Q').mean()
        # return self.OUTPUT_D
    
    def write_output_M(self,OUTPUT):
        """ This function averages the output NetCDF in time  (Monthly)

        Returns:
             self.RESTART  ::  xarray structure
        """
        self.OUTPUT_M = OUTPUT
        self.OUTPUT_M = self.OUTPUT_M.resample(time='1M').mean()
        return self.OUTPUT_M
    
    def write_output_Y(self,OUTPUT):
        """ This function averages the output NetCDF in time  (Yearly)

        Returns:
             self.RESTART  ::  xarray structure
        """
        self.OUTPUT_Y = OUTPUT
        self.OUTPUT_Y = self.OUTPUT_Y.resample(time='1Y').mean()
        return self.OUTPUT_Y
    
    def write_output_SEAS(self,OUTPUT):
        """ This function averages the output NetCDF in time  (Seasonnal)

        Returns:
             self.RESTART  ::  xarray structure
        """
        self.OUTPUT_SEAS = OUTPUT
        self.OUTPUT_SEAS = self.OUTPUT_SEAS.resample(time='1Q').mean()
        return self.OUTPUT_SEAS
    #----------------------------------------------
    # Initializes the restart xarray dataset
    #----------------------------------------------
    def init_restart_dataset(self):
        """ This function creates the restart file 
            
        Returns:
            
            self.RESTART  ::  xarray structure"""
        
        self.RESTART = xr.dataset()
        self.RESTART.coords['time'] = self.FORCING.coords['time'][-1]
        self.RESTART.coords['lat'] = self.FORCING.coords['lat']
        self.RESTART.coords['lon'] = self.FORCING.coords['lon']
        self.RESTART.coords['layer'] = np.arange(nlmax)
    
        print('Restart dFORCINGset ... ok \n')
        print('--------------------------------------------------------------\n')
        
        return self.RESTART
  

    #==============================================================================
    # This function creates the global numpy arrays which store the profiles.
    # The global array is filled with the local results from the workers. Finally,
    # the arrays are assigned to the RESTART dataset and is stored to disc (see COSIPY.py)
    #==============================================================================
    def create_global_restart_arrays(self):
        self.RES_NLAYERS = np.full((self.ny,self.nx), np.nan)
        self.RES_NEWSNOWHEIGHT = np.full((self.ny, self.nx), np.nan)
        self.RES_NEWSNOWTIMESTAMP = np.full((self.ny, self.nx), np.nan)
        self.RES_OLDSNOWTIMESTAMP = np.full((self.ny, self.nx), np.nan)
        self.RES_LAYER_HEIGHT = np.full((self.ny,self.nx,nlmax), np.nan)
        self.RES_LAYER_RHO = np.full((self.ny,self.nx,nlmax), np.nan)
        self.RES_LAYER_T = np.full((self.ny,self.nx,nlmax), np.nan)
        self.RES_LAYER_LWC = np.full((self.ny,self.nx,nlmax), np.nan)
        self.RES_LAYER_IF = np.full((self.ny,self.nx,nlmax), np.nan)


    #----------------------------------------------
    # Initializes the local restart xarray dataset
    #----------------------------------------------
    def create_local_restart_dataset(self):
        """ This function creates the result dataset for a grid point 
        Args:
            
            self.FORCING    ::  self.FORCING structure 
            
        Returns:
            
            self.RESTART  ::  one-dimensional self.RESULT structure"""
    
        self.RESTART = xr.dataset()
        self.RESTART.coords['time'] = self.FORCING.coords['time'][-1]
        self.RESTART.coords['lat'] = self.FORCING.coords['lat']
        self.RESTART.coords['lon'] = self.FORCING.coords['lon']
        self.RESTART.coords['layer'] = np.arange(nlmax)
        
        self.add_variable_along_scalar(self.RESTART, np.full((1), np.nan), 'NLAYERS', '-', 'Number of layers')
        self.add_variable_along_scalar(self.RESTART, np.full((1), np.nan), 'NEWSNOWHEIGHT', 'm .w.e', 'New snow height')
        self.add_variable_along_scalar(self.RESTART, np.full((1), np.nan), 'NEWSNOWTIMESTAMP', 's', 'New snow timestamp')
        self.add_variable_along_scalar(self.RESTART, np.full((1), np.nan), 'OLDSNOWTIMESTAMP', 's', 'Old snow timestamp')

        self.add_variable_along_layer(self.RESTART, np.full((self.RESTART.coords['layer'].shape[0]), np.nan), 'LAYER_HEIGHT', 'm', 'Layer height')
        self.add_variable_along_layer(self.RESTART, np.full((self.RESTART.coords['layer'].shape[0]), np.nan), 'LAYER_RHO', 'kg m^-3', 'Density of layer')
        self.add_variable_along_layer(self.RESTART, np.full((self.RESTART.coords['layer'].shape[0]), np.nan), 'LAYER_T', 'K', 'Layer temperature')
        self.add_variable_along_layer(self.RESTART, np.full((self.RESTART.coords['layer'].shape[0]), np.nan), 'LAYER_LWC', '-', 'Layer liquid water content')
        self.add_variable_along_layer(self.RESTART, np.full((self.RESTART.coords['layer'].shape[0]), np.nan), 'LAYER_IF', '-', 'Layer ice fraction')


        return self.RESTART


    #==============================================================================
    # This function assigns the local results from the workers to the global
    # numpy arrays. The y and x values are the lat/lon indices.
    #==============================================================================
    def copy_local_restart_to_global(self,y,x,local_restart):
        self.RES_NLAYERS[y,x] = local_restart.NLAYERS
        self.RES_NEWSNOWHEIGHT[y,x] = local_restart.NEWSNOWHEIGHT
        self.RES_NEWSNOWTIMESTAMP[y,x] = local_restart.NEWSNOWTIMESTAMP
        self.RES_OLDSNOWTIMESTAMP[y,x] = local_restart.OLDSNOWTIMESTAMP
        self.RES_LAYER_HEIGHT[y,x,:] = local_restart.LAYER_HEIGHT 
        self.RES_LAYER_RHO[y,x,:] = local_restart.LAYER_RHO
        self.RES_LAYER_T[y,x,:] = local_restart.LAYER_T
        self.RES_LAYER_LWC[y,x,:] = local_restart.LAYER_LWC
        self.RES_LAYER_IF[y,x,:] = local_restart.LAYER_IF

    
    #==============================================================================
    # This function adds the global numpy arrays to the RESULT dataset which will
    # be written to disc.
    #==============================================================================
    def write_restart_to_file(self):
        self.add_variable_along_latlon(self.RESTART, self.RES_NLAYERS, 'NLAYERS', '-', 'Number of layers')
        self.add_variable_along_latlon(self.RESTART, self.RES_NEWSNOWHEIGHT, 'new_snow_height', 'm .w.e', 'New snow height')
        self.add_variable_along_latlon(self.RESTART, self.RES_NEWSNOWTIMESTAMP, 'new_snow_timestamp', 's', 'New snow timestamp')
        self.add_variable_along_latlon(self.RESTART, self.RES_OLDSNOWTIMESTAMP, 'old_snow_timestamp', 's', 'Old snow timestamp')
        self.add_variable_along_latlonlayer(self.RESTART, self.RES_LAYER_HEIGHT, 'LAYER_HEIGHT', 'm', 'Height of each layer')
        self.add_variable_along_latlonlayer(self.RESTART, self.RES_LAYER_RHO, 'LAYER_RHO', 'kg m^-3', 'Layer density')
        self.add_variable_along_latlonlayer(self.RESTART, self.RES_LAYER_T, 'LAYER_T', 'K', 'Layer temperature')
        self.add_variable_along_latlonlayer(self.RESTART, self.RES_LAYER_LWC, 'LAYER_LWC', '-', 'Layer liquid water content')
        self.add_variable_along_latlonlayer(self.RESTART, self.RES_LAYER_IF, 'LAYER_IF', '-', 'Layer ice fraction')


    # TODO: Make it Pythonian - Finish the getter/setter functions
    @property
    def RAIN(self):
        return self.__RAIN
    @property
    def SNOWFALL(self):
        return self.__SNOWFALL
    @property
    def LWin(self):
        return self.__LWin
    @property
    def LWout(self):
        return self.__LWout
    @property
    def H(self):
        return self.__H
    @property
    def LE(self):
        return self.__LE
    @property
    def B(self):
        return self.__B
    @property
    def QRR(self):
        return self.__QRR
    @property
    def MB(self):
        return self.__MB
    
    
    @RAIN.setter
    def RAIN(self, x):
        self.__RAIN = x
    @SNOWFALL.setter
    def SNOWFALL(self, x):
        self.__SNOWFALL = x
    @LWin.setter
    def LWin(self, x):
        self.__LWin = x
    @LWout.setter
    def LWout(self, x):
        self.__LWout = x
    @H.setter
    def H(self, x):
        self.__H = x
    @LE.setter
    def LE(self, x):
        self.__LE = x
    @B.setter
    def B(self, x):
        self.__B = x
    @QRR.setter
    def QRR(self, x):
        self.__QRR = x
    @MB.setter
    def MB(self, x):
        self.__MB = x


    #==============================================================================
    # The following functions return the OUTPUT, RESTART and GRID structures
    #==============================================================================
    #----------------------------------------------
    # Getter/Setter functions 
    #----------------------------------------------
    def get_output(self):
        return self.OUTPUT
    
    def get_output_daily(self):
        return self.OUTPUT_D
    
    def get_output_monthly(self):
        return self.OUTPUT_M
    
    def get_output_yearly(self):
        return self.OUTPUT_Y
    
    def get_output_seas(self):
        return self.OUTPUT_SEAS
    
    def get_restart(self):
        return self.RESTART

    def get_grid_restart(self):
        return self.GRID_RESTART

    #==============================================================================
    # Auxiliary functions for writing variables to NetCDF files
    #==============================================================================
    def add_variable_along_scalar(self, ds, var, name, units, long_name):
        """ This function self.adds missing variables to the self.FORCING class """
        ds[name] = var
        ds[name].attrs['units'] = units
        ds[name].attrs['long_name'] = long_name
        ds[name].encoding['_FillValue'] = -9999
        return ds

    def add_variable_along_latlon(self, ds, var, name, units, long_name):
        """ This function self.adds missing variables to the self.FORCING class """
        ds[name] = ((northing,easting), var)
        ds[name].attrs['units'] = units
        ds[name].attrs['long_name'] = long_name
        ds[name].encoding['_FillValue'] = -9999
        return ds
    
    def add_variable_along_time(self, ds, var, name, units, long_name):
        """ This function self.adds missing variables to the self.FORCING class """
        ds[name] = xr.DataArray(var, coords=[('time', ds.time)])
        ds[name].attrs['units'] = units
        ds[name].attrs['long_name'] = long_name
        ds[name].encoding['_FillValue'] = -9999
        return ds
    
    def add_variable_along_latlontime(self, ds, var, name, units, long_name):
        """ This function self.adds missing variables to the self.FORCING class """
        ds[name] = (('time',northing,easting), var)
        ds[name].attrs['units'] = units
        ds[name].attrs['long_name'] = long_name
        ds[name].encoding['_FillValue'] = -9999
        return ds
    
    def add_variable_along_latlonlayertime(self, ds, var, name, units, long_name):
        """ This function self.adds missing variables to the self.FORCING class """
        ds[name] = (('time',northing,easting,'layer'), var)
        ds[name].attrs['units'] = units
        ds[name].attrs['long_name'] = long_name
        ds[name].encoding['_FillValue'] = -9999
        return ds
    
    def add_variable_along_latlonlayerdate(self, ds, var, name, units, long_name):
        """ This function self.adds missing variables to the self.FORCING class """
        ds[name] = (('date',northing,easting,'layer'), var)
        ds[name].attrs['units'] = units
        ds[name].attrs['long_name'] = long_name
        ds[name].encoding['_FillValue'] = -9999
        return ds
    
    def add_variable_along_latlonlayer(self, ds, var, name, units, long_name):
        """ This function self.adds missing variables to the self.FORCING class """
        ds[name] = ((northing,easting,'layer'), var)
        ds[name].attrs['units'] = units
        ds[name].attrs['long_name'] = long_name
        ds[name].encoding['_FillValue'] = -9999
        return ds
    
    def add_variable_along_layertime(self, ds, var, name, units, long_name):
        """ This function self.adds missing variables to the self.FORCING class """
        ds[name] = (('time','layer'), var)
        ds[name].attrs['units'] = units
        ds[name].attrs['long_name'] = long_name
        ds[name].encoding['_FillValue'] = -9999
        return ds
    
    def add_variable_along_layer(self, ds, var, name, units, long_name):
        """ This function self.adds missing variables to the self.FORCING class """
        ds[name] = (('layer'), var)
        ds[name].attrs['units'] = units
        ds[name].attrs['long_name'] = long_name
        ds[name].encoding['_FillValue'] = -9999
        return ds
    #==============================================================================
    # Auxiliary functions for time averaging the output 
    #==============================================================================
    def add_variable_along_layer(self, ds, var, name, units, long_name):
        """ This function self.adds missing variables to the self.FORCING class """
        ds[name] = (('layer'), var)
        ds[name].attrs['units'] = units
        ds[name].attrs['long_name'] = long_name
        ds[name].encoding['_FillValue'] = -9999
        return ds