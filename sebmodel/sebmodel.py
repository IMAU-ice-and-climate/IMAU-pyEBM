#!/usr/bin/env python3

''' Third-party import '''
import pandas as pd
from datetime import datetime
import logging
import os
from itertools import product
from dask.diagnostics import ProgressBar
from dask.distributed import progress, wait, as_completed,  Client, LocalCluster
import tracemalloc
import linecache
import yaml
from dask import compute, delayed
from tornado import gen
from dask_jobqueue import SLURMCluster

''' Local imports '''
from info import * 
# from input import * 
from globals import *
from IO import *
from sebmodel_core import * 

def main():
    """ SEB runscript, based on COSIPY runscript https://cryo-tools.org/tools/cosipy/
    Params 
    ======
    chstation: name of the station (command line -aws)
    version: 1D or 2D (command line -v)
    convert: switch to convert info file from fortran SEB model (command line -c)
    Returns
    ======
    """
       
    print('AWS is = ',chstation)

    output_netcdf = chstation + '_ALL.nc'
    data_path = '../output/' + chstation + '/'
    
    #------------------------------------------
    # Create input and output dataset
    #------------------------------------------ 
    print(os.getcwd())
    IO = IOClass()
    FORCING = IO.create_forcing_file() 
    FORCING.head()
    
    # Create global output dataset
    OUTPUT = IO.create_output_file() 

    #----------------------------------------------
    # Calculation - Multithreading using all cores  
    #----------------------------------------------
    
    # Auxiliary variables for futures
    futures = []

    # Measure time
    start_time = datetime.now()
    
    #-----------------------------------------------
    # Create a client for distributed calculations
    #-----------------------------------------------
    if version == '2D':
        if (slurm_use):
            with SLURMCluster(job_name=name, cores=cores, processes=cores, memory=memory, 
                              shebang=shebang, job_extra=slurm_parameters, 
                              local_directory=slurm_local_directory) as cluster:
                cluster.scale(nodes*cores)   
                print(cluster.job_script())
                print("You are using SLURM!\n")
                print(cluster)
                run_sebmodel_2D(cluster, IO, FORCING, OUTPUT, futures)
        else:
            with LocalCluster(scheduler_port=8786, n_workers=None, threads_per_worker=1, silence_logs=True) as cluster:
                print(cluster)
                run_sebmodel_2D(cluster, IO, FORCING, OUTPUT, futures)
            
    else:  
        run_sebmodel_1D(IO, FORCING, OUTPUT)

    print('\n')
    print('--------------------------------------------------------------')
    print('Write results ...')
    print('-------------------------------------------------------------- \n')
    start_writing = datetime.now()

    #-----------------------------------------------
    # Write results and restart files
    #-----------------------------------------------  
    encoding = dict()
    for var in IO.get_output().data_vars:
        dataMin = IO.get_output()[var].min(skipna=True).values
        dataMax = IO.get_output()[var].max(skipna=True).values
        dtype = 'int16'
        FillValue = -9999 
        encoding[var] = dict(zlib=True, complevel=compression_level,_FillValue=FillValue)
  
    IO.average_output(OUTPUT)   
    IO.get_output().to_netcdf(os.path.join(data_path,output_netcdf), encoding=encoding, mode = 'w')
    IO.get_output_daily().to_netcdf(os.path.join(data_path,chstation + '_DAY.nc'), mode = 'w')
    IO.get_output_monthly().to_netcdf(os.path.join(data_path,chstation + '_MONTH.nc'), mode = 'w')
    IO.get_output_yearly().to_netcdf(os.path.join(data_path,chstation + '_YEAR.nc'),  mode = 'w')
    IO.get_output_seas().to_netcdf(os.path.join(data_path,chstation + '_SEAS.nc'),  mode = 'w')
    
    #-----------------------------------------------
    # Stop time measurement
    #-----------------------------------------------
    duration_run = datetime.now() - start_time
    duration_run_writing = datetime.now() - start_writing

    #-----------------------------------------------
    # Print out some information
    #-----------------------------------------------
    print("\t Time required tor write restart and output files: %4g minutes %2g seconds \n" % (duration_run_writing.total_seconds()//60.0,duration_run_writing.total_seconds()%60.0))
    print("\t Total run duration: %4g minutes %2g seconds \n" % (duration_run.total_seconds()//60.0,duration_run.total_seconds()%60.0))
    print('--------------------------------------------------------------')
    print('\t SIMULATION WAS SUCCESSFUL')
    print('--------------------------------------------------------------')         
       
def run_sebmodel_2D(cluster, IO, FORCING, OUTPUT, futures):

    with Client(cluster) as client:
        print('--------------------------------------------------------------')
        print('\t Starting clients and submit jobs ... \n')
        print('-------------------------------------------------------------- \n')

        print(cluster)
        print(client)

        # Get dimensions of the whole domain
        ny = FORCING.dims['y']
        nx = FORCING.dims['x']

        # Get some information about the cluster/nodes
        total_grid_points = FORCING.dims['y']*FORCING.dims['x']
        total_cores = cores*nodes
        points_per_core = total_grid_points // total_cores
        print(total_grid_points, total_cores, points_per_core)

        # Distribute data and model to workers
        start_res = datetime.now()

        for y,x in product(range(FORCING.dims['y']),range(FORCING.dims['x'])):
            
            mask = FORCING.mask.isel(y=y, x=x)
            # Only run SEB model on glacier points
            if (mask.values==1):
                futures.append(client.submit(sebmodel_core, FORCING.isel(y=y, x=x), y, x))
                    
        # Finally, do the calculations and print the progress
        progress(futures)

        # Create numpy arrays which aggregates all local outputs
        IO.create_global_output_arrays()

        #---------------------------------------
        # Assign local outputs to global 
        #---------------------------------------
        for future in as_completed(futures):

            # Get the outputs from the workers
            # a = future.result()
            (indY,indX,Sin,Sout,Lin,Loutobs,Loutmod,
            SH,LE,GH,Restsource,Source,sumdivs,T,P,WS,q,T0,q0, 
            T2m,q2m,WS10m,z0m,z,dz,temp,dens,
            kice,cpice,rhocp,energy, lid,
            mass, grainsize, water, ice, dsdz, refrfrac,
            icemelt, icemeltmdt, dsnowacc, hsnowmod, runoff, runoffdt, surfwater,
            melt, meltdt, surfmelt, surfmeltdt, sumdrift, subl, subldt ,precip,
            precipdt,dens_lay1,temp_lay1,dz_lay1,sumwater,topwater,air_content,effective_air_content,errorflag) = future.result()
            
            IO.copy_local_to_global(
                indY,indX,Sin,Sout,Lin,Loutobs,Loutmod,
                SH,LE,GH,Restsource,Source,sumdivs,T,P,WS,q,T0,q0,
                T2m,q2m,WS10m,z0m,z,dz,temp,dens,kice,cpice,rhocp,energy,
                lid, mass, grainsize, water, ice, dsdz,refrfrac,
                icemelt, icemeltmdt, dsnowacc, hsnowmod, runoff, runoffdt, surfwater,
                melt, meltdt, surfmelt, surfmeltdt, sumdrift, subl, subldt,
                precip,precipdt,dens_lay1,temp_lay1,dz_lay1,sumwater,topwater,air_content,effective_air_content,errorflag
                )

            # IO.copy_local_restart_to_global(indY,indX,local_restart)

            # Write outputs to file
            IO.write_outputs_to_file()

        # Measure time
        end_res = datetime.now()-start_res 
        print("\t Time required to do calculations: %4g minutes %2g seconds \n" % 
              (end_res.total_seconds()//60.0,end_res.total_seconds()%60.0))
        
def run_sebmodel_1D(IO, FORCING, OUTPUT):
    start_res = datetime.now()
    IO.create_global_output_arrays()

    (indY,indX,Sin,Sout,Lin,Loutobs,Loutmod,
    SH,LE,GH,Restsource,Source,sumdivs,T,P,WS,q,T0,q0, 
    T2m,q2m,WS10m,z0m,z,dz,temp,dens,
    kice,cpice,rhocp,energy, lid,
    mass, grainsize, water, ice, dsdz, refrfrac,
    icemelt, icemeltmdt, dsnowacc, hsnowmod, runoff, runoffdt, surfwater,
    melt, meltdt, surfmelt, surfmeltdt, sumdrift, subl, subldt, precip,
    precipdt, dens_lay1, temp_lay1, dz_lay1, sumwater, topwater,air_content,effective_air_content, errorflag) = \
        sebmodel_core(FORCING.isel(lat=0, lon=0), 0, 0) 

    IO.copy_local_to_global(indY,indX,Sin,Sout,Lin,Loutobs,Loutmod,
        SH,LE,GH,Restsource,Source,sumdivs,T,P,WS,q,T0,q0,
        T2m,q2m,WS10m,z0m,z,dz,temp,dens,kice,cpice,rhocp,energy,
        lid, mass, grainsize, water, ice, dsdz,refrfrac,
        icemelt, icemeltmdt, dsnowacc, hsnowmod, runoff, runoffdt, surfwater,
        melt, meltdt, surfmelt, surfmeltdt, sumdrift, subl, subldt,
        precip, precipdt, dens_lay1, temp_lay1, dz_lay1, sumwater, topwater,air_content,effective_air_content, errorflag)      

    IO.write_outputs_to_file()

    end_res = datetime.now()-start_res 
    print("\t Time required to do calculations: %4g minutes %2g seconds \n" % 
          (end_res.total_seconds()//60.0,end_res.total_seconds()%60.0))


''' MODEL EXECUTION '''

if __name__ == "__main__":
    
   main()
