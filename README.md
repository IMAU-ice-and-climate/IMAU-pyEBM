# IMAU-pyEBM
The IMAU Python SEB model for snow and ice 


### Python installation
See [Anaconda Installation](https://docs.anaconda.com/anaconda/) to install Anaconda. Then, [verify your installation](https://docs.anaconda.com/anaconda/install/verify-install/).



Open a terminal (macOS) or an anaconda prompt (Windows):
* Windows: Start -> search “Anaconda Prompt”
* macOS: Launchpad -> Terminal

Create a new conda environment named "sebmodel" with python installed, and activate it:
```bash
conda create -n sebmodel python
conda activate sebmodel
```
Note that Python version 3.9.2 might be required for compatilibity with the numba package (on March 2022).

Install the required Python packages:

```bash
conda install pandas distributed xarray numba bottleneck netcdf4 dask-distributed
```
This can take a few minutes.

### Download the SEB model
In the terminal, navigate to a directory _INSTALLDIR_ of your choice.
```bash
cd INSTALLDIR
```

* For developpers (using git), clone the repository.
  In the terminal, run:
  ```bash
  git clone https://github.com/mvantiggelen/seb-model_python.git
  ```

* For users, download the [.zip file](https://github.com/mvantiggelen/seb-model_python/archive/refs/heads/main.zip) and unzip in _INSTALLDIR_:
  ```bash
  wget "https://github.com/IMAU-ice-and-climate/IMAU-pyEBM/archive/refs/heads/main.zip"
  ```

### Preparing forcing and parameter files
Make sure that you have the following directory structure :
```
seb-model_python-main
├── input
│   ├── ant_aws14
│   │   ├── ant_aws14_ebm.py
│   │   ├── ant_aws14_HOUR-EBM.txt
│   ├── grl_2Dtest
│   │   ├── grl_2Dtest.py
│   │   ├── grl_2Dtest_HOUR-EBM.nc
├── output
│   ├── ant_aws14
│   ├── grl_2Dtest
├── sebmodel
```
In the 'input' directory, there should be one folder per station. 
In each station folder, there should be 2 files:

* awsid_ebm.py: parameters

* awsid_HOUR-EBM.txt (1D) or awsid_HOUR-EBM.nc (2D): atmospheric forcing 

Each of these files needs to be in a very specific format for the SEB model to run.

### Running the SEB model from the terminal
Go to the sebmodel directory

```bash
cd sebmodel/
```

Run the preparation script:
```bash
python sebmodel_pre.py --aws ant_aws07 --v 1D 
```
You can change the name of the station (aws_ant07) and the type of version (1D or 2D).
This script doe snothing more than copy the parameter file awsid_ebm.py from the input folder to the info.py file in the source code directory. 


Now, run the SEB model:
```bash
python sebmodel.py
```

If everything works, you should see the something like:


```bash
(sebmodel) UU-CQ0WFQCMYW:sebmodel Tigge006$ python sebmodel_pre.py --aws ant_aws14 --v 1D
chstation is  ant_aws14
version is  1D
convert info file is  False
(sebmodel) UU-CQ0WFQCMYW:sebmodel Tigge006$ python sebmodel.py 
AWS is =  ant_aws14
/Users/Tigge006/surfdrive/04_Scripts/SEB_model/IMAU-pyEBM/sebmodel

 Maximum available time interval from 2009-01-21T22:00 until 2022-11-27T23:00. Time steps: 121394 


--------------------------------------------------------------
Checking input FORCING .... 

Temperature FORCING (T) present
Specific humidity FORCING (q) present
Incoming hortwave FORCING (Sin) present
Outgoing shortwave FORCING (Sout) present
Wind velocity FORCING (WS) present
Precipitation FORCING (precip) present
WARNING Precip in input out of range MAX: 74.64 MIN: 0.00 

Pressure FORCING (P) present
Albedo FORCING (alb) present
Incoming longwave FORCING (Lin) present
Outgoing longwave FORCING (Lout) present
Serie FORCING (Serie) present
Aerodynamic roughness FORCING (z0m) present
Height of T-q sensor FORCING (zt) present
Height of WS sensor FORCING (zm) present


Output dataset ... ok
Starting initgrid
initial number of layers is:  29 29 1.6104124773293442 8.84747999408851 1
Starting initsnow
Starting initgrains
Starting get_errorflag
Starting checkdata
END data checking, start EBM calculation
1 resizegrid 28 0 -4.630686355675629 -2.3153431778378146 -1490.3394073517702 321.83985113246257
REDEFGRID: added layer at bottom 34 2.0 7.7508149209270645 5.501610087699024 4.749199894689264 0.5016199647765767 2 576.0305129737757
0 % done
2 % done
.
.
.
97 % done
99 % done
	 Time required to do calculations:    0 minutes 10.6027 seconds 



--------------------------------------------------------------
Write results ...
-------------------------------------------------------------- 

	 Time required tor write restart and output files:    0 minutes 10.0734 seconds 

	 Total run duration:    0 minutes 20.6799 seconds 

--------------------------------------------------------------
	 SIMULATION WAS SUCCESSFUL
--------------------------------------------------------------

```

The output files are written in the 'output' folder:

```
seb-model_python-main
├── input
│   ├── ant_aws14
│   │   ├── ant_aws14_ebm.py
│   │   ├── ant_aws14_HOUR-EBM.txt
├── output
│   ├── ant_aws14
│   │   ├── ant_aws14_ALL.nc
│   │   ├── ant_aws14_DAY.nc
│   │   ├── ant_aws14_MONTH.nc
│   │   ├── ant_aws14_SEAS.nc
│   │   ├── ant_aws14_YEAR.nc
├── sebmodel
```
The raw (hourly or 30min) output is stored in the 'ant_aws14_ALL.nc' file. 

The other files contain the same variables but averaged per day, month, seasonnaly or yearly. 

### Opening the output data

The output file are in netCDF format. They can be opened using the 'xarray' module in python.

For instance, make a new .py file which contains something like this:

```python 

import xarray as xr

file = 'IMAU-pyEBM/output/ant_aws14/ant_aws14_ALL.nc' # change this !
ds = xr.open_dataset(file)
ds['melt'].plot()

```

### Additionnal information

For similar SEB models, please check:

* COSIPY: https://github.com/cryotools/cosipy
* CryoGrid Community model: https://github.com/CryoGrid/CryoGridCommunity_source
* EB_AUTO Spreadsheet Energy Balance Model: https://github.com/atedstone/ebmodel
* Distributed Energy Balance Model : https://github.com/regine/meltmodel
* GEUS Surface Energy Balance and Firn Model: https://github.com/BaptisteVandecrux/SEB_Firn_model

If you have any questions, feel free to contact Maurice van Tiggelen, m.vantiggelen@uu.nl 