# -*- coding: cp1252 -*-
"""
Bias correction calculations - test run for use in the downscaling program
Author: Richard Fewster
Start Date: 30/04/2020
Most Recent Update: 30/04/2020

BCor_Temp = (CMIP_fut_tmp - CMIP_hist_tmp) + CRU_tmp

BCor_Pre = (CRU_pre / CMIP_hist_pre) * CMIP_fut_pre

"""
#######################################################################################################
"""
STEP 7: BIAS CORRECTION
"""
########################################################################################################

import xarray as xr
import datetime as dt
import numpy as np
import pandas as pd
import re
from netCDF4 import Dataset, date2index, num2date, date2num

"""
(4.1) Import CRU TS 4.03 data files
"""
# Load in observational temperature dataset
CRU_tmp_dset = xr.open_mfdataset(r"G:\Climate_Data\3_Observational_data\CRU data\cru_ts4.03.1901.2018.tmp.dat.nc", combine='by_coords')
# Load in observational precipitation dataset
CRU_pre_dset = xr.open_mfdataset(r"G:\Climate_Data\3_Observational_data\CRU data\cru_ts4.03.1901.2018.pre.dat.nc", combine='by_coords')

"""
(7.x) Create dummy variables for test run
"""
#d1
CRU_tmp_d1 = CRU_tmp_dset-2 # works
#Subsetting to period and time of interest

lat =CRU_tmp_d1['lat']# unnecessary
lon =CRU_tmp_d1['lon']# unnecessary
time =CRU_tmp_d1['time'] # unnecessary

# To subset easily, use .sel and specify limits from each variable 
CRU_fut_tmp = CRU_tmp_d1.sel(lat=slice(50., 90.), lon=slice(50., 90.), time=slice("1961-01-16", "1990-12-16"))

#d2
CRU_tmp_d1 = CRU_tmp_dset-1 # works
#Subsetting to period and time of interest
# To subset easily, use .sel and specify limits from each variable 
CRU_hist_tmp = CRU_tmp_d1.sel(lat=slice(50., 90.), lon=slice(50., 90.), time=slice("1961-01-16", "1990-12-16"))

#OG
# To subset easily, use .sel and specify limits from each variable 
CRU_tmp = CRU_tmp_dset.sel(lat=slice(50., 90.), lon=slice(50., 90.), time=slice("1961-01-16", "1990-12-16"))

"""
Bias correction

BCor_Temp = (CMIP_fut_tmp - CMIP_hist_tmp) + CRU_tmp

BCor_Pre = (CRU_pre / CMIP_hist_pre) * CMIP_fut_pre
"""
BCor_Temp = (CRU_fut_tmp - CRU_hist_tmp) + CRU_tmp
xxxx

"""
(X.1) Exporting the results to netcdf format
"""

print('Data export to NetCDF...')
# Export the land only data in Kelvin
DIR = r'G:\Climate_Data\4_Python_Scripting\2_RF\4_Bias_correction_testing\Sel_CRU_tmp_d1.nc'
Sel_CRU_tmp_d1.to_netcdf(DIR)
