# -*- coding: cp1252 -*-
"""
Bias correction calculations - test run for use in the downscaling program
Author: Richard Fewster
Start Date: 30/04/2020
Most Recent Update: 04/05/2020

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
(1) TEMPERATURE
(1.1) Import CRU TS 4.03 tmp data files
"""
# Load in observational temperature dataset
CRU_tmp_dset = xr.open_mfdataset(r"G:\Climate_Data\3_Observational_data\CRU data\cru_ts4.03.1901.2018.tmp.dat.nc", combine='by_coords')

"""
(1.2) Create dummy variables for test run
"""
#OG
# To subset easily, use .sel and specify limits from each variable 
CRU_tmp = CRU_tmp_dset.sel(lat=slice(50., 90.), lon=slice(50., 90.), time=slice("1961-01-16", "1990-12-16"))

#dummy 1 - historical
CRU_tmp_d1 = CRU_tmp_dset-1 
#Subsetting to period and time of interest
CRU_hist_tmp = CRU_tmp_d1.sel(lat=slice(50., 90.), lon=slice(50., 90.), time=slice("1961-01-16", "1990-12-16"))

#dummy 2 - future
CRU_tmp_d2 = CRU_tmp_dset-2
#Subsetting to period and time of interest
CRU_fut_tmp = CRU_tmp_d2.sel(lat=slice(50., 90.), lon=slice(50., 90.), time=slice("1961-01-16", "1990-12-16"))

"""
(1.3) Temperature Bias correction

BCor_Temp = (CMIP_fut_tmp - CMIP_hist_tmp) + CRU_tmp
"""
# Replace CRU dummy variables with CMIP equivalents
BCor_Temp = (CRU_fut_tmp - CRU_hist_tmp) + CRU_tmp

"""
(1.4) Exporting the tmp results to netcdf format
"""
print('Temperature data export to NetCDF...')
# Export the original, sliced CRU data
DIR = r'G:\Climate_Data\4_Python_Scripting\2_RF\4_Bias_correction_testing\CRU_tmp.nc'
CRU_tmp.to_netcdf(DIR)

DIR = r'G:\Climate_Data\4_Python_Scripting\2_RF\4_Bias_correction_testing\BCor_Temp.nc'
BCor_Temp.to_netcdf(DIR)

###################################################################
###################################################################
###################################################################

"""
(2) PRECIPITATION
(2.1) Import CRU TS 4.03 pre data files
"""
# Load in observational precipitation dataset
CRU_pre_dset = xr.open_mfdataset(r"G:\Climate_Data\3_Observational_data\CRU data\cru_ts4.03.1901.2018.pre.dat.nc", combine='by_coords')

"""
(2.2) Create dummy variables for test run
"""
#OG
# To subset easily, use .sel and specify limits from each variable 
CRU_pre = CRU_pre_dset.sel(lat=slice(50., 90.), lon=slice(50., 90.), time=slice("1961-01-16", "1990-12-16"))

#dummy 1 - historical
CRU_pre_d1 = CRU_pre_dset+1
#Subsetting to period and time of interest
CRU_hist_pre = CRU_pre_d1.sel(lat=slice(50., 90.), lon=slice(50., 90.), time=slice("1961-01-16", "1990-12-16"))

#dummy 2 - future
CRU_pre_d2 = CRU_pre_dset+2
#Subsetting to period and time of interest
CRU_fut_pre = CRU_pre_d2.sel(lat=slice(50., 90.), lon=slice(50., 90.), time=slice("1961-01-16", "1990-12-16"))

"""
(1.3) Precipitation Bias correction

BCor_Pre = (CRU_pre / CMIP_hist_pre) * CMIP_fut_pre
"""
# Replace CRU dummy variables with CMIP equivalents
alpha = CRU_pre / CRU_hist_pre
# Set limits for alpha
alpha = xr.where(a < 0.25, 0.25, a)
alpha_ltd = xr.where(alpha > 4.0, 4.0, a)
# Apply the limited alpha coefficient to bias correct future precipitation
BCor_Pre = alpha_ltd * CRU_fut_pre

"""
(1.4) Exporting the pre results to netcdf format
"""
print('Precipitation data export to NetCDF...')
# Export the original, sliced CRU data
DIR = r'G:\Climate_Data\4_Python_Scripting\2_RF\4_Bias_correction_testing\CRU_pre.nc'
CRU_pre.to_netcdf(DIR)

DIR = r'G:\Climate_Data\4_Python_Scripting\2_RF\4_Bias_correction_testing\BCor_pre.nc'
BCor_Pre.to_netcdf(DIR)
