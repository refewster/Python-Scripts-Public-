# -*- coding: cp1252 -*-
########################################################################################################
"""
Program: Monthly averages

Author: Richard Fewster
Start Date: 13/04/2020
Most Recent Update: 13/05/2020

"""
"""
STEP 1: IMPORT LIBRARIES AND DATA FILES
"""
########################################################################################################
"""
(1.1) Import required libraries.
"""
print('(1.1) Importing libraries...')

import xarray as xr
import datetime as dt
import numpy as np
import pandas as pd
import re
from netCDF4 import Dataset, date2index, num2date, date2num
import scipy
import gridfill
import iris

print('Import complete')

"""
(1.2) Create a list of required netCDF files.
"""
print('(1.2) Importing data files...')
# CMIP land-sea mask file
mask_file = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\2_CanESM5\land\sftlf_fx_CanESM5_hist-volc_r1i1p1f1_gn.nc"
# Temperature files
tas_file_hist = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\2_CanESM5\tmp\historical\tas_Amon_CanESM5_historical_r1i1p1f1_gn_185001-201412.nc"
tas_file_hist = xr.open_dataset(tas_file_hist)

tas_file_ssp1 = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\2_CanESM5\tmp\ssp1_26\tas_Amon_CanESM5_ssp126_r1i1p1f1_gn_201501-210012.nc"
tas_file_ssp1 = xr.open_dataset(tas_file_ssp1)
print('Step 1 complete')

########################################################################################################
"""
STEP 2: IMPORT CMIP CLIMATE DATA AND CALCULATE CLIMATE AVERAGES.
"""
########################################################################################################

"""
(2.1) Slice CMIP data to desired time period and study area.
"""
# Temperature slices
tas_hist_slice = tas_file_hist.sel(time=slice("1961-01-16", "1990-12-16"), lat=slice(50., 90.))
tas_ssp1_slice = tas_file_ssp1.sel(time=slice("2090-01-16", "2099-12-16"), lat=slice(50., 90.)) 

"""
(2.2) Average climate values for study time period
"""
# Time variable already returning monthly values, need to find way to average these for each month (e.g. all Januaries, Februaries, etc.)

# Temperature
tas_hist_mean = tas_hist_slice['tas'].mean('time',keep_attrs=True)

# Group by month
tas_hist_mth = tas_hist_slice['tas'].groupby("time.month").mean('time',keep_attrs=True)


"""
(2.2) Convert climate data into desired units
"""
# Convert from Kelvin to Celsius
tas_hist_mean_monthly_C = tas_hist_mth-273.15


print('Step 2 Complete')


DIR = r'G:\Climate_Data'
#tas_hist_mth.to_netcdf(DIR+'hist_mth.nc')
#tas_hist_slice.to_netcdf(DIR+'hist_time.nc')

######################################
"""
Precipitation
"""
pre_file_hist = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\pre\historical\*.nc", combine='by_coords')
pre_hist_slice = pre_file_hist.sel(time=slice("1961-01-16", "1990-12-16"), lat=slice(50., 90.)) 
# Convert from mm/second to mm per day
# 60 x 60 x 24 = 86,400 (one day)
pre_day = pre_hist_slice['pr'] * 86400.

# Calculate the number of days in each month
# adapted from http://xarray.pydata.org/en/stable/examples/monthly-means.html
dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}
# Define new function get_dpm
def get_dpm(time, calendar='standard'):
    """
    return a array of days per month corresponding to the months provided in `months`
    """
    month_length = np.zeros(len(time), dtype=np.int)
    cal_days = dpm[calendar]
    for i, (month, year) in enumerate(zip(time.month, time.year)):
        month_length[i] = cal_days[month]
    return month_length

# Use get_dpm to calculate the number of days in each month in slice
# Time coordinates match those of the precipitation data
month_length = xr.DataArray(get_dpm(pre_day.time.to_index(), calendar='noleap'),
                            coords=[pre_day.time], name='month_length')

# Multiply the days in the month by the precipitation value for that month
pre_mth= month_length * pre_day

# Average monthly precipitation values by month
pre_hist_mth = pre_mth.groupby("time.month").mean('time',keep_attrs=True)
pre_hist_mth = pre_hist_mth.rename('Average Monthly Precipitation')

# Output results
DIR = r'G:\Climate_Data'
pre_day.to_netcdf(DIR+'hist_ds.nc')
pre_mth.to_netcdf(DIR+'hist_pre_mth.nc')
pre_hist_mth.to_netcdf(DIR+'hist_pre_hist_mth.nc')



