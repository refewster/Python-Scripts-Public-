# -*- coding: cp1252 -*-
"""
Program Purpose: To mask out ocean grid cells in surface temperature netcdfs from the NorESM2-MM netcdf.
Author: Richard Fewster
Start Date: 28/04/2020
Most Recent Update: 28/04/2020
"""
print('Start.')
"""
1) Import required libraries.
"""
import xarray as xr
import datetime as dt
import numpy as np
import pandas as pd
import re
from netCDF4 import Dataset, date2index, num2date, date2num

# Make sure files are odered correctly in folders before proceeding.

"""
2) Create a list of required netCDF files.
"""
mask_file = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\land\sftlf_fx_NorESM2-MM_historical_r1i1p1f1_gn.nc"
climate_file_hist = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\pre\historical\*.nc", combine='by_coords')
# to combine all netcdf files together use xr.open_mfdataset(path, combine = 'by coords') <- reorders the arrays before concatenating.
# use * at end of path name to merge all files ending in .nc in the specifed folder
# need dask library installed
climate_file_ssp1 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\pre\ssp1_26\*.nc", combine='by_coords')
climate_file_ssp2 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\pre\ssp2_45\*.nc", combine='by_coords')
climate_file_ssp3 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\pre\ssp3_70\*.nc", combine='by_coords')
climate_file_ssp5 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\pre\ssp5_85\*.nc", combine='by_coords')

"""
3.1) Firstly, we will use xarray to assign the land cover percentage data to a new object.
"""
mask_dset = xr.open_dataset(mask_file) #Use xarray to open the mask dataset
land_perc = mask_dset['sftlf'] # assign the land percentage variable to a new object
# sftlf is the standardised term for land percentage cover
print('Land_perc:', land_perc)

print('Max land area:', land_perc.data.max(), '%') # check that max land area is 100 %
print('Min land area:', land_perc.data.min(), '%') # check that min land area is 0 %

"""
3.2) Second we will load in the climate data and assign precipitation (pr) to a new object.
New addition - slice to a specific time period (Check that this is producing correct numbers). 
"""
# open climate file as a new dset
climate_dset = climate_file_ssp1
# assign time to a new variable
time = climate_dset['time']
# Create a time slice - Set dates of interest
# .sel(time = slice( start date, end date)
time_slice = climate_dset.sel(time=slice("2090-01-16", "2099-12-16"))
# Average all tas values in that time period (object = time_slice)
clim_slice = time_slice['pr'].mean('time',keep_attrs=True)
# View some results
clim_slice[:3]

#convert to mm per year:
pr_yr = clim_slice*31556926
# convert from seconds into years
# 60 x 60 x 24 = 86,400 (one day)
# 86,400 x 365 = 31556926 (one year)
pr_month = pr_yr/12
# Divide by 12 to get estimate for per month

"""
3.3) Masking process
Masking out ocean (i.e. selecting only grid cells with > 50 % land)
"""
clim_land = pr_yr.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
#numpy added a np.where function that allows us to simply use a logical command
print(clim_land)

"""
4) Exporting the results to netcdf format
"""
# Export the land only data
land_dataDIR_pr = r'G:\Climate_Data\4_Python_Scripting\2_RF\3_NorESM2_MM\Applying_land_sea_mask\clim_pr_ssp1_2090_2099_land.nc'
clim_land.to_netcdf(land_dataDIR_pr)

#########################################
print('End.')
