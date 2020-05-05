# -*- coding: cp1252 -*-
########################################################################################################
"""
Program: Downscaling Program (DRAFT)

Author: Richard Fewster
Start Date: 30/04/2020
Most Recent Update: 05/05/2020

Program Purpose: To downscale and bias-correct CMIP6 climate outputs to the resolution of a higher-resolution observational dataset.
Specifically:
i) Regrid coarse-scale CMIP6 output files onto the higher-resolution grid of the CRU TS 4.03 dataset.
ii) Extrapolate terrrestrial climate over the ocean, to generate terrrestrial seasonality along coastlines.
iii) Bias-correct future climate outputs by applying the future climate anomaly to the observational data.
iv) Finally, use the observational dataset's land-sea mask to remove any grid cells located in oceanic regions.

"""
########################################################################################################

print('Downscaling Program.')
print('Start.')

########################################################################################################
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

"""
(1.2) Create a list of required netCDF files.
"""
print('(1.2) Importing data files...')
# CMIP land-sea mask file
mask_file = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\land\sftlf_fx_NorESM2-MM_historical_r1i1p1f1_gn.nc"

# to combine all netcdf files together use xr.open_mfdataset(path, combine = 'by coords') <- reorders the arrays before concatenating.
# use * at end of path name to merge all files ending in .nc in the specifed folder

# Temperature files
tas_file_hist = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\tmp\historical\*.nc", combine='by_coords')
tas_file_ssp1 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\tmp\ssp1_26\*.nc", combine='by_coords')
tas_file_ssp2 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\tmp\ssp2_45\*.nc", combine='by_coords')
tas_file_ssp3 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\tmp\ssp3_70\*.nc", combine='by_coords')
tas_file_ssp5 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\tmp\ssp5_85\*.nc", combine='by_coords')

# Precipitation files
pre_file_hist = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\pre\historical\*.nc", combine='by_coords')
pre_file_ssp1 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\pre\ssp1_26\*.nc", combine='by_coords')
pre_file_ssp2 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\pre\ssp2_45\*.nc", combine='by_coords')
pre_file_ssp3 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\pre\ssp3_70\*.nc", combine='by_coords')
pre_file_ssp5 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\pre\ssp5_85\*.nc", combine='by_coords')

########################################################################################################
"""
STEP 2: IMPORT CMIP CLIMATE DATA AND CALCULATE CLIMATE AVERAGES.
"""
########################################################################################################

"""
(2.1) Slice CMIP data to desired time period and study area.
"""
print('(2.1) Slicing CMIP files...')
# Create a time slice - Set dates of interest
# .sel(time = slice( start date, end date)

# Temperature slices
tas_ssp5_slice = tas_file_hist.sel(time=slice("2090-01-16", "2099-12-16"), lat=slice(50., 90.)) 
tas_ssp5_slice = tas_file_ssp1.sel(time=slice("2090-01-16", "2099-12-16"), lat=slice(50., 90.)) 
tas_ssp5_slice = tas_file_ssp2.sel(time=slice("2090-01-16", "2099-12-16"), lat=slice(50., 90.)) 
tas_ssp5_slice = tas_file_ssp3.sel(time=slice("2090-01-16", "2099-12-16"), lat=slice(50., 90.)) 
tas_ssp5_slice = tas_file_ssp5.sel(time=slice("2090-01-16", "2099-12-16"), lat=slice(50., 90.)) 

# Precipitation files
pre_ssp5_slice = pre_file_hist.sel(time=slice("2090-01-16", "2099-12-16"), lat=slice(50., 90.)) 
pre_ssp5_slice = pre_file_ssp1.sel(time=slice("2090-01-16", "2099-12-16"), lat=slice(50., 90.)) 
pre_ssp5_slice = pre_file_ssp2.sel(time=slice("2090-01-16", "2099-12-16"), lat=slice(50., 90.)) 
pre_ssp5_slice = pre_file_ssp3.sel(time=slice("2090-01-16", "2099-12-16"), lat=slice(50., 90.)) 
pre_ssp5_slice = pre_file_ssp5.sel(time=slice("2090-01-16", "2099-12-16"), lat=slice(50., 90.)) 

"""
(2.2) Average climate values for each month
"""
# Time variable already returning monthly values, need to find way to average these for each month (e.g. all Januarys, Februarys, etc.)

# Temperature
#tas_hist_mean_monthly_C =
#tas_ssp1_mean_monthly_C =
#tas_ssp2_mean_monthly_C =
#tas_ssp3_mean_monthly_C =
#tas_ssp5_mean_monthly_C =

#Precipitation
#pre_hist_mean_monthly =
#pre_ssp1_mean_monthly =
#pre_ssp2_mean_monthly =
#pre_ssp3_mean_monthly =
#pre_ssp5_mean_monthly =

"""
(2.2) Convert climate data into desired units
"""
# Convert from Kelvin to Celsius
# e.g. tas_hist_mean_monthly_C = tas_hist_mean_monthly_K-273.15
# e.g. tas_ssp1_mean_monthly_C = tas_ssp1_mean_monthly_K-273.15
# e.g. tas_ssp2_mean_monthly_C = tas_ssp2_mean_monthly_K-273.15
# e.g. tas_ssp3_mean_monthly_C = tas_ssp3_mean_monthly_K-273.15
# e.g. tas_ssp5_mean_monthly_C = tas_ssp5_mean_monthly_K-273.15

# Convert from mm/second to mm per month
# 60 x 60 x 24 = 86,400 (one day)
# 86,400 x 365 = 31,556,926 (one year)
# 31,556,926 / 12 = 2,629,743.8 (estimate for one month)
#pre_hist_mean_monthly_mm = pre_hist_mean_monthly * 2,629,743.8 
#pre_ssp1_mean_monthly_mm = pre_ssp1_mean_monthly * 2,629,743.8 
#pre_ssp2_mean_monthly_mm = pre_ssp2_mean_monthly * 2,629,743.8 
#pre_ssp3_mean_monthly_mm = pre_ssp3_mean_monthly * 2,629,743.8 
#pre_ssp5_mean_monthly_mm = pre_ssp5_mean_monthly * 2,629,743.8 

# Better way of doing this?? ^^

########################################################################################################
"""
STEP 3: MASK CMIP OCEANIC CELLS 
"""
########################################################################################################

"""
(3.1) Use xarray to assign the land cover percentage data for the CMIP model to a new object.
"""
print('(3.1) Assign CMIP land-sea mask...')
# sftlf is the standardised variable name for land percentage cover
mask_dset = xr.open_dataset(mask_file) #Use xarray to open the mask dataset
land_perc = mask_dset['sftlf'] # assign the land percentage variable to a new object


print('Max land area (CMIP mask):', land_perc.data.max(), '%') # check that max land area is 100 %
print('Min land area (CMIP mask):', land_perc.data.min(), '%') # check that min land area is 0 %

"""
(3.2) Mask out ocean in CMIP datasets (i.e. selecting only grid cells with > 50 % land)
"""
print('(3.2) Apply land-sea mask...')
#numpy includes a np.where function that allows us to simply use a logical command

# Mask out temperature data
tas_hist_land_C = tas_hist_mean_monthly_C.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
tas_ssp1_land_C = tas_ssp1_mean_monthly_C.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
tas_ssp2_land_C = tas_ssp2_mean_monthly_C.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
tas_ssp3_land_C = tas_ssp3_mean_monthly_C.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
tas_ssp5_land_C = tas_ssp5_mean_monthly_C.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %

# Mask out preciptiation data
pre_hist_land_C = pre_hist_mean_monthly_mm.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
pre_ssp1_land_C = pre_ssp1_mean_monthly_mm.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
pre_ssp2_land_C = pre_ssp2_mean_monthly_mm.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
pre_ssp3_land_C = pre_ssp3_mean_monthly_mm.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
pre_ssp5_land_C = pre_ssp5_mean_monthly_mm.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %


########################################################################################################
"""
STEP 4: IMPORT MODERN OBSERVATIONAL DATASET (CRU TS 4.03)
"""
########################################################################################################

"""
(4.1) Import CRU TS 4.03 data files
"""
print('(4.1) Import CRU TS 4.03 datasets...')

# Load in observational temperature dataset
CRU_tmp_dset = xr.open_mfdataset(r"G:\Climate_Data\3_Observational_data\CRU data\cru_ts4.03.1901.2018.tmp.dat.nc", combine='by_coords')
# Load in observational precipitation dataset
CRU_pre_dset = xr.open_mfdataset(r"G:\Climate_Data\3_Observational_data\CRU data\cru_ts4.03.1901.2018.pre.dat.nc", combine='by_coords')

"""
(4.2) Slicing observational data.
"""
# REMEMBER: CRU datasets use -179.75 -> 179.75 for lon
CRU_tmp_slice = CRU_tmp_dset.sel(time=slice("1961-01-16", "1990-12-16"), lat=slice(50., 90.)) # Slice to match study region of SSP files
CRU_pre_slice = CRU_pre_dset.sel(time=slice("1961-01-16", "1990-12-16"), lat=slice(50., 90.)) # Slice to match study region of SSP files


########################################################################################################
"""
STEP 5: REGRIDDING TO HIGHER RESOLUTION
"""
########################################################################################################








########################################################################################################
"""
STEP 6: EXTRAPOLATION 
"""
########################################################################################################









########################################################################################################
"""
STEP 7: BIAS CORRECTION
"""
########################################################################################################
"""
(7.3) Variable formatting
"""
#CMIP_fut_tmp = 
#CMIP_hist_tmp
#CRU_tmp =

#CMIP_fut_pre = 
#CMIP_hist_pre =
#CRU_pre = 


"""
(7.2) Temperature Bias correction

BCor_Temp = (CMIP_fut_tmp - CMIP_hist_tmp) + CRU_tmp
"""
# Replace CRU dummy variables with CMIP equivalents
SSP1_BCor_Temp = (SSP1_fut_tmp - CMIP_hist_tmp) + CRU_tmp
SSP2_BCor_Temp = (SSP2_fut_tmp - CMIP_hist_tmp) + CRU_tmp
SSP3_BCor_Temp = (SSP3_fut_tmp - CMIP_hist_tmp) + CRU_tmp
SSP5_BCor_Temp = (SSP5_fut_tmp - CMIP_hist_tmp) + CRU_tmp

"""
(7.3) Precipitation Bias correction

BCor_Pre = (CRU_pre / CMIP_hist_pre) * CMIP_fut_pre
"""
# Calculate alpha...
a = CRU_pre / CMIP_hist_pre
# Set limits for alpha...
a = xr.where(a < 0.25, 0.25, a)
a_ltd = xr.where(a > 4.0, 4.0, a)

# Apply the limited alpha coefficient to bias correct future precipitation
SSP1_BCor_Pre = a_ltd * SSP1_fut_pre
SSP2_BCor_Pre = a_ltd * SSP2_fut_pre
SSP3_BCor_Pre = a_ltd * SSP3_fut_pre
SSP5_BCor_Pre = a_ltd * SSP5_fut_pre

########################################################################################################
"""
STEP 8: APPLY OBSERVATIONAL LAND-SEA MASK
"""
########################################################################################################

# (8.1) Create CRU land-sea mask


# (8.2) Apply CRU land sea mask to downscaled CMIP data

# Mask temperature files...
# SSP1_BCor_Temp_land = 
# SSP2_BCor_Temp_land = 
# SSP3_BCor_Temp_land = 
# SSP5_BCor_Temp_land = 

# Mask precipitation files...
# SSP1_BCor_Pre_land =
# SSP2_BCor_Pre_land = 
# SSP3_BCor_Pre_land = 
# SSP5_BCor_Pre_land =

########################################################################################################
"""
STEP XX: OUTPUT RESULTS 
"""
########################################################################################################

"""
(X.1) Exporting the results to netcdf format
"""
print('Data export to NetCDF...')
# Temperature files
DIR = r'G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\downscaled_outputs\NorESM2_MM_downscaled_tas_'
SSP1_BCor_Temp_land.to_netcdf(DIR+'SSP1.nc')
SSP2_BCor_Temp_land.to_netcdf(DIR+'SSP2.nc')
SSP3_BCor_Temp_land.to_netcdf(DIR+'SSP3.nc')
SSP5_BCor_Temp_land.to_netcdf(DIR+'SSP5.nc')

# Precipitation files
DIR = r'G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\downscaled_outputs\NorESM2_MM_downscaled_pre_'
SSP1_BCor_Pre_land.to_netcdf(DIR+'SSP1.nc')
SSP2_BCor_Pre_land.to_netcdf(DIR+'SSP2.nc')
SSP3_BCor_Pre_land.to_netcdf(DIR+'SSP3.nc')
SSP5_BCor_Pre_land.to_netcdf(DIR+'SSP5.nc')

"""
(X.2) Output monthly climate values for each 0.5 degree grid cell as .csv
"""
#print('Data export to .csv...')








#########################################
print('End.')
