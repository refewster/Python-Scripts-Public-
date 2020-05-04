# -*- coding: cp1252 -*-
########################################################################################################
"""
Program: Downscaling Program (DRAFT)

Author: Richard Fewster
Start Date: 30/04/2020
Most Recent Update: 30/04/2020

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

mask_file = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\land\sftlf_fx_NorESM2-MM_historical_r1i1p1f1_gn.nc"
climate_file_hist = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\tmp\historical\*.nc", combine='by_coords')
# to combine all netcdf files together use xr.open_mfdataset(path, combine = 'by coords') <- reorders the arrays before concatenating.
# use * at end of path name to merge all files ending in .nc in the specifed folder
climate_file_ssp1 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\tmp\ssp1_26\*.nc", combine='by_coords')
climate_file_ssp2 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\tmp\ssp2_45\*.nc", combine='by_coords')
climate_file_ssp3 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\tmp\ssp3_70\*.nc", combine='by_coords')
climate_file_ssp5 = xr.open_mfdataset(r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\tmp\ssp5_85\*.nc", combine='by_coords')

########################################################################################################
"""
STEP 2: IMPORT CMIP CLIMATE DATA AND CALCULATE CLIMATE AVERAGES.
"""
########################################################################################################

"""
(2.1) Load in the CMIP climate data and assign the time variabe to a new object.
"""
print('(2) Processing CMIP files...')

# open climate file as a new dset
CMIP_dset = climate_file_ssp5
# assign time to a new variable
time =CMIP_dset['time']

"""
(2.2) Slice to a specific time period. 
"""
# Create a time slice - Set dates of interest
# .sel(time = slice( start date, end date)
time_slice = CMIP_dset.sel(time=slice("2090-01-16", "2099-12-16"))
# Average all tas values in that time period (object = time_slice)
# Produces a slice of mean surface temperature in kelvin
CMIP_slice_K = time_slice['tas'].mean('time',keep_attrs=True)
# View some results
#clim_slice_K[:3]

"""
(2.3) Convert temperature data from Kelvin to Celsius
"""
# Convert from Kelvin to Celsius
CMIP_slice_C = CMIP_slice_K-273.15

########################################################################################################
"""
STEP 3: MASK CMIP OCEANIC CELLS 
"""
########################################################################################################

"""
(3.1) Use xarray to assign the land cover percentage data for the CMIP model to a new object.
"""
print('(3.1) Assign CMIP land-sea mask...')

mask_dset = xr.open_dataset(mask_file) #Use xarray to open the mask dataset
land_perc = mask_dset['sftlf'] # assign the land percentage variable to a new object
# sftlf is the standardised term for land percentage cover

print('Max land area (CMIP mask):', land_perc.data.max(), '%') # check that max land area is 100 %
print('Min land area (CMIP mask):', land_perc.data.min(), '%') # check that min land area is 0 %

"""
(3.2) Mask out ocean in CMIP datasets (i.e. selecting only grid cells with > 50 % land)
"""
print('(3.2) Apply CMIP land-sea mask...')

CMIP_land_K = CMIP_slice_K.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
#numpy added a np.where function that allows us to simply use a logical command

CMIP_land_C = CMIP_slice_C.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
#numpy added a np.where function that allows us to simply use a logical command

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









########################################################################################################
"""
STEP 8: APPLY OBSERVATIONAL LAND-SEA MASK
"""
########################################################################################################










########################################################################################################
"""
STEP XX: OUTPUT RESULTS 
"""
########################################################################################################

"""
(X.1) Exporting the results to netcdf format
"""
print('Data export to NetCDF...')
# Export the land only data in Kelvin
land_dataDIR_K = r'G:\Climate_Data\4_Python_Scripting\2_RF\3_NorESM2_MM\Applying_land_sea_mask\clim_tas_ssp5_2090_2099_K_land.nc'
CMIP_land_K.to_netcdf(land_dataDIR_K)

# Export the land only data in Celsius
land_dataDIR_C = r'G:\Climate_Data\4_Python_Scripting\2_RF\3_NorESM2_MM\Applying_land_sea_mask\clim_tas_ssp5_2090_2099_C_land.nc'
CMIP_land_C.to_netcdf(land_dataDIR_C)

"""
(X.2) Output climate values as .csv
"""
#print('Data export to .csv...')

#########################################
print('End.')
