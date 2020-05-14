# -*- coding: cp1252 -*-
########################################################################################################
"""
Program: Downscaling Program (DRAFT)

Author: Richard Fewster
Start Date: 30/04/2020
Most Recent Update: 12/05/2020

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
# Temperature
# GroupBy subdivides dataset into months before averaging. This code prodcues monthly mean temperatures.
tas_hist_mean_monthly_K = tas_hist_slice['tas'].groupby("time.month").mean('time',keep_attrs=True)
tas_ssp1_mean_monthly_K = tas_ssp1_slice['tas'].groupby("time.month").mean('time',keep_attrs=True)

"""
(2.2) Convert climate data into desired units
"""
# Convert from Kelvin to Celsius
tas_hist_mean_monthly_C = tas_hist_mean_monthly_K-273.15
tas_ssp1_mean_monthly_C = tas_ssp1_mean_monthly_K-273.15

print('Step 2 Complete')

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
mask_dset = xr.open_dataset(mask_file)#Use xarray to open the mask dataset
mask_dset_slice = mask_dset.sel(lat=slice(50., 90.)) 
land_perc = mask_dset_slice['sftlf'] # assign the land percentage variable to a new object


print('Max land area (CMIP mask):', land_perc.data.max(), '%') # check that max land area is 100 %
print('Min land area (CMIP mask):', land_perc.data.min(), '%') # check that min land area is 0 %

"""
(3.2) Mask out ocean in CMIP datasets (i.e. selecting only grid cells with > 50 % land)
"""
print('(3.2) Apply land-sea mask...')
#numpy includes a np.where function that allows us to simply use a logical command

# Mask out temperature data
tas_hist_land_C = tas_hist_mean_monthly_C.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
#tas_ssp1_land_C = tas_ssp1_mean_monthly_C.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %

# Mask out preciptiation data
#pre_hist_land_C = pre_hist_mean_monthly_mm.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
#pre_ssp1_land_C = pre_ssp1_mean_monthly_mm.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %

print('Step 3 complete')

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
CRU_tmp_file = r"G:\Climate_Data\3_Observational_data\CRU data\CRU_TS_404\cru_ts4.04.1901.2019.tmp.dat.nc"
CRU_tmp_dset = xr.open_mfdataset(CRU_tmp_file, combine='by_coords')
# Load in observational precipitation dataset
# CRU_pre_file =r"G:\Climate_Data\3_Observational_data\CRU data\CRU_TS_404\cru_ts4.04.1901.2019.pre.dat.nc"
#CRU_pre_dset = xr.open_mfdataset(CRU_pre_file, combine='by_coords')

"""
"""
# REMEMBER: CRU datasets use -179.75 -> 179.75 for lon
CRU_tmp_slice = CRU_tmp_dset.sel(time=slice("1961-01-16", "1990-12-16"), lat=slice(50., 90.)).drop(labels='stn') # Slice to match study region of SSP files
# Drop 'stn' variable

#CRU_pre_slice = CRU_pre_dset.sel(time=slice("1961-01-16", "1990-12-16"), lat=slice(50., 90.)) # Slice to match study region of SSP files

print('Step 4 complete')

########################################################################################################
"""
STEP 5: EXTRAPOLATION 
Need to remove NaNs for interpolation step below.
Use Poisson Equation solver with overrelaxation to extrapolate terrestrial data over ocean.
"""
########################################################################################################
# 

# Use grid fill to extrapolate over ocean
# Can  install gridfill manually from https://github.com/ajdawson/gridfill

# Create dictionary with settings for Poisson grid filling. Definitions:
# eps = Tolerance for determining the solution complete. [1e-3]
# relax = Relaxation constant. Usually 0.45 <= *relax* <= 0.6. Defaults to 0.6. [0.6]
# itermax = Maximum number of iterations of the relaxation scheme. Defaults to 100 iterations. [2000]
# initzonal = If *False* missing values will be initialized to zero, if *True* missing values will be initialized to the zonal mean. Defaultsto *False*. [True]
# cyclic = Set to *False* if the x-coordinate of the grid is not cyclic, set to *True* if it is cyclic. Defaults to *False*. [Not used]
# verbose = If *True* information about algorithm performance will be printed to stdout, if *False* nothing is printed. Defaults to *False*. [True]

# Try iris instead of xarray:
tas_hist_land_Ciris = tas_hist_land_C.to_iris() # convert to Iris cube
# Need to add a cyclic attribute to the iris cube itself
tas_hist_land_Ciris.coord('longitude').circular = True

tas_hist_land_Ciris_backfilled = gridfill.fill_cube(tas_hist_land_Ciris, 1e-3, 0.6, 2000, initzonal=True, verbose=True)

print('Step 5 complete')

########################################################################################################
"""
STEP 6: REGRIDDING TO HIGHER RESOLUTION
Base new grid on CRU grid (0.5 x 0.5 deg). Interpolate using bicubic spline.
"""
# Giving up with xarray, use iris
# ("1961-01-16", "1990-12-16"), lat=slice(50., 90.))
CRU_tmp_array = iris.load(CRU_tmp_file)[1] # only tmp not stn
# Select dates (not actually required for regridding, but may be useful later?)
import datetime
date1 = datetime.datetime.strptime('19610116T0000Z','%Y%m%dT%H%MZ')
date2 = datetime.datetime.strptime('19901216T0000Z','%Y%m%dT%H%MZ')
date_range = iris.Constraint(time=lambda cell: date1 <= cell.point <= date2 )
# Select latitudes
lat1 = 50.
lat2 = 90.
lat_range = iris.Constraint(latitude=lambda cell: lat1 <= cell.point <= lat2 )
CRU_tmp_array_t_slice = CRU_tmp_array.extract(date_range)
CRU_tmp_array_slice = CRU_tmp_array_t_slice.extract(lat_range)


print("(6) Performing linear interpolation with Iris...")
# Bilinear interpolation using CRU_tmp_array_slice grid
tas_hist_land_Ciris_backfilled_high = tas_hist_land_Ciris_backfilled.regrid(CRU_tmp_array_slice, iris.analysis.Linear())
print("Step 6 complete")

########################################################################################################
"""
STEP 7: APPLY OBSERVATIONAL LAND-SEA MASK
"""
########################################################################################################
"""
(7.1) Create CRU land sea mask
"""
## Convert back to xarray (if wanted):
tas_hist_land_C_backfilled_high = xr.DataArray.from_iris(tas_hist_land_Ciris_backfilled_high)
CRU_tmp_xr = xr.DataArray.from_iris(CRU_tmp_array_slice)

# Make CRU mask same shape as CMIP tas array
CRU_mask_xr = CRU_tmp_xr.groupby("time.month").mean('time',keep_attrs=True)

"""
(7.2) Apply CRU land sea mask to downscaled CMIP data
"""
print("Applying observational land sea mask to downscaled CMIP data...")

# The fill value for missing values in the CRU data is -999. This line selects only those which are greater than that value.
tas_hist_land_C_backfilled_high_masked = tas_hist_land_C_backfilled_high.where(CRU_mask_xr.data >-998) # Mask precipitation files...
# Provides a runtime warning - solve this? Masking does appear to have worked though.


#ssp 1


print("Step 7 complete")


xx


########################################################################################################
"""
STEP 8: BIAS CORRECTION
"""
########################################################################################################
"""
(7.3) Variable formatting
"""
#CMIP_fut_tmp = tas_SSP1_land_backfilled_high_C
#CMIP_hist_tmp = tas_hist_land_Ciris_backfilled_high
#CRU_tmp = CRU_tmp_slice

"""
(7.2) Temperature Bias correction

BCor_Temp = (CMIP_fut_tmp - CMIP_hist_tmp) + CRU_tmp
"""
# Replace CRU dummy variables with CMIP equivalents
SSP1_BCor_Temp = (SSP1_fut_tmp - CMIP_hist_tmp) + CRU_tmp

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
tas_hist_land_backfilled_high_C.to_netcdf(DIR+'hist.nc')
SSP1_BCor_Temp_land.to_netcdf(DIR+'SSP1.nc')


"""
(X.2) Output monthly climate values for each 0.5 degree grid cell as .csv
"""
#print('Data export to .csv...')








#########################################
print('End.')
