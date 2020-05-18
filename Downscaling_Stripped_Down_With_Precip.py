# -*- coding: cp1252 -*-
########################################################################################################
"""
Program: Downscaling Program (DRAFT) + precip

Author: Richard Fewster
Start Date: 30/04/2020
Most Recent Update: 18/05/2020

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
# For models with multiple decadal files, put an asterisk after the folder name e.g. "..historical\*.nc"
tas_file_hist = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\2_CanESM5\tmp\historical\tas_Amon_CanESM5_historical_r1i1p1f1_gn_185001-201412.nc"
tas_file_hist = xr.open_mfdataset(tas_file_hist)
tas_file_ssp1 = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\2_CanESM5\tmp\ssp1_26\tas_Amon_CanESM5_ssp126_r1i1p1f1_gn_201501-210012.nc"
tas_file_ssp1 = xr.open_mfdataset(tas_file_ssp1)

# Precipitation files
# For models with multiple decadal files, put an asterisk after the folder name e.g. "..historical\*.nc"
pre_file_hist = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\2_CanESM5\pre\historical\pr_Amon_CanESM5_historical_r1i1p1f1_gn_185001-201412.nc"
pre_file_hist = xr.open_mfdataset(pre_file_hist)
pre_file_ssp1 = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\2_CanESM5\pre\ssp1_26\pr_Amon_CanESM5_ssp126_r1i1p1f1_gn_201501-210012.nc"
pre_file_ssp1 = xr.open_mfdataset(pre_file_ssp1)

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
tas_hist_slice = tas_file_hist.sel(time=slice("1961-01-16", "1990-12-16"))
tas_ssp1_slice = tas_file_ssp1.sel(time=slice("2090-01-16", "2099-12-16")) 

# Precipitation slices
pre_hist_slice = pre_file_hist.sel(time=slice("1961-01-16", "1990-12-16"))
pre_ssp1_slice = pre_file_ssp1.sel(time=slice("2090-01-16", "2099-12-16")) 

"""
(2.2) Convert climate data into desired units and average climate values for study time period
"""
# Temperature
# GroupBy subdivides dataset into months before averaging. This code produces monthly mean temperatures.
tas_hist_mean_monthly_K = tas_hist_slice['tas'].groupby("time.month").mean('time',keep_attrs=True)
tas_ssp1_mean_monthly_K = tas_ssp1_slice['tas'].groupby("time.month").mean('time',keep_attrs=True)

# Convert from Kelvin to Celsius
tas_hist_mean_monthly_C = tas_hist_mean_monthly_K-273.15
tas_ssp1_mean_monthly_C = tas_ssp1_mean_monthly_K-273.15


# Precipitation
# Convert from mm/second to mm per day
# 60 x 60 x 24 = 86,400 (one day)
pre_hist_day = pre_hist_slice['pr'] * 86400.
pre_ssp1_day = pre_ssp1_slice['pr'] * 86400.

# Calculate the number of days in each month
# adapted from http://xarray.pydata.org/en/stable/examples/monthly-means.html
# Build a dictionary of calendars
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
pre_hist_month_length = xr.DataArray(get_dpm(pre_hist_day.time.to_index(), calendar='noleap'), coords=[pre_hist_day.time], name='month_length')
pre_ssp1_month_length = xr.DataArray(get_dpm(pre_ssp1_day.time.to_index(), calendar='noleap'), coords=[pre_ssp1_day.time], name='month_length')
# Multiply the days in the month by the precipitation value for that month
pre_hist_mth= pre_hist_month_length * pre_hist_day
pre_ssp1_mth= pre_ssp1_month_length * pre_ssp1_day

# Average monthly precipitation values by month
pre_hist_mean_mth = pre_hist_mth.groupby("time.month").mean('time',keep_attrs=True)
pre_ssp1_mean_mth = pre_ssp1_mth.groupby("time.month").mean('time',keep_attrs=True)

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

# Mask out preciptiation data
pre_hist_land = pre_hist_mean_mth.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
pre_ssp1_land = pre_ssp1_mean_mth.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %

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
CRU_pre_file =r"G:\Climate_Data\3_Observational_data\CRU data\CRU_TS_404\cru_ts4.04.1901.2019.pre.dat.nc"
CRU_pre_dset = xr.open_mfdataset(CRU_pre_file, combine='by_coords')

print('Step 4 complete')

########################################################################################################
"""
STEP 5: EXTRAPOLATION 
Need to remove NaNs for interpolation step below.
Use Poisson Equation solver with overrelaxation to extrapolate terrestrial data over ocean.
"""
########################################################################################################

xxxxxxxxxxxxxxxx Got up to here
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
tas_ssp1_land_Ciris = tas_ssp1_land_C.to_iris() # convert to Iris cube
# Need to add a cyclic attribute to the iris cube itself
tas_hist_land_Ciris.coord('longitude').circular = True
tas_ssp1_land_Ciris.coord('longitude').circular = True

tas_hist_land_Ciris_backfilled = gridfill.fill_cube(tas_hist_land_Ciris, 1e-3, 0.6, 2000, initzonal=True, verbose=True)
tas_ssp1_land_Ciris_backfilled = gridfill.fill_cube(tas_ssp1_land_Ciris, 1e-3, 0.6, 2000, initzonal=True, verbose=True)
# Increase iterations to ~7000 for NorESM2
print('Step 5 complete')

########################################################################################################
"""
STEP 6: REGRIDDING TO HIGHER RESOLUTION
Base new grid on CRU grid (0.5 x 0.5 deg). Interpolate using bicubic spline.
"""
"""
(6.1) Convert observational dataset to iris and reshape
"""
# Convert to iris for this step
# ("1961-01-16", "1990-12-16"), lat=slice(50., 90.))
CRU_tmp_array = iris.load(CRU_tmp_file)[1] # only tmp not stn
# Select dates (not actually required for regridding, but may be useful later?)
import datetime
date1 = datetime.datetime.strptime('19610116T0000Z','%Y%m%dT%H%MZ')
date2 = datetime.datetime.strptime('19901216T0000Z','%Y%m%dT%H%MZ')
date_range = iris.Constraint(time=lambda cell: date1 <= cell.point <= date2 )
CRU_tmp_array_slice = CRU_tmp_array.extract(date_range)

"""
(6.2) Perform bilinear interpolation
"""
print("(6) Performing linear interpolation with Iris...")
# Bilinear interpolation using CRU_tmp_array_slice grid
tas_hist_land_Ciris_backfilled_high = tas_hist_land_Ciris_backfilled.regrid(CRU_tmp_array_slice, iris.analysis.Linear())
tas_ssp1_land_Ciris_backfilled_high = tas_ssp1_land_Ciris_backfilled.regrid(CRU_tmp_array_slice, iris.analysis.Linear())

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
tas_ssp1_land_C_backfilled_high = xr.DataArray.from_iris(tas_ssp1_land_Ciris_backfilled_high)
#CRU_tmp_xr = xr.DataArray.from_iris(CRU_tmp_array_slice)


# Make CRU mask same shape as CMIP tas array
#CRU_mask_xr = CRU_tmp_xr.groupby("time.month").mean('time',keep_attrs=True)

# Iris slightly dodgy with variable selection - not an issue for previous step
# Reload CRU data as an xarray and select tmp variable
#Reshape to match shape of CMIP tas array
CRU_tmp_slice_xr = CRU_tmp_dset.sel(time=slice("1961-01-16", "1990-12-16"))
CRU_tmp_xr = CRU_tmp_slice_xr['tmp'].groupby("time.month").mean('time',keep_attrs=True)


"""
(7.2) Apply CRU land sea mask to downscaled CMIP data
"""
print("Applying observational land sea mask to downscaled CMIP data...")

# The fill value for missing values in the CRU data is -999. This line selects only those which are greater than that value.
tas_hist_land_C_backfilled_high_masked = tas_hist_land_C_backfilled_high.where(CRU_tmp_xr.data >-998) # Mask precipitation files...
# Provides a runtime warning - solve this? Masking does appear to have worked though.
# Possible that fill value is higher value (inf?)
#np.warnings.filterwarnings('ignore')
#np.seterr(divide='ignore', invalid='ignore')

tas_ssp1_land_C_backfilled_high_masked = tas_ssp1_land_C_backfilled_high.where(CRU_tmp_xr.data >-998) # Mask precipitation files...

print("Step 7 complete")

########################################################################################################
"""
STEP 8: BIAS CORRECTION
"""
########################################################################################################
"""
(8.1) Variable formatting
"""
CMIP_fut_tmp = tas_ssp1_land_C_backfilled_high_masked
CMIP_hist_tmp = tas_hist_land_C_backfilled_high_masked
CRU_tmp = CRU_tmp_xr

"""
(8.2) Temperature Bias correction

BCor_Temp = (CMIP_fut_tmp - CMIP_hist_tmp) + CRU_tmp
"""

SSP1_BCor_tmp = (CMIP_fut_tmp - CMIP_hist_tmp) + CRU_tmp
SSP1_BCor_tmp = SSP1_BCor_tmp.rename('near-surface temperature (degrees Celsius)') # Rename variable
"""
(8.3) Precipitation Bias correction

BCor_Pre = (CRU_pre / CMIP_hist_pre) * CMIP_fut_pre
"""
# Calculate alpha...
#a = CRU_pre / CMIP_hist_pre
# Set limits for alpha...
#a = xr.where(a < 0.25, 0.25, a)
#a_ltd = xr.where(a > 4.0, 4.0, a)

# Apply the limited alpha coefficient to bias correct future precipitation
#SSP1_BCor_pre = a_ltd * SSP1_fut_pre


"""
(8.4) Subset to the region of interest
"""
# Select only grid cells within the latitudinal bands
#SSP1_BCor_tmp_lat_slice = SSP1_BCor_tmp.sel(lat=slice(50., 90.))
#CMIP_hist_tmp_lat_slice = CMIP_hist_tmp.sel(lat=slice(50., 90.))
#CRU_tmp_lat_slice = CRU_tmp.sel(lat=slice(50., 90.))

print("Step 8 complete")

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
DIR = r'G:\Climate_Data\1_CMIP_DATA\2_CMIP6\2_CanESM5\downscaled_outputs\tas'

#SSP1_BCor_tmp_lat_slice.to_netcdf(DIR+'CanESM5_downscaled_tas_SSP1_BCor.nc')
#CMIP_hist_tmp_lat_slice.to_netcdf(DIR+'CanESM5_downscaled_tas_historical.nc')
#CRU_tmp_lat_slice.to_netcdf(DIR+'CRU_tmp.nc')


CMIP_fut_tmp.to_netcdf(DIR+'fut_cmip.nc')
CMIP_hist_tmp.to_netcdf(DIR+'hist_cmip.nc')
CRU_tmp.to_netcdf(DIR+'cru_tmp_full.nc')
SSP1_BCor_tmp.to_netcdf(DIR+'bcor_fut_cmip.nc')



"""
(X.2) Output monthly climate values for each 0.5 degree grid cell as .csv
"""
#print('Data export to .csv...')








#########################################
print('End.')
