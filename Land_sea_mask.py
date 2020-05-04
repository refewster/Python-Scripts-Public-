# -*- coding: cp1252 -*-
"""
Program Purpose: To mask out ocean grid cells in surface temperature netcdfs from the CanESM5 netcdf.
Author: Richard Fewster
Start Date: 17/04/2020
Most Recent Update: 17/04/2020
GitHub version
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

"""
2) Create a list of required netCDF files.
"""
mask_file = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\2_CanESM5\land\sftlf_fx_CanESM5_hist-volc_r1i1p1f1_gn.nc"
climate_file = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\2_CanESM5\tmp\historical\tas_Amon_CanESM5_historical_r1i1p1f1_gn_185001-201412.nc"

"""
Identify the dimensions of each dataset:
"""
mask_data = Dataset(mask_file)
climate_data = Dataset(climate_file)

print('Mask file dimensions:', mask_data.dimensions)
print('Climate file dimensions:',climate_data.dimensions)

"""
Identify the variables stored within each dataset:
"""
print('Mask file variables:', mask_data.variables)
print('Climate file variables:', climate_data.variables)

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
3.2) Second we will load in the climate data and assign temperature (tas) to a new object.
"""
climate_dset = xr.open_dataset(climate_file)
clim = climate_dset['tas'].mean('time', keep_attrs=True) # create a new object with the temperature data for each grid cell
# keep_attrs maintains attribute data (e.g. units, variable names etc) in new object

"""
3.3) Masking process
The process below is a loop, that technically works but is not easy to read. It is preferable to use the second approach
- a one line command. 
"""
#nlats, nlons = clim.data.shape
#for y in range(nlats):
 #   for x in range(nlons):
  #      if land_perc.data[y, x] > 50:
    #        clim.data[y, x] = np.nan

"""
Masking out ocean (i.e. selecting only grid cells with > 50 % land)
"""
clim_land = clim.where(land_perc.data > 50.) # selects all grid cells where land % is less than 50 %
#numpy added a np.where function that allows us to simply use a logical command
print(clim_land)

"""
Masking out land (i.e. selecting only grid cells with < 50 % land)
"""
clim_ocean = clim.where(land_perc.data < 50.) # selects all grid cells where land % is less than 50 %
#numpy added a np.where function that allows us to simply use a logical command
# can select a different value for this if necessary
print(clim_ocean)

"""
4) Exporting the results to netcdf format
"""
# Export the land only data
land_dataDIR = r'G:\Climate_Data\4_Python_Scripting\2_RF\2_CanESM5\Applying_land_sea_mask\clim_land.nc'
clim_land.to_netcdf(land_dataDIR)

# Export the ocean only data
ocean_dataDIR = r'G:\Climate_Data\4_Python_Scripting\2_RF\2_CanESM5\Applying_land_sea_mask\clim_ocean.nc'
clim_ocean.to_netcdf(ocean_dataDIR)

#########################################
print('End.')
