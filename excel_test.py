# -*- coding: cp1252 -*-
########################################################################################################
"""
Program: Excel extraction

Author: Richard Fewster
Start Date: 19/05/2020
Most Recent Update: 19/05/2020

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
import os
print('Import complete')

"""
(1.2) Create a list of required netCDF files.
"""
print('(1.2) Importing data files...')
# CMIP land-sea mask file
mask_file = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\land\sftlf_fx_NorESM2-MM_historical_r1i1p1f1_gn.nc"
# Temperature files
# For models with multiple decadal files, put an asterisk after the folder name e.g. "..historical\*.nc"
tas_file_ssp1 = r"G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\downscaled_outputs\NorESM2_downscaled_monthly_tas_bcor_ssp1.nc"
ssp1 = xr.open_mfdataset(tas_file_ssp1, combine='by_coords')
ssp1_lat = ssp1.sel(lat=slice(50., 90.))

# Export path for temperature files
tmp_DIR = r'G:\Climate_Data\1_CMIP_DATA\2_CMIP6\1_NorESM2_MM\downscaled_outputs\excel_test\test_'

df= ssp1_lat.to_dataframe() # produces a dataframe with 12 rows for each month for each site.
# slicing to one month creates just one list of x,y,z - append on new months after? 


jan = ssp1.sel(month=slice(1))
df= jan.to_dataframe()

jan = ssp1.sel(month=1)
jan_df= jan.to_dataframe()
jan_df = jan_df.drop(["month", "height"], axis=1).rename(columns={"mean monthly near-surface temperature (degrees Celsius)": "Jan_MMT"})

feb = ssp1.sel(month=2)
feb_df= feb.to_dataframe()
feb_df = feb_df.drop(["month", "height"], axis=1).rename(columns={"mean monthly near-surface temperature (degrees Celsius)": "Feb_MMT"})

mar = ssp1.sel(month=3)
mar_df= mar.to_dataframe()
mar_df = mar_df.drop(["month", "height"], axis=1).rename(columns={"mean monthly near-surface temperature (degrees Celsius)": "Mar_MMT"})

apr = ssp1.sel(month=4)
apr_df= apr.to_dataframe()
apr_df = apr_df.drop(["month", "height"], axis=1).rename(columns={"mean monthly near-surface temperature (degrees Celsius)": "Apr_MMT"})

may = ssp1.sel(month=5)
may_df= may.to_dataframe()
may_df = may_df.drop(["month", "height"], axis=1).rename(columns={"mean monthly near-surface temperature (degrees Celsius)": "May_MMT"})

jun = ssp1.sel(month=6)
jun_df= jun.to_dataframe()
jun_df = jun_df.drop(["month", "height"], axis=1).rename(columns={"mean monthly near-surface temperature (degrees Celsius)": "Jun_MMT"})

jul = ssp1.sel(month=7)
jul_df= jul.to_dataframe()
jul_df = jul_df.drop(["month", "height"], axis=1).rename(columns={"mean monthly near-surface temperature (degrees Celsius)": "Jul_MMT"})

aug = ssp1.sel(month=8)
aug_df= aug.to_dataframe()
aug_df = aug_df.drop(["month", "height"], axis=1).rename(columns={"mean monthly near-surface temperature (degrees Celsius)": "Aug_MMT"})

sep = ssp1.sel(month=9)
sep_df= sep.to_dataframe()
sep_df = sep_df.drop(["month", "height"], axis=1).rename(columns={"mean monthly near-surface temperature (degrees Celsius)": "Sep_MMT"})

oct = ssp1.sel(month=10)
oct_df= oct.to_dataframe()
oct_df = oct_df.drop(["month", "height"], axis=1).rename(columns={"mean monthly near-surface temperature (degrees Celsius)": "Oct_MMT"})

nov = ssp1.sel(month=11)
nov_df= nov.to_dataframe()
nov_df = nov_df.drop(["month", "height"], axis=1).rename(columns={"mean monthly near-surface temperature (degrees Celsius)": "Nov_MMT"})

dec = ssp1.sel(month=12)
dec_df= dec.to_dataframe()
dec_df = dec_df.drop(["month", "height"], axis=1).rename(columns={"mean monthly near-surface temperature (degrees Celsius)": "Dec_MMT"})
import sys
#while True:
answer = input('Do you want to continue?:')
if answer.lower().startswith("y"):
      print("ok, carry on then")
      df_col = pd.concat([jan_df,feb_df, mar_df, apr_df, may_df, jun_df, jul_df, aug_df, sep_df, oct_df, nov_df, dec_df], axis=1)
      df_col = df_col.reset_index()
      df_col.index = df_col.index + 1
      df_col.to_csv(tmp_DIR)
elif answer.lower().startswith("n"):
      print("ok, sayonnara")
      sys.exit()




# Find way to add id column?
