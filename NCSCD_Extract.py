# -*- coding: cp1252 -*-
########################################################################################################
"""
Program: SOC data extract from netCDF
Description: This program reads in NCSCD netCDF files and outputs an excel file with the desired data. 
Author: Richard Fewster (gy15ref@leeds.ac.uk)
Start Date: 26/05/2020
Most Recent Update: 26/05/2020

"""
########################################################################################################
########################################################################################################
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
import scipy
import gridfill
import iris

print('Import complete')

"""
(1.2) Import required netCDF files and setup export directory.
"""
print('(1.2) Importing data files...')

histel = r"G:\GIS\SOC\CIRCUMPOLAR\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_WGS84_histel_pct_05deg.nc"
histel = xr.open_mfdataset(histel, combine='by_coords')

histosol = r"G:\GIS\SOC\CIRCUMPOLAR\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_WGS84_histosol_pct_05deg.nc"
histosol = xr.open_mfdataset(histosol, combine='by_coords')

# either opening method works fine! (changed it to test something).
SOCC_30= xr.open_dataset(r"G:\GIS\SOC\CIRCUMPOLAR\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_WGS84_SOCC30_05deg.nc")
SOCC_100= xr.open_dataset(r"G:\GIS\SOC\CIRCUMPOLAR\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_WGS84_SOCC100_05deg.nc")
SOCC_200= xr.open_dataset(r"G:\GIS\SOC\CIRCUMPOLAR\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_WGS84_SOCC200_05deg.nc")
SOCC_300= xr.open_dataset(r"G:\GIS\SOC\CIRCUMPOLAR\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_WGS84_SOCC300_05deg.nc")

print('(1.3) Setup export directory...')

# Export path 
DIR = r'G:\GIS\PEATLAND_MAP\NCSCD_data_'

########################################################################################################
"""
STEP 2: OUTPUT RESULTS 
"""
########################################################################################################


answer = input('(OPTIONAL) Export data files to .csv?:')
if answer.lower().startswith("y"):
      print("Data export to .csv...")
# Prevent warnings from flashing up - turn off/on as desired
      np.warnings.filterwarnings('ignore')
# Temperature files
# Historical
      histel_df= histel.to_dataframe() # turn this data into a pandas dataframe
      histel_df = histel_df.rename(columns={"NCSCDv2": "HISTEL_PERC"}) # drop unnecessary columns, rename variable columns to month
      print('#')
      histosol_df= histosol.to_dataframe() # turn this data into a pandas dataframe
      histosol_df = histosol_df.rename(columns={"NCSCDv2": "HISTOSOL_PERC"}) # drop unnecessary columns, rename variable columns to month
      print('##')      
      peat_df = pd.concat([histel_df, histosol_df], axis=1) # add each variable as a column
      peat_df = peat_df.reset_index() # add id column
      peat_df.index = peat_df.index + 1 # start id index at 1, not 0
      peat_df.to_csv(DIR+'PEAT_PRESENCE.csv')
      print('PEAT_PRESENCE.csv complete')



      SOCC_30_df=SOCC_30.to_dataframe() # turn this data into a pandas dataframe
      SOCC_30_df = SOCC_30_df.rename(columns={"NCSCDv2": "SOCC_30"}) # drop unnecessary columns, rename variable columns to month
      print('##')      
      SOCC_100_df=SOCC_100.to_dataframe() # turn this data into a pandas dataframe
      SOCC_100_df = SOCC_100_df.rename(columns={"NCSCDv2": "SOCC_100"}) # drop unnecessary columns, rename variable columns to month
      print('####')      
      SOCC_200_df=SOCC_200.to_dataframe() # turn this data into a pandas dataframe
      SOCC_200_df = SOCC_200_df.rename(columns={"NCSCDv2": "SOCC_200"}) # drop unnecessary columns, rename variable columns to month
      print('#####')      
      SOCC_300_df=SOCC_300.to_dataframe() # turn this data into a pandas dataframe
      SOCC_300_df = SOCC_300_df.rename(columns={"NCSCDv2": "SOCC_300"}) # drop unnecessary columns, rename variable columns to month
      print('######')      
      SOCC_df = pd.concat([SOCC_30_df, SOCC_100_df, SOCC_200_df, SOCC_300_df], axis=1) # add each variable as a column
      SOCC_df = SOCC_df.reset_index() # add id column
      SOCC_df.index = SOCC_df.index + 1 # start id index at 1, not 0
      SOCC_df.to_csv(DIR+'SOCC.csv')
      # For some reason SOCC_30 isnt sorted in the same way as the others - leading to duplication
      print('SOCC.csv complete')

 # Turn warnings back on
      np.warnings.filterwarnings('default')
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #    
     
elif answer.lower().startswith("n"):
    pass
else:
        print("Enter either yes/no")

print('Data Export complete')
################################################################################################################################
print('End.')
import sys
sys.exit()
