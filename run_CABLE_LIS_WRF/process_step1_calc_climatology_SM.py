#!/usr/bin/python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"

from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.interpolate import griddata
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from common_utils import *

# Create new nc file to save climatology SM 
file_path = "/g/data/w97/mm3972/model/cable/runs/Land_drought_rainfall_runs/spinup_run_1970_1999/outputs"

for year in np.arange(1970,2000):
    file_in   = f'{file_path}/cable_out_{year}.nc'
    
    SoilMoist = off_file.variables['SoilMoist'][:]
    SoilMoistIce s