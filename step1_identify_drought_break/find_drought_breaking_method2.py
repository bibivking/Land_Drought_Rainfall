#!/usr/bin/python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"

import os
import gc
import sys
import glob
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timedelta
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import matplotlib.ticker as mticker
from multiprocessing import Pool
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature, OCEAN
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from common_utils import *

def find_drought_breaking_rainfall_pixel(precip_series, spi_series, spi_thres=-1, daily_rain_thres=5,
    drought_period=30, drought_period_tot_rain_thres=20, recover_period=30, spi_increase_thres=1):
    
    ntime = len(precip_series)

    if np.all(np.isnan(spi_series)):
        
        return np.full(ntime, np.nan)

    else:
        
        # condition 1: the rainfall intensity > 5 mm/day
        condition1 = precip_series > daily_rain_thres

        # condition 2: the day before the rainfall is in drought (<-1)
        condition2 = spi_series < spi_thres

        # condition 3 : total rainfall of perivous 30 days < 20 mm and no rainfall in the past 30 days > 5 mm/day
        condition3 = np.zeros(ntime, dtype=bool)

        # Handle the first 30 time points

        condition3[:drought_period]     = False

        # For the remaining time points, check the condition
        for t in np.arange(drought_period, ntime):
            if (( np.max(precip_series[t-drought_period:t]) < daily_rain_thres) 
            & (np.sum(precip_series[t-drought_period:t]) < drought_period_tot_rain_thres)):
                condition3[t] = True
            else:
                condition3[t] = False

        # condition 4 : drought is broken (SPI>-1) or SPI inceases by > 1 within 30 days after the rain 
        condition4 = np.zeros(ntime, dtype=bool)

        # Handle the last 14 time points separately
        condition4[:ntime-recover_period] = False

        # For the remaining time points, check the condition
        for t in range(ntime - recover_period-1):
            if ((np.max(spi_series[t+1:t+1+recover_period]) > spi_thres) 
            | (np.max(spi_series[t+1:t+1+recover_period])-spi_series[t-1] > spi_increase_thres)):
                condition4[t] = True

        conditions = (condition1) & (condition2) & (condition3) & (condition4)

        gc.collect()

        return conditions

def find_drought_breaking_rainfall_parallel(output_file, rain_file, spi_file, mask_file, spi_thres=-1, daily_rain_thres=5,
                                            drought_period=30, drought_period_tot_rain_thres=20, recover_period=30, 
                                            spi_increase_thres=1):

    # days from 1 Jan 1950
    date_s = (datetime(2000, 1, 1) - datetime(1950, 1, 1)).days

    # Read rainfall data
    with nc.Dataset(rain_file, mode='r') as f_rain:
        precip     = f_rain.variables['precip'][date_s:,10:,:-45]
    gc.collect()

    # Read SPI data
    with nc.Dataset(spi_file, mode='r') as f_spi:
        lat        = f_spi.variables['lat'][:]
        lon        = f_spi.variables['lon'][:]
        spi        = f_spi.variables['spi_pearson_90'][date_s:,:,:]
        time       = nc.num2date(f_spi.variables['time'][date_s:],f_spi.variables['time'].units,
                                 only_use_cftime_datetimes=False,only_use_python_datetimes=True)
        time_val   = f_spi.variables['time'][date_s:]
        ntime      = len(time)
        nlat       = len(lat)
        nlon       = len(lon)
        time_units = f_spi.variables['time'].units
    gc.collect()

    # Read land sea mask    
    with nc.Dataset(mask_file, mode='r') as f_land:
        landsea    = f_land.variables['landsea'][:]
        landsea_3d = np.repeat(landsea[np.newaxis, :, :], ntime, axis=0)
    gc.collect()

    spi_land    = np.where(landsea_3d == 0, spi, np.nan)
    precip_land = np.where(landsea_3d == 0, precip, np.nan)

    # Prepare arguments for each variable and soil layer
    args_list = []
    for i in np.arange(nlat):
        for j in np.arange(nlon):  # Assuming 6 soil layers
            # print('spi_land[:,i,j] is', spi_land[:, i, j])
            args_list.append((precip_land[:, i, j], 
                              spi_land[:, i, j],
                              spi_thres, 
                              daily_rain_thres, 
                              drought_period, 
                              drought_period_tot_rain_thres, 
                              recover_period, 
                              spi_increase_thres))

    # Use multiprocessing to process each soil layer of each variable
    with Pool() as pool:
        results = pool.starmap(find_drought_breaking_rainfall_pixel, args_list)
    
    print('shape of results is', np.shape(results))
    results_array    = np.array(results)
    reshaped_results = results_array.reshape(nlat, nlon, ntime)

    drought_breakings = np.zeros((ntime,nlat,nlon))

    # Organize regridded data and finalize LIS file
    for i in np.arange(nlat):
        for j in np.arange(nlon):
            drought_breakings[:,i,j] = reshaped_results[i,j,:]
    
    # Save the result to a new NetCDF file
    with nc.Dataset(output_file, 'w', format='NETCDF4') as f_out:

        # Create dimensions
        f_out.createDimension('time', len(time))
        f_out.createDimension('lat', len(lat))
        f_out.createDimension('lon', len(lon))

        # Create variables and assign data
        f_out.createVariable('time', 'f8', ('time',))[:] = time_val
        f_out.variables['time'].units                    = time_units
        f_out.createVariable('lat', 'f4', ('lat',))[:]   = lat
        f_out.createVariable('lon', 'f4', ('lon',))[:]   = lon

        drought_break_days_out       = f_out.createVariable('is_drought_break', 'i4', ('time', 'lat', 'lon'), fill_value=-9999.)
        drought_break_days_out[:]    = drought_breakings
        drought_break_days_out.units = "1: Yes; 0: No"
        drought_break_days_out.long_name = (
            f"condition 1: the rainfall intensity > {daily_rain_thres} mm/day"
            f"condition 2: the day before the rainfall is in drought (<{daily_rain_thres})"
            f"condition 3: total rainfall of perivous {drought_period} days < {drought_period_tot_rain_thres} mm and no rainfall in the past {drought_period} days > {daily_rain_thres} mm/day"
            f"condition 4: drought is broken (SPI> {daily_rain_thres}) or SPI inceases by > {spi_increase_thres} within {drought_period} days after the rain "
        )

    return


if __name__ == "__main__":

    rain_file = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/agcd_v1_precip_total_r005_daily_1950_2023.nc'
    spi_file  = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_90_reorder_nan_filled.nc'  # Replace with the path to your file 
    mask_file = '/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_mask.nc'
    output_file = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_90_reorder_nan_filled_drought_breaking_date_method2_condition2.nc'
    
    # condition 1 & 3
    daily_rain_thres              = 5
    
    # condition 2
    spi_thres                     = -1

    # condition 3
    drought_period                = 15 #30
    drought_period_tot_rain_thres = 10 #20

    # condition 4
    recover_period                = 30
    spi_increase_thres            = 1

    find_drought_breaking_rainfall_parallel(output_file, rain_file, spi_file, mask_file, spi_thres, daily_rain_thres,
                                            drought_period, drought_period_tot_rain_thres, recover_period, spi_increase_thres)
