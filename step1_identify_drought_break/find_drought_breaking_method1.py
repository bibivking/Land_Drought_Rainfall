#!/usr/bin/python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"

import os
import gc
import sys
import glob
from tqdm import tqdm
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timedelta
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import matplotlib.ticker as mticker
from common_utils import *
from multiprocessing import Pool

def find_drought_breaking_date(spi_series, threshold, min_period):
    
    # Initialize the output array
    result = np.zeros_like(spi_series, dtype=np.int8)

    # Iterate over time starting from day 30 (to have a full 30-day history)
    for t in range(min_period, spi_series.shape[0]):
        
        # Extract the 30-day window for the current time step
        spi_window = spi_series[t - min_period:t]
        
        # Condition 1: SPI in the last 30 days didn't increase and < threshold
        condition_1 = np.all(np.diff(spi_window) <= 0, axis=0) & (spi_window[-1] < threshold)
        
        # Condition 2: SPI increases today compared to yesterday
        condition_2 = spi_series[t] > spi_series[t - 1]
        
        # Combine conditions
        result[t] = np.where(condition_1 & condition_2, 1, 0)

    return result

def process_grid_point(args):
    """
    Process a single (lat, lon) grid point to calculate accumulated days.
    """
    spi_series, threshold, min_period = args
    if not np.isnan(spi_series).all():  # Skip grid points with only NaN values
        return find_drought_breaking_date(spi_series, threshold, min_period)
    else:
        return np.full(len(spi_series),np.nan)  # Return None for NaN-only series

def find_drought_breaking_date_map(file_path, output_file, threshold=-1, min_period=30):

    # Load the NetCDF file
    with nc.Dataset(file_path, mode='r') as f:
        spi        = f.variables['spi_pearson_90'][:]
        time       = f.variables['time'][:]
        lat        = f.variables['lat'][:]
        lon        = f.variables['lon'][:]
        time_units = f.variables['time'].units

    # Handle missing values
    spi = np.where(spi == -9999., np.nan, spi)

    # Prepare arguments for parallel processing
    args_list   = [(spi[:, i, j], threshold, min_period) for i in range(len(lat)) for j in range(len(lon))]

    # Use multiprocessing to calculate accumulated days for each grid point
    with Pool() as pool:
        results = list(tqdm(pool.imap(process_grid_point, args_list), total=len(args_list)))

    # Prepare an empty array to store the results
    drought_break_days = np.zeros_like(spi, dtype=np.int8)

    # Fill the drought_break_days array with results
    for idx, result in enumerate(results):
        if result is not None:
            i, j = divmod(idx, len(lon))
            drought_break_days[:, i, j] = result

    drought_break_days = np.where(np.isnan(drought_break_days), -9999., drought_break_days)

    # Save the result to a new NetCDF file
    with nc.Dataset(output_file, 'w', format='NETCDF4') as f_out:
        # Create dimensions
        f_out.createDimension('time', len(time))
        f_out.createDimension('lat', len(lat))
        f_out.createDimension('lon', len(lon))

        # Create variables and assign data
        f_out.createVariable('time', 'f8', ('time',))[:] = time
        f_out.variables['time'].units                    = time_units
        f_out.createVariable('lat', 'f4', ('lat',))[:]   = lat
        f_out.createVariable('lon', 'f4', ('lon',))[:]   = lon

        drought_break_days_out       = f_out.createVariable('drought_break_days', 'i4', ('time', 'lat', 'lon'),-9999.)
        drought_break_days_out[:]    = drought_break_days
        drought_break_days_out.units = "1: Yes; 0: No"
        drought_break_days_out.long_name = (
            f"Drought-breaking date indicator: SPI < {threshold} for at least {min_period} days, "
            "with an increase on the current day."
        )

if __name__ == "__main__":

    file_path   = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_90_reorder_nan_filled.nc'  # Replace with the path to your file
    output_file = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_90_reorder_nan_filled_drought_breaking_date.nc'
    threshold   = -1
    min_period  = 30

    find_drought_breaking_date_map(file_path, output_file, threshold=threshold, min_period=min_period)
