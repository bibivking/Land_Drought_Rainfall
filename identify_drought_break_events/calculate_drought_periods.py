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
from common_utils import *
from multiprocessing import Pool

def calculate_accumulated_days(spi_series, threshold, min_period):
    """
    Function to calculate accumulated days for a given SPI series.
    """
    condition               = spi_series < threshold
    accumulated_days_series = np.zeros_like(spi_series, dtype=int)

    count = 0
    for i, val in enumerate(condition):
        if val:  # SPI is below the threshold
            count += 1
        else:
            if count >= min_period:
                accumulated_days_series[i - count:i] = np.arange(1, count + 1)
            count = 0

    # Handle cases where the period extends to the last time step
    if count >= min_period:
        accumulated_days_series[-count:] = np.arange(1, count + 1)

    return accumulated_days_series

def process_grid_point(args):
    """
    Process a single (lat, lon) grid point to calculate accumulated days.
    """
    spi_series, threshold, min_period = args
    if not np.isnan(spi_series).all():  # Skip grid points with only NaN values
        return calculate_accumulated_days(spi_series, threshold, min_period)
    else:
        return None  # Return None for NaN-only series

def calculate_drought_periods(file_path, output_file, threshold=-1, min_period=30):
    """
    Function to calculate drought periods and save the results in a new NetCDF file.
    """

    # Load the NetCDF file
    f = nc.Dataset(file_path, mode='r')

    # Extract the necessary variables
    spi  = f.variables['spi_pearson_90'][:]
    time = f.variables['time'][:]
    lat  = f.variables['lat'][:]
    lon  = f.variables['lon'][:]

    # set nan value
    spi  = np.where(spi==-9999., np.nan, spi)

    # Prepare an empty array to store the results
    accumulated_days = np.zeros_like(spi)

    # Prepare arguments for parallel processing
    args_list = [(spi[:, i, j], threshold, min_period) for i in range(lat.shape[0]) for j in range(lon.shape[0])]

    # Use multiprocessing to calculate accumulated days for each grid point
    with Pool() as pool:
        results = pool.map(process_grid_point, args_list)

    # Fill the accumulated_days array with results
    for idx, result in enumerate(results):
        if result is not None:
            i, j = divmod(idx, lon.shape[0])
            accumulated_days[:, i, j] = result

    # Save the result to a new NetCDF file
    f_out = nc.Dataset(output_file, 'w', format='NETCDF4')

    # Create dimensions
    f_out.createDimension('time', len(time))
    f_out.createDimension('lat', len(lat))
    f_out.createDimension('lon', len(lon))

    # Create variables
    time_out             = f_out.createVariable('time', 'f8', ('time',))
    lat_out              = f_out.createVariable('lat', 'f4', ('lat',))
    lon_out              = f_out.createVariable('lon', 'f4', ('lon',))
    accumulated_days_out = f_out.createVariable('accumulated_days', 'i4', ('time', 'lat', 'lon'))

    # Assign data to variables
    time_out.standard_name   = "time"
    time_out.long_name       = "time"
    time_out.bounds          = "time_bnds"
    time_out.axis            = "T"
    time_out.units           = "days since 1850-01-01"
    time_out.calendar        = "standard"
    time_out[:]              = time

    lon_out.standard_name    = "longitude"
    lon_out.long_name        = "longitude"
    lon_out.units            = "degrees_east"
    lon_out.axis             = "X"
    lon_out.bounds           = "lon_bnds"
    lon_out[:]               = lon

    lat_out.standard_name    = "latitude"
    lat_out.long_name        = "latitude"
    lat_out.units            = "degrees_north"
    lat_out.axis             = "Y"
    lat_out.bounds           = "lat_bnds"
    lat_out[:]               = lat

    # Assign attributes to the new variable
    accumulated_days_out.units     = "days"
    accumulated_days_out.long_name = f"Accumulated days with SPI < {threshold} for at least {min_period} days"
    accumulated_days_out[:]        = accumulated_days

    # Close the NetCDF files
    f_out.close()
    f.close()

    print(f"Accumulated days for SPI < {threshold} and period >= {min_period} days saved to {output_file}")
    return

if __name__ == "__main__":

    file_path   = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_90_reorder_nan_filled.nc'  # Replace with the path to your file
    output_file = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_90_reorder_nan_filled_drought_periods.nc'
    threshold   = -1
    min_period  = 30

    calculate_drought_periods(file_path, output_file, threshold=threshold, min_period=min_period)
