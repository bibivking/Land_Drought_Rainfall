#!/usr/bin/python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"

import os
import numpy as np
from scipy.interpolate import CubicSpline
from netCDF4 import Dataset

# Load the NetCDF file
MODIS_LAI_path = '/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily/'
input_file     = f"{MODIS_LAI_path}MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2002-2024.nc"
output_file    = f"{MODIS_LAI_path}remove_high_frequent_varibility/MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2002-2024_cubic_spline_smooth.nc"

# Open the dataset
ds_in    = Dataset(input_file, 'r')
time     = ds_in.variables['time'][:]
lai_data = ds_in.variables['Lai_500m'][:, :, :]  # Assuming 'LAI' is the variable name and shape is (time, lat, lon)
lat      = ds_in.variables['latitude'][:]
lon      = ds_in.variables['longitude'][:]

# set nan as default value
lai_data = np.where(lai_data == 255, np.nan, lai_data)

# Prepare an empty array for the smoothed data
smoothed_lai = np.full_like(lai_data, np.nan)  # Initialize with NaNs for consistency

# Convert time to a numeric index for interpolation
time_numeric = np.arange(len(time))

# Loop through each (lat, lon) grid cell and apply cubic spline
for i in range(lai_data.shape[1]):      # Loop over latitude
    for j in range(lai_data.shape[2]):  # Loop over longitude
        y = lai_data[:, i, j]
        # Check if all values are NaN
        if not np.all(np.isnan(y)):
            # Mask NaN values for interpolation
            mask = ~np.isnan(y)
            # Fit a cubic spline to non-NaN data
            cs = CubicSpline(time_numeric[mask], y[mask])
            # Interpolate over the entire time range
            smoothed_lai[:, i, j] = cs(time_numeric)
        else:
            smoothed_lai[:, i, j] = np.nan  # Leave as NaN if all values are NaN

# ture nan to 255
smoothed_lai = np.where(np.isnan(smoothed_lai), 255, smoothed_lai)

# delete output file if it exists
if os.path.exists(output_file):
    os.remove(output_file)

# Write the smoothed data to a new NetCDF file
ds_out = Dataset(output_file, 'w')

# Create dimensions
ds_out.createDimension('time', len(time))
ds_out.createDimension('latitude', len(lat))
ds_out.createDimension('longitude', len(lon))

# Create variables and set attributes
time_out = ds_out.createVariable('time', 'f4', ('time',))
lat_out  = ds_out.createVariable('latitude', 'f4', ('latitude',))
lon_out  = ds_out.createVariable('longitude', 'f4', ('longitude',))
lai_out  = ds_out.createVariable('Lai_500m', 'f4', ('time', 'latitude', 'longitude'))

# Copy data into new variables
time_out[:]      = time
lat_out[:]       = lat
lon_out[:]       = lon
lai_out[:, :, :] = smoothed_lai

# Copy attributes if needed
time_out.setncatts(ds_in.variables['time'].__dict__)
lat_out.setncatts(ds_in.variables['latitude'].__dict__)
lon_out.setncatts(ds_in.variables['longitude'].__dict__)

lai_attributes = ds_in.variables['Lai_500m'].__dict__.copy()
lai_attributes.pop('_FillValue', None)  # Remove '_FillValue' if present
lai_out.setncatts(lai_attributes)
lai_out.setncattr('long_name', 'LAI smoothed with 31-day cubic spline')

print(f"Smoothed data saved to {output_file}")
ds_in.close()
ds_out.close()
