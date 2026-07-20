#!/usr/bin/python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"

import numpy as np
from scipy.interpolate import CubicSpline
from netCDF4 import Dataset, num2date
import xarray as xr
import os
import matplotlib.pyplot as plt

def is_leap_year(year):
    return (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0))

def remove_high_frequent_varibility(var_name, file_in, file_out, window, year_s, year_e):

    # Read LAI input file
    f = Dataset(file_in, 'r')
    time = num2date(f.variables['time'][:], f.variables['time'].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    LAI = f.variables[var_name][:]

    # Read climatology files
    file_clim_leap = '/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_climatology/MCD15A3H.061_clim_5km_leapyear_31day_smooth_SG_filter_window=15timestep_order=2.nc'
    f_leap = Dataset(file_clim_leap, 'r')
    LAI_leap = f_leap.variables['LAI'][:]

    file_clim_common = '/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_climatology/MCD15A3H.061_clim_5km_commonyear_31day_smooth_SG_filter_window=15timestep_order=2.nc'
    f_common = Dataset(file_clim_common, 'r')
    LAI_common = f_common.variables['LAI'][:]

    # Calculate anomaly
    ts = 0
    for year in np.arange(year_s, year_e):
        if is_leap_year(year):
            LAI[ts:ts+366, :, :] -= LAI_leap[:, :, :]
            ts += 366
        else:
            LAI[ts:ts+365, :, :] -= LAI_common[:, :, :]
            ts += 365

    # Calculate smooth
    LAI_anomaly = xr.DataArray(
        LAI,
        dims=["time", "lat", "lon"],
        coords={"time": time, "lat": lat, "lon": lon},
        name="LAI"
    )

    # Calculate the 15-day rolling mean along the 'time' dimension
    LAI_anomaly_smooth = LAI_anomaly.rolling(time=window, center=True, min_periods=1).mean() #

    # Add climatology back
    ts = 0
    for year in np.arange(year_s, year_e):
        if is_leap_year(year):
            LAI_anomaly_smooth[ts:ts+366, :, :] += LAI_leap[:, :, :]
            ts += 366
        else:
            LAI_anomaly_smooth[ts:ts+365, :, :] += LAI_common[:, :, :]
            ts += 365

    LAI_anomaly_smooth

    # Plot the anomalies and smoothed anomalies
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot(LAI_anomaly[1000:1365, 238, 793], color='red', label='Anomaly 1')
    ax.plot(LAI_anomaly_smooth[:, 238, 793], color='blue', label='Smooth 1')
    ax.plot(LAI_anomaly[1000:1365, 131, 736], color='orange', label='Anomaly 2')
    ax.plot(LAI_anomaly_smooth[1000:1365, 131, 736], color='green', label='Smooth 2')

    ax.legend()
    plt.tight_layout()
    plt.savefig('./plots/check_smoothed_anomaly.png', dpi=100)

    # Write the smoothed data to a new NetCDF file
    if os.path.exists(file_out):
        os.remove(file_out)

    ds_out = Dataset(file_out, 'w')

    # Create dimensions
    ds_out.createDimension('time', len(time))
    ds_out.createDimension('latitude', len(lat))
    ds_out.createDimension('longitude', len(lon))

    # Create variables and set attributes
    time_out = ds_out.createVariable('time', 'f4', ('time',))
    lat_out = ds_out.createVariable('latitude', 'f4', ('latitude',))
    lon_out = ds_out.createVariable('longitude', 'f4', ('longitude',))
    lai_out = ds_out.createVariable('LAI', 'f4', ('time', 'latitude', 'longitude'))

    # Copy data into new variables
    time_out[:] = f.variables['time'][:]
    lat_out[:] = lat[:]
    lon_out[:] = lon[:]
    lai_out[:, :, :] = LAI_anomaly_smooth.values

    # Set attributes
    time_out.setncatts(f.variables['time'].__dict__)
    lat_out.setncatts(f.variables['latitude'].__dict__)
    lon_out.setncatts(f.variables['longitude'].__dict__)

    lai_attributes = f.variables[var_name].__dict__.copy()
    lai_attributes.pop('_FillValue', None)  # Remove '_FillValue' if present
    lai_out.setncatts(lai_attributes)
    lai_out.setncattr('long_name', 'LAI smoothed with 31-day rolling mean')

    print(f"Smoothed data saved to {file_out}")
    f.close()
    f_leap.close()
    f_common.close()
    ds_out.close()

if __name__ == "__main__":

    file_in  = "/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily/MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2000-2023_SG_filter_window=15timestep_order=2.nc"
    file_out = '/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily/remove_high_frequent_varibility_method2_smooth_anomaly/MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2000-2023_remove_high_freq_SG_filter_window=15timestep_order=2.nc'
    window   = 31
    year_s   = 2000
    year_e   = 2023
    var_name = 'LAI'
    remove_high_frequent_varibility(var_name, file_in, file_out, window, year_s, year_e)
