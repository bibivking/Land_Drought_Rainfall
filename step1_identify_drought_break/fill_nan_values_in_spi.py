import os
import sys
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
from common_utils import *

def fill_nan_value_in_spi(var_name, file_in, file_out, file_mask):

    # Read SPI
    f_spi      = nc.Dataset(file_in, 'r')
    spi        = f_spi.variables[var_name][:,10:,:-45]
    lat_spi    = f_spi.variables['lat'][10:]
    lon_spi    = f_spi.variables['lon'][:-45]
    time       = f_spi.variables['time'][:]
    ntime      = len(spi[:,0,0])

    # Read mask
    f_mask     = nc.Dataset(file_mask, 'r')
    landsea    = f_mask.variables['landsea'][:]

    # Expand the array to 3D by adding a time axis at the beginning
    landsea_3d = np.expand_dims(landsea, axis=0)  # Shape becomes (1, lat, lon)

    # Repeat along the new time axis
    landsea_3d = np.repeat(landsea_3d, ntime, axis=0)  # Shape becomes (time, lat, lon)

    # mask out sea pixel
    spi_tmp    = np.where(landsea_3d == 1, -9999., spi)

    # set nan value as 0
    spi_tmp    = np.where(np.isnan(spi_tmp), 0, spi_tmp)

    # delete output file if it exists
    if os.path.exists(file_out):
        os.remove(file_out)

    # Write the smoothed data to a new NetCDF file
    with Dataset(file_out, 'w') as ds_out:
        # Create dimensions
        ds_out.createDimension('time', ntime)
        ds_out.createDimension('lat', len(lat_spi))
        ds_out.createDimension('lon', len(lon_spi))

        # Create variables and set attributes
        time_out = ds_out.createVariable('time', 'f4', ('time',))
        lat_out  = ds_out.createVariable('lat', 'f4', ('lat',))
        lon_out  = ds_out.createVariable('lon', 'f4', ('lon',))
        spi_out  = ds_out.createVariable(var_name, 'f4', ('time', 'lat', 'lon'), fill_value=-9999.)

        # Copy data into new variables
        time_out[:]      = time
        lat_out[:]       = lat_spi
        lon_out[:]       = lon_spi
        spi_out[:, :, :] = spi_tmp

        # Copy attributes if needed
        time_attributes = f_spi.variables['time'].__dict__.copy()
        time_attributes.pop('_FillValue', None)  # Remove '_FillValue' if present
        time_out.setncatts(time_attributes)

        lat_attributes  = f_spi.variables['lat'].__dict__.copy()
        lat_attributes.pop('_FillValue', None)  # Remove '_FillValue' if present
        lat_out.setncatts(lat_attributes)

        lon_attributes  = f_spi.variables['lon'].__dict__.copy()
        lon_attributes.pop('_FillValue', None)  # Remove '_FillValue' if present
        lon_out.setncatts(lon_attributes)

        spi_attributes = f_spi.variables[var_name].__dict__.copy()
        spi_attributes.pop('_FillValue', None)  # Remove '_FillValue' if present
        spi_out.setncatts(spi_attributes)
        spi_out.setncattr('long_name', f'{var_name} smoothed with 31-day cubic spline')

    return

if __name__ == "__main__":

    file_in     = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_90_reorder.nc'
    file_out    = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_90_reorder_nan_filled.nc'
    # file_spi_90 = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_90_reorder.nc'
    file_mask   = '/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_CSIRO_AU_NAT_ELEV_DLCM_mask.nc'
    var_name    = 'spi_pearson_90'

    fill_nan_value_in_spi(var_name, file_in, file_out, file_mask)
