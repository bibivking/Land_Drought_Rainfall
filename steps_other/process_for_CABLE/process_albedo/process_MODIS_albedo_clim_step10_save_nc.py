import os
import gc
from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
import argparse
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import matplotlib.ticker as mticker
from multiprocessing import Pool, cpu_count
from common_utils import *
from process_MODIS_albedo_time_varying_step8_gap_fill import gap_fill, fill_gap_for_point

def days_of_year(year):
    if year%4==0:
        return 366
    else:
        return 365

def save_climatology_common_for_CABLE(albedo_band):

    '''
    Input:
        MCD15A3H.061_500m_aid0001_albedo_regridded_daily_2003-2023_climatology_common.nc
    Output:
        MCD15A3H.061_2003-2023_climatology_5km_common_year.nc
        MCD15A3H.061_2003-2023_climatology_5km_leap_year.nc
    Function:
        1. Scaling the MODIS albedo by 0.1, and replace missing values to 0.001 (=C%albedo_THRESH in CABLE code)
        2. Produce new netcdf files
        3. Note that climatology doesn't use the 11-day smoothed daily albedo
    '''

    # ======== Read MODIS albedo ========
    file_path    = "/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_climatology/"
    var_name     = "Albedo_"+albedo_band
    # Scaling and reset default values
    file_in      = file_path+f"/MCD43A3.061_500m_aid0001_{albedo_band}_2003-2023_climatology_common.nc"
    f_in         = nc.Dataset(file_in,'r')
    albedo_in    = f_in.variables[var_name][:,:,:]

    # print('np.unique(f_in.variables[var_name][:,:,:])',np.unique(f_in.variables[var_name][:,:,:]))
    # print('np.unique(f_in.variables[var_name][:,:,:].data)',np.unique(f_in.variables[var_name][:,:,:].data))
    lat_in       = f_in.variables['latitude'][:]
    lon_in       = f_in.variables['longitude'][:]
    nlat         = len(lat_in)
    nlon         = len(lon_in)

    f_in.close()

    # ======== saving common year ========
    # create file and write global attributes
    file_com            = file_path+f"MCD43A3.061_2003-2023_{albedo_band}_climatology_common.nc"
    f_com               = nc.Dataset(file_com, 'w', format='NETCDF4')
    f_com.history       = "Created by: %s" % (os.path.basename(__file__))
    f_com.creation_date = "%s" % (datetime.now())

    # set dimensions
    f_com.createDimension('time', None)
    f_com.createDimension('lat', nlat)
    f_com.createDimension('lon', nlon)
    f_com.Conventions      = "CF-1.0"

    time_com               = f_com.createVariable('time', 'f8', ('time'))
    time_com.long_name     = "Time"
    time_com.standard_name = "time"
    time_com.calendar      = "proleptic_gregorian"
    time_com.units         = "days since 2001-01-01 00:00:00"
    time_com[:]            = np.arange(0, 365, 1)

    lat_com                = f_com.createVariable('lat', 'f4', ('lat'))
    lat_com.long_name      = "Latitude"
    lat_com.standard_name  = "latitude"
    lat_com.axis           = "Y"
    lat_com.units          = "degrees_north"
    lat_com[:]             = lat_in

    lon_com                = f_com.createVariable('lon', 'f4', ('lon'))
    lon_com.long_name      = "Longitude"
    lon_com.standard_name  = "longitude"
    lon_com.axis           = "X"
    lon_com.units          = "degrees_east"
    lon_com[:]             = lon_in

    var_com                = f_com.createVariable(var_name, 'i2', ('time', 'lat', 'lon'),fill_value=32767)
    var_com.long_name      = f"Albedo_{albedo_band}"
    var_com.standard_name  = f"Albedo_{albedo_band}"
    var_com.add_offset     = 0
    var_com.units          = "reflectance, no units"
    var_com.scale_factor   = 0.001
    var_com.missing_value  = 32767
    var_com[:]             = albedo_in

    # print('np.unique(f_com.variables[var_name][:,:,:])',np.unique(f_com.variables[var_name][:,:,:]))
    # print('np.unique(f_com.variables[var_name][:,:,:].data)',np.unique(f_com.variables[var_name][:,:,:].data))

    f_com.close()

    return

def gap_fill_climatology_common(albedo_band, window):

    '''
    Gap fill common year and make leap year
    '''
    file_path    = "/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_climatology/"
    
    scale_factor   = 0.001
    missing_value  = 32767
    var_name       = "Albedo_"+albedo_band

    f_gap_fill     = nc.Dataset(file_path+f"MCD43A3.061_2003-2023_{albedo_band}_climatology_common.nc", 'r+')
    var            = f_gap_fill.variables[var_name][:,:,:]
    lat_in         = f_gap_fill.variables['lat'][:]
    lon_in         = f_gap_fill.variables['lon'][:]

    ntime          = var.shape[0]
    nlat           = var.shape[1]
    nlon           = var.shape[2]

    var_in             = np.zeros((ntime+31*2,nlat,nlon))
    var_in[:31,:,:]    = var[-31:,:,:]
    var_in[31:-31,:,:] = var[:,:,:]
    var_in[-31:,:,:]   = var[:31,:,:]
    var_in             = np.where(var_in > 1, np.nan, var_in)
    var_in             = np.where(var_in < 0, np.nan, var_in)

    var_gap_filled     = gap_fill(var_in, window)

    var_out            = var_gap_filled[31:-31,:,:]

    gc.collect()

    var_out            = np.where(np.isnan(var_out), missing_value*scale_factor, var_out)
    var_out            = np.where(var_out>1, missing_value*scale_factor, var_out)
    var_out            = np.where(var_out<0, missing_value*scale_factor, var_out)

    f_gap_fill.variables[var_name][:,:,:] = var_out

    f_gap_fill.close()

def save_climatology_leap_for_CABLE(albedo_band):

    # 
    scale_factor   = 0.001
    missing_value  = 32767

    # READ IN
    file_path      = "/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_climatology/"
    var_name       = "Albedo_"+albedo_band

    f_com          = nc.Dataset(file_path+f"MCD43A3.061_2003-2023_{albedo_band}_climatology_common.nc", 'r')
    var_com        = f_com.variables[var_name][:,:,:] # when add .data, the default value will become -1, I don't know why
    lat_in         = f_com.variables['lat'][:]
    lon_in         = f_com.variables['lon'][:]
    f_com.close()

    ntime          = var_com.shape[0]
    nlat           = var_com.shape[1]
    nlon           = var_com.shape[2]

    # ======== saving leap year ========
    # create file and write global attributes
    file_leap            = file_path+f"MCD43A3.061_2003-2023_{albedo_band}_climatology_leap.nc"
    f_leap               = nc.Dataset(file_leap, 'w', format='NETCDF4')
    f_leap.history       = "Created by: %s" % (os.path.basename(__file__))
    f_leap.creation_date = "%s" % (datetime.now())

    # set dimensions
    f_leap.createDimension('time', None)
    f_leap.createDimension('lat', nlat)
    f_leap.createDimension('lon', nlon)
    f_leap.Conventions   = "CF-1.0"

    time_leap               = f_leap.createVariable('time', 'f8', ('time'))
    time_leap.long_name     = "Time"
    time_leap.standard_name = "time"
    time_leap.calendar      = "proleptic_gregorian"
    time_leap.units         = "days since 2000-01-01 00:00:00"
    time_leap[:]            = np.arange(0, 366, 1)

    lat_leap                = f_leap.createVariable('lat', 'f4', ('lat'))
    lat_leap.long_name      = "Latitude"
    lat_leap.standard_name  = "latitude"
    lat_leap.axis           = "Y"
    lat_leap.units          = "degrees_north"
    lat_leap[:]             = lat_in

    lon_leap                = f_leap.createVariable('lon', 'f4', ('lon'))
    lon_leap.long_name      = "Longitude"
    lon_leap.standard_name  = "longitude"
    lon_leap.axis           = "X"
    lon_leap.units          = "degrees_east"
    lon_leap[:]             = lon_in

    albedo_leap_in             = np.zeros((366, nlat, nlon), dtype=var_com.dtype) # defaults is float64, so need  dtype=var_com.dtype
    albedo_leap_in[:59,:,:]    = var_com[:59,:,:]
    albedo_leap_in[59,:,:]     = np.mean(var_com[58:60,:,:],axis=0)
    albedo_leap_in[60:366,:,:] = var_com[59:365,:,:]

    # if don't add thsee two lines, the pixel with missing value will be -1
    albedo_leap_in             = np.where(albedo_leap_in>1, missing_value*scale_factor, albedo_leap_in)
    albedo_leap_in             = np.where(albedo_leap_in<0, missing_value*scale_factor, albedo_leap_in)

    var_leap                = f_leap.createVariable(var_name, 'i2', ('time', 'lat', 'lon'),fill_value=32767)
    var_leap.long_name      = f"Albedo_{albedo_band}"
    var_leap.standard_name  = f"Albedo_{albedo_band}"
    var_leap.add_offset     = 0
    var_leap.units          = "reflectance, no units"
    var_leap.scale_factor   = 0.001
    var_leap.missing_value  = 32767
    var_leap[:]             = albedo_leap_in

    f_leap.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Save climatology to nc file for CABLE')
    parser.add_argument('--albedo_band', type=str, default="BSA_nir", help='Albedo band to process, BSA_nir, WSA_nir, BSA_vis, WSA_vis')
    parser.add_argument('--window', type=int, default=0, help='Window size for gap filling (default: 90 days)')

    args = parser.parse_args()

    # save_climatology_common_for_CABLE(args.albedo_band)
    # gap_fill_climatology_common(args.albedo_band, args.window)
    save_climatology_leap_for_CABLE(args.albedo_band)
