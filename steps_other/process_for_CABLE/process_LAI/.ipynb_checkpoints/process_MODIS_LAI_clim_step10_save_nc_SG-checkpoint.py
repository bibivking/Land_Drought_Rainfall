import os
from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import matplotlib.ticker as mticker
from multiprocessing import Pool, cpu_count
from common_utils import *

def days_of_year(year):
    if year%4==0:
        return 366
    else:
        return 365

def save_climatology_for_CABLE():
    
    '''
    Input: 
        MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2003-2023_climatology_common.nc
    Output: 
        MCD15A3H.061_2003-2023_climatology_5km_common_year.nc
        MCD15A3H.061_2003-2023_climatology_5km_leap_year.nc
    Function:
        1. Scaling the MODIS LAI by 0.1, and replace missing values to 0.001 (=C%LAI_THRESH in CABLE code)
        2. Produce new netcdf files 
        3. Note that climatology doesn't use the 11-day smoothed daily LAI 
    '''

    # ======== Read MODIS LAI ========
    scale_factor = 0.1
    file_path    = "/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_climatology/"

    # Scaling and reset default values
    file_in      = file_path+"MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2003-2023_climatology_common_SG_filter_window=15timestep_order=2.nc"
    f_in         = nc.Dataset(file_in,'r')
    lai_in       = f_in.variables['Lai_500m'][:,:,:] #*scale_factor
    lai_in       = np.where(lai_in>20., 0.001, lai_in) # set to the value of C%LAI_THRESH in CABLE
    lat_in       = f_in.variables['latitude'][:]
    lon_in       = f_in.variables['longitude'][:]
    nlat         = len(lat_in)
    nlon         = len(lon_in)
    f_in.close()
	
    # ======== saving common year ========
    # create file and write global attributes
    file_com            = file_path+"MCD15A3H.061_2003-2023_climatology_5km_common_year_SG_filter_window=15timestep_order=2.nc"
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

    var_com                = f_com.createVariable('LAI', 'f4', ('time', 'lat', 'lon'),fill_value=-9999.)
    var_com.long_name      = "MCD15A3H MODIS/Terra Gridded 500M Leaf Area Index LAI (4-day composite)"
    var_com.standard_name  = "LAI"
    var_com.units          = "m^2/m^2"
    var_com[:]             = lai_in

    f_com.close()

    # ======== saving leap year ========
    # create file and write global attributes
    file_leap            = file_path+"MCD15A3H.061_2003-2023_climatology_5km_leap_year_SG_filter_window=15timestep_order=2.nc"
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

    lai_leap_in             = np.zeros((366, nlat, nlon))
    lai_leap_in[:59,:,:]    = lai_in[:59,:,:]
    lai_leap_in[59,:,:]     = (lai_in[58,:,:]+lai_in[59,:,:])/2.
    lai_leap_in[60:366,:,:] = lai_in[59:365,:,:]

    var_leap                = f_leap.createVariable('LAI', 'f4', ('time', 'lat', 'lon'),fill_value=-9999.)
    var_leap.long_name      = "MCD15A3H MODIS/Terra Gridded 500M Leaf Area Index LAI (4-day composite)"
    var_leap.standard_name  = "LAI"
    var_leap.units          = "m^2/m^2"
    var_leap[:]             = lai_leap_in
    f_leap.close()

    return

if __name__ == "__main__":

    save_climatology_for_CABLE()