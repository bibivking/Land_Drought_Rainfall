import os
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import argparse
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import matplotlib.ticker as mticker
from common_utils import *

def days_of_year(year):
    if year%4==0:
        return 366
    else:
        return 365

def save_daily_varying_for_CABLE(window):

    '''
    Input:
        MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2000-2023_xxday_smooth.nc
    Output:
        MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2000/.../2023_xxday_smooth.nc

    Function:
        1. Scaling the MODIS LAI by 0.1, and replace missing values to 0.001 (=C%LAI_THRESH in CABLE code)
        2. Produce new netcdf files
    '''

    scale_factor = 0.1
    file_path    = "/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily/remove_high_frequent_varibility_method2_smooth_anomaly/"

    # Scaling and reset default values
    file_in      = '/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily/remove_high_frequent_varibility_method2_smooth_anomaly/MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2000-2023_remove_high_freq.nc'
    f_in         = nc.Dataset(file_in,'r')
    lai_in       = f_in.variables['LAI'][:,:,:] #*scale_factor
    # lai_in       = np.where(lai_in>20., 0.001, lai_in) # set to the value of C%LAI_THRESH in CABLE
    lat_in       = f_in.variables['latitude'][:]
    lon_in       = f_in.variables['longitude'][:]
    nlat         = len(lat_in)
    nlon         = len(lon_in)
    f_in.close()

    yr_s         = 0
    yr_e         = 0

    # ======== saving daily varying ========
    # create file and write global attributes
    for yr in np.arange(2000, 2024, 1):

        yr_s = yr_e
        yr_e = yr_s+days_of_year(yr)

        file_out            = f"{file_path}MCD15A3H.061_500m_aid0001_LAI_regridded_daily_{str(yr)}_{window}day_smooth.nc"
        f_out               = nc.Dataset(file_out, 'w', format='NETCDF4')
        f_out.history       = "Created by: %s" % (os.path.basename(__file__))
        f_out.creation_date = "%s" % (datetime.now())

        # set dimensions
        f_out.createDimension('time', None)
        f_out.createDimension('lat', nlat)
        f_out.createDimension('lon', nlon)
        f_out.Conventions      = "CF-1.0"

        time_out               = f_out.createVariable('time', 'f8', ('time'))
        time_out.long_name     = "Time"
        time_out.standard_name = "time"
        time_out.calendar      = "proleptic_gregorian"
        time_out.units         = "days since 2000-01-01 00:00:00"
        time_out[:]            = np.arange(yr_s, yr_e, 1)

        lat_out                = f_out.createVariable('lat', 'f4', ('lat'))
        lat_out.long_name      = "Latitude"
        lat_out.standard_name  = "latitude"
        lat_out.axis           = "Y"
        lat_out.units          = "degrees_north"
        lat_out[:]             = lat_in

        lon_out                = f_out.createVariable('lon', 'f4', ('lon'))
        lon_out.long_name      = "Longitude"
        lon_out.standard_name  = "longitude"
        lon_out.axis           = "X"
        lon_out.units          = "degrees_east"
        lon_out[:]             = lon_in

        var_out                = f_out.createVariable('LAI', 'f4', ('time', 'lat', 'lon'),fill_value=-9999.)
        var_out.long_name      = "MCD15A3H MODIS/Terra Gridded 500M Leaf Area Index LAI (4-day composite)"
        var_out.standard_name  = "LAI"
        var_out.units          = "m^2/m^2"
        var_out[:]             = lai_in[yr_s:yr_e,:,:]

        f_out.close()

        time_out= None
        lat_out = None
        lon_out = None
        var_out = None

    return

if __name__ == "__main__":

    window   = 31
    save_daily_varying_for_CABLE(window)
