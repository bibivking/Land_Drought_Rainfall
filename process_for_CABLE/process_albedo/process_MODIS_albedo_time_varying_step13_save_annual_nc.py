import os
import argparse
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
import matplotlib.ticker as mticker
from multiprocessing import Pool, cpu_count
from common_utils import *

def days_of_year(year):
    if year%4==0:
        return 366
    else:
        return 365

def save_daily_varying_for_CABLE(albedo_band, window):
    
    '''
    Input: 
        e.g. MCD43A3.061_500m_aid0001_BSA_nir_Feb2000-Jan2024_albedo_regridded_daily_gapfill_31day_smooth.nc
    Output: 
        e.g. MCD43A3.061_albedo_BSA_nir_5km_2000...2003_31day_smooth.nc
        
    Function:
        1. Gap fill 2000 by climatology
        2. Produce annually netcdf files 
    '''

    # Setting
    var_name     = "Albedo_"+albedo_band 
    file_path    = "/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/"

    # Read input 
    file_in      = f"{file_path}MCD43A3.061_500m_aid0001_{albedo_band}_Feb2000-Jan2024_albedo_regridded_daily_gapfill_{window}day_smooth.nc"
    f_in         = nc.Dataset(file_in,'r')
    albedo_in    = f_in.variables[var_name][:,:,:]
    lat_in       = f_in.variables['latitude'][:]
    lon_in       = f_in.variables['longitude'][:]
    nlat         = len(lat_in)
    nlon         = len(lon_in)
    f_in.close()

    # Reset default values
    albedo_in    = np.where(albedo_in>1., 0, albedo_in) 
    albedo_in    = np.where(albedo_in<0., 0, albedo_in) 

    # Gap fill 2000 by climatology 
    day_out_tot = 0
    for yr in np.arange(2000, 2024, 1):
        day_out_tot = day_out_tot + days_of_year(yr)
    albedo_out = np.zeros((day_out_tot, nlat, nlon))

    # Read climatology albedo to fill 2000
    file_clim    = f"/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_climatology/MCD43A3.061_clim_{albedo_band}_leap_{window}day_smooth.nc"
    f_clim       = nc.Dataset(file_clim,'r')
    albedo_clim  = f_clim.variables[var_name][:,:,:]
    f_clim.close()
    
    # Use climatology fill from 1st Jan to 9th March 2000
    albedo_out[:69,:,:] = albedo_clim[:69,:,:]
    # Read time-varying from 10th March 2020 to 31st Dec 2023
    albedo_out[69:,:,:] = albedo_in[:-16,:,:]

    # ======== saving daily varying ========
    # create file and write global attributes
    yr_s         = 0
    yr_e         = 0
    for yr in np.arange(2000, 2024, 1):
        yr_s = yr_e
        yr_e = yr_s+days_of_year(yr)
            
        file_out            = f"{file_path}MCD43A3.061_albedo_{albedo_band}_5km_{yr}_{window}day_smooth.nc"
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

        var_out                = f_out.createVariable(var_name, 'f4', ('time', 'lat', 'lon'),fill_value=-9999.)
        var_out.long_name      = f"MCD43A3.061 for aid0001, {albedo_band}"
        var_out.standard_name  = var_name
        var_out.units          = "reflectance, no units"
        var_out[:]             = albedo_out[yr_s:yr_e,:,:]

        f_out.close()

        time_out= None
        lat_out = None
        lon_out = None
        var_out = None

    return

if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Save time-varying albedo to nc file for CABLE')
    parser.add_argument('--albedo_band', type=str, default="BSA_nir", help='Albedo band to process, BSA_nir, WSA_nir, BSA_vis, WSA_vis')
    parser.add_argument('--window', type=int, default=0, help='Window size for gap filling')

    args = parser.parse_args()

    save_daily_varying_for_CABLE(args.albedo_band, args.window)
