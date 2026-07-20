import os
import gc
import argparse
import copy
import logging
import netCDF4 as nc
import numpy as np
from multiprocessing import Pool, cpu_count
from common_utils import *

def is_leap(year):
    if year%4==0:
        return True
    else:
        return False

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Gap-fill MODIS albedo data by climatology.')
    parser.add_argument('file_path', type=str, help='Path to the netCDF file')
    parser.add_argument('--albedo_band', type=str, default="BSA_nir", help='Albedo band to process, BSA_nir, WSA_nir, BSA_vis, WSA_vis')
    
    args = parser.parse_args()

    # ==== Here is a mistake ====
    # The original default value for albedo is 255, but I changed it to LAI's in 
    # producing clim and time-varying gap-filled files. So keep default value as LAI
    scale_factor   = 0.001
    missing_value  = 32767
    # ===========================

    albedo_band    = args.albedo_band
    var_name       = "Albedo_"+albedo_band
    file_path      = args.file_path 

    # Get year
    split_path     = file_path.split('_')[-5]
    year           = int(split_path[-4:])
    
    if year == 2024:
        year = 2023

    print(split_path, year)

    # Read climatology
    if is_leap(year):
        clim_in_file = f"/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_climatology/MCD43A3.061_clim_{albedo_band}_leap_31day_smooth.nc"
    else:             
        clim_in_file = f"/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_climatology/MCD43A3.061_clim_{albedo_band}_common_31day_smooth.nc"
    f_in             = nc.Dataset(clim_in_file, 'r')
    var_clim_tmp     = f_in.variables[var_name][:,:,:] # do not need to time scale_factor, since f_gap_fill.variables by default times scale_factor
    var_clim_tmp     = np.where(((var_clim_tmp>1.) | (var_clim_tmp<0.)), np.nan, var_clim_tmp)
    f_in.close()

    # Read nearby gap filled file
    f_gap_fill     = nc.Dataset(file_path, 'r+')
    var            = f_gap_fill.variables[var_name][:,:,:] # do not need to time scale_factor, since f_gap_fill.variables by default times scale_factor
    var            = np.where(((var>1.) | (var<0.)), np.nan, var)
    ntime, nlat, nlon= var.shape[0], var.shape[1], var.shape[2]

    # Connect climatology input to the same as file needs filling 
    var_clim       = np.zeros((ntime, nlat, nlon))

    if year == 2000:
        # 24 Feb 2000 ~ 31 Dec 2000
        var_clim[:,:,:]   = var_clim_tmp[54:,:,:]
    elif year == 2023:
        # 1 Jan 2023 ~ 31 Jan 2024
        var_clim[:365,:,:]= var_clim_tmp
        var_clim[365:,:,:]= var_clim_tmp[:31,:,:]
    else:
        var_clim = var_clim_tmp

    print('np.sum(np.isnan(var))',np.sum(np.isnan(var)))
    var_out  = np.where(np.isnan(var), var_clim, var)
    print('np.sum(np.isnan(var_out))',np.sum(np.isnan(var_out)))
    var_out  = np.where(np.isnan(var_out), missing_value*scale_factor, var_out)
    f_gap_fill.variables[var_name][:,:,:] = var_out
    # print("np.unique(var_out[0,:,:])",np.unique(var_out[0,:,:]))

    f_gap_fill.close()
