import argparse
import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata

def is_leap(year):
    """Check if a year is a leap year."""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)

def regrid_data(lat_in, lon_in, lat_out, lon_out, input_data, method='linear',threshold=None):

    if len(np.shape(lat_in)) == 1:
        lon_in_2D, lat_in_2D = np.meshgrid(lon_in,lat_in)
        lon_in_1D            = np.reshape(lon_in_2D,-1)
        lat_in_1D            = np.reshape(lat_in_2D,-1)
    elif len(np.shape(lat_in)) == 2:
        lon_in_1D            = np.reshape(lon_in,-1)
        lat_in_1D            = np.reshape(lat_in,-1)
    else:
        print("ERROR: lon_in has ", len(np.shape(lat_in)), "dimensions")

    if len(np.shape(lat_out)) == 1:
        lon_out_2D, lat_out_2D = np.meshgrid(lon_out,lat_out)
    elif len(np.shape(lat_out)) == 2:
        lon_out_2D            = lon_out
        lat_out_2D            = lat_out
    else:
        print("ERROR: lon_out has ", len(np.shape(lat_in)), "dimensions")

    value_tmp = np.reshape(input_data,-1)

    # Check NaN - input array shouldn't have NaN
    if threshold is None:
        mask_values     = ~np.isnan(value_tmp)
    else:
        #print("np.all([~np.isnan(value_tmp),value_tmp>threshold],axis=0) ",np.all([~np.isnan(value_tmp),value_tmp>threshold],axis=0))
        mask_values     = np.all([~np.isnan(value_tmp),value_tmp>threshold],axis=0)

    value     = value_tmp[mask_values]
    # ======= CAUTION =======
    lat_in_1D = lat_in_1D[mask_values]  # here I make nan in values as the standard
    lon_in_1D = lon_in_1D[mask_values]
    #print("shape value = ", np.shape(value))
    #print("shape lat_in_1D = ", np.shape(lat_in_1D))
    #print("shape lon_in_1D = ", np.shape(lon_in_1D))
    # =======================
    #print("value =",value)
    #print("lat_in_1D =",lat_in_1D)
    #print("lon_in_1D =",lon_in_1D)
    Value = griddata((lon_in_1D, lat_in_1D), value, (lon_out_2D, lat_out_2D), method=method)

    return Value

def use_SM_ST_after_90year_spinup(leap_year=None, year=None):

    file_in = "/g/data/w97/mm3972/model/cable/runs/Land_drought_rainfall_runs/prepare_SM_ST_for_offline_run/outputs/cable_out_1999.nc"

    with nc.Dataset(file_in, 'r') as f_in:
        SM_in  = f_in.variables['SoilMoist'][:,:, :, :]
        ST_in  = f_in.variables['SoilTemp'][:,:, :, :]
        lat_in = f_in.variables['latitude'][:, :]
        lon_in = f_in.variables['longitude'][:, :]
        SM_in  = np.where(SM_in<0, np.nan, SM_in)
        ST_in  = np.where(ST_in<0, np.nan, ST_in)

    # Open the original NetCDF file
    if leap_year==None:
        file_out = f'/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_time_varying_{year}.nc'
    elif leap_year=='True':
        file_out = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_clim_leap.nc'
    elif leap_year=='False':
        file_out = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_clim_common.nc'

    # Open the source NetCDF file in read mode
    f_out    = nc.Dataset(file_out, 'r+')
    lat_out  = f_out.variables['latitude'][:]
    lon_out  = f_out.variables['longitude'][:]
    watr_out = f_out.variables['watr'][:]

    for m in np.arange(12):
        for s in np.arange(6):
            SoilMoist_tmp = regrid_data(lat_in, lon_in, lat_out, lon_out, SM_in[m,s,:,:], method='nearest',threshold=None)
            f_out.variables['SoilMoist'][m,s,:,:] = np.where(SoilMoist_tmp > watr_out[s,:,:], SoilMoist_tmp, watr_out[s,:,:]+0.001) # to make sure SM is always larger than watr
            f_out.variables['SoilTemp'][m,s,:,:]  = regrid_data(lat_in, lon_in, lat_out, lon_out, ST_in[m,s,:,:], method='nearest',threshold=None)

    f_out.close()

    return

def increase_hydraulic_conductivity(leap_year=None, year=None):

    # Open the original NetCDF file
    if leap_year==None:
        file_out = f'/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_time_varying_{year}.nc'
    elif leap_year=='True':
        file_out = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_clim_leap.nc'
    elif leap_year=='False':
        file_out = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_clim_common.nc'

    # Open the source NetCDF file in read mode
    f_out    = nc.Dataset(file_out, 'r+')
    f_out.variables['hyds'][:]      = f_out.variables['hyds'][:]*10.
    f_out.variables['hyds_vec'][:]  = f_out.variables['hyds_vec'][:]*10.

    f_out.close()

    return

def reduce_soil_moisture_in_gridinfo(reduce_percent=10, leap_year=None, year=None):

    # Shink to
    shink_to = (100-reduce_percent)/100

    # Reduce soil column moisture
    # Open the original NetCDF file
    if leap_year==None:
        file_out = f'/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/reduce_{reduce_percent}percent_SM/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_time_varying_{year}.nc'
    elif leap_year=='True':
        file_out = f'/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/reduce_{reduce_percent}percent_SM/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_clim_leap.nc'
    elif leap_year=='False':
        file_out = f'/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/reduce_{reduce_percent}percent_SM/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_clim_common.nc'

    with nc.Dataset(file_out, 'r+') as f_out:
        # Open the source NetCDF file in read mode
        watr_out = f_out.variables['watr'][:]
        for m in np.arange(12):
            for s in np.arange(6):
                SoilMoist_tmp = f_out.variables['SoilMoist'][m,s,:,:].data*shink_to
                f_out.variables['SoilMoist'][m,s,:,:] = np.where(SoilMoist_tmp > watr_out[s,:,:], SoilMoist_tmp, watr_out[s,:,:]+0.001) # to make sure SM is always larger than watr

    return

def reduce_soil_moisture_in_restart(reduce_percent=10):
   
    # Shink to
    shink_to = (100-reduce_percent)/100

    # Reduce aquifer moisture
    file_restart = f'/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/reduce_{reduce_percent}percent_SM/restart_1970_one_year_5km_run_with_SM_ST_spinup_read_from_gridinfo.nc'
    
    with nc.Dataset(file_restart, 'r+') as f_restart:
        f_restart.variables['wbtot0'][:] = f_restart.variables['wbtot0'][:]*shink_to 
        f_restart.variables['GWwb'][:]   = f_restart.variables['GWwb'][:]*shink_to
        
        watr_hys                         = f_restart.variables['watr_hys'][:]

        tmp                              = f_restart.variables['wb'][:]*shink_to  
        f_restart.variables['wb'][:]     = np.where(tmp>watr_hys, tmp, watr_hys)
        
        tmp                              = f_restart.variables['wb_hys'][:]*shink_to   
        f_restart.variables['wb_hys'][:] = np.where(tmp>watr_hys, tmp, watr_hys)

    return

if __name__ == "__main__":

    reduce_percent = 10 

    # leap_year = "True"
    # reduce_soil_moisture_in_gridinfo(reduce_percent=reduce_percent, leap_year=leap_year)
    # print('Finish leap year')

    # leap_year = "False"
    # reduce_soil_moisture_in_gridinfo(reduce_percent=reduce_percent, leap_year=leap_year)
    # print('Finish common year')

    # for year in np.arange(2000,2024,1):
    #     reduce_soil_moisture_in_gridinfo(reduce_percent=reduce_percent, year=year)
    #     print(f'Finish {year}')

    reduce_soil_moisture_in_restart(reduce_percent)
    print('Finish restart')