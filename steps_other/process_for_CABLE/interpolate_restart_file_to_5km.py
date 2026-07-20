import argparse
import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata
from multiprocess import Pool
from common_utils import *

def interpolate_to_5km( value, lat_in, lon_in, lat_out_grid, lon_out_grid, lat_out, lon_out, method='nearest'):

    value_grid = regrid_data(lat_in, lon_in, lat_out_grid, lon_out_grid, value, method)
    mp         = len(lat_out)
    value_out  = np.zeros(mp)

    # Set the output lat and lon
    lon_out_2D, lat_out_2D = np.meshgrid(lon_out_grid, lat_out_grid)

    for i in np.arange(mp):
        mask         = ( abs(lat_out_2D - lat_out[i]) < 0.0001) & ( abs(lon_out_2D - lon_out[i]) < 0.0001)

        if np.any(mask == True):
            value_out[i] = value_grid[mask][0]
        else:
            print('i, lat_out[i], lon_out[i]', i, lat_out[i], lon_out[i], 'does not exist')

    return value_out

def extract_value(args):
    i, lat_out, lon_out, lat_out_2D, lon_out_2D, value_grid = args

    # Create the mask for the current lat/lon pair
    mask = (np.abs(lat_out_2D - lat_out[i]) < 0.0001) & (np.abs(lon_out_2D - lon_out[i]) < 0.0001)

    if np.any(mask):
        return value_grid[mask][0]  # Return the first matching value
    else:
        print(f'i, lat_out[i], lon_out[i]: {i}, {lat_out[i]}, {lon_out[i]} does not exist')
        return np.nan  # Return NaN if not found

def interpolate_to_5km_parallel(value, lat_in, lon_in, lat_out_grid, lon_out_grid, lat_out, lon_out, method='nearest'):

    value_grid = regrid_data(lat_in, lon_in, lat_out_grid, lon_out_grid, value, method)

    # Set the output lat and lon
    lon_out_2D, lat_out_2D = np.meshgrid(lon_out_grid, lat_out_grid)
    mp = len(lat_out)

    # Prepare arguments for parallel processing
    args = [(i, lat_out, lon_out, lat_out_2D, lon_out_2D, value_grid) for i in range(mp)]

    # Use multiprocess Pool for parallel processing
    with Pool() as pool:
        value_out = pool.map(extract_value, args)

    return np.array(value_out)

if __name__ == "__main__":

    file_10km      = '/g/data/w97/mm3972/model/cable/runs/Land_drought_rainfall_runs/prepare_SM_ST_for_offline_run/restarts/restart_1999.nc'
    file_5km_grid  = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_clim_common.nc'
    file_5km       = '/g/data/w97/mm3972/model/cable/runs/Land_drought_rainfall_runs/spinup_run_1970_1999/restarts/restart_1970_one_year_5km_run_with_SM_ST_spinup_read_from_gridinfo.nc'

    soil_carbon_pools  = 2
    plant_carbon_pools = 3
    nsoil              = 6

    var_names_1d_mp = [ 'wbtot0', 'GWwb']
    # 'trad', 'tss', 'rtsoil', 'runoff', 'rnof1', 'rnof2', 'ghflux', 'sghflux', 'ga',
    # 'sublayer_dz', 'dgdtg', 'fev', 'fes', 'fhs',

    var_names_2d_soil_mp = ['wb', 'wbice', 'wb_hys', 'smp_hys', 'hys_fac']
    # 'tgg', 'tggsn', 'wb_hys',
    # 'smp_hys', 'ssat_hys', 'watr_hys', 'hys_fac'

    f_5km    = nc.Dataset(file_5km, 'r+')
    lat_5km  = f_5km.variables['latitude'][:]
    lon_5km  = f_5km.variables['longitude'][:]

    f_5km_grid    = nc.Dataset(file_5km_grid, 'r')
    lat_5km_grid  = f_5km_grid.variables['latitude'][:]
    lon_5km_grid  = f_5km_grid.variables['longitude'][:]

    f_10km   = nc.Dataset(file_10km, 'r')
    lat_10km = f_10km.variables['latitude'][:]
    lon_10km = f_10km.variables['longitude'][:]

    for var_name in var_names_1d_mp:
        print('processing',var_name)
        value_10km = f_10km.variables[var_name][:].data
        f_5km.variables[var_name][:] = interpolate_to_5km_parallel(value_10km, lat_10km, lon_10km, lat_5km_grid, lon_5km_grid, lat_5km, lon_5km, )

    for var_name in var_names_2d_soil_mp:
        print('processing',var_name)
        for s in np.arange(nsoil):
            value_10km = f_10km.variables[var_name][:].data  # Add this line to retrieve the correct value for each soil layer
            f_5km.variables[var_name][s,:] = interpolate_to_5km_parallel(value_10km[s,:], lat_10km, lon_10km, lat_5km_grid, lon_5km_grid, lat_5km, lon_5km, )

    cplant_10km = f_10km.variables['cplant'][:].data
    for p in np.arange(plant_carbon_pools):
        print('processing cplant')
        f_5km.variables['cplant'][p,:] = interpolate_to_5km_parallel(cplant_10km[p,:], lat_10km, lon_10km, lat_5km_grid, lon_5km_grid, lat_5km, lon_5km,)

    csoil_10km = f_10km.variables['csoil'][:].data
    for p in np.arange(soil_carbon_pools):
        print('processing csoil')
        f_5km.variables['csoil'][p,:] = interpolate_to_5km_parallel(csoil_10km[p,:], lat_10km, lon_10km, lat_5km_grid, lon_5km_grid, lat_5km, lon_5km, )

    f_5km.close()
    f_5km_grid.close()
    f_10km.close()
