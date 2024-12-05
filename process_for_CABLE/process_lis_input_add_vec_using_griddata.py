'''
The problem with this file is that it runs too slow, 
so it is not practical to use for a lot of WRF simulations.
'''

import argparse
import gc
import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata
import logging
from multiprocessing import Pool
import time
import xesmf as xe
import xarray as xr

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger()

def regrid_data(lat_in, lon_in, lat_out, lon_out, input_data, method='linear', threshold=None):
    
    start = time.time()
    
    if len(np.shape(lat_in)) == 1:
        lon_in_2D, lat_in_2D = np.meshgrid(lon_in, lat_in)
    else:
        lon_in_2D, lat_in_2D = lon_in, lat_in

    lon_in_1D = lon_in_2D.ravel()
    lat_in_1D = lat_in_2D.ravel()
    value_tmp = input_data.ravel()

    # Mask invalid values
    if threshold is None:
        mask_values = ~np.isnan(value_tmp)
    else:
        mask_values = np.all([~np.isnan(value_tmp), value_tmp > threshold], axis=0)

    lon_in_1D, lat_in_1D, value = lon_in_1D[mask_values], lat_in_1D[mask_values], value_tmp[mask_values]

    # Interpolation
    if len(np.shape(lat_in)) == 1:
        lon_out_2D, lat_out_2D = np.meshgrid(lon_out, lat_out)
    else:
        lon_out_2D = lon_out
        lat_out_2D = lat_out

    Value = griddata((lon_in_1D, lat_in_1D), value, (lon_out_2D, lat_out_2D), method=method)
    gc.collect()

    end = time.time()
    print(f"Execution time of regrid_data: {end - start} seconds")
    
    return Value

def process_soil_layer(layer, gridinfo_file, gridinfo_var_name, lis_input_file, lis_var_name, mask_lis, lat_lis, lon_lis):

    start = time.time()
    
    # Read gridinfo file
    with nc.Dataset(gridinfo_file, 'r') as f_gridinfo:
        lat_grid = f_gridinfo.variables['lat'][:]
        lon_grid = f_gridinfo.variables['lon'][:]
        var_in   = f_gridinfo.variables[gridinfo_var_name][layer, :, :]

    # Regrid data for the soil layer
    regridded_data = regrid_data(lat_grid, lon_grid, lat_lis, lon_lis, var_in)
    regridded_data = np.where(mask_lis == 1, regridded_data, -9999.)
    gc.collect()

    end = time.time()
    print(f"Execution time of process_soil_layer: {end - start} seconds")

    return layer, regridded_data

def finalize_lis_file(lis_input_file, lis_var_name, soil_depth_data):
    
    start = time.time()
    lis_none_vec = {'SAND_VEC': 'SAND',
                    'CLAY_VEC': 'CLAY', 
                    'SILT_VEC': 'SILT',
                    'OCSOIL_VEC': 'OCSOIL',
                    'BULKDENSITY_VEC': 'BULKDENSITY' }
    
    with nc.Dataset(lis_input_file, 'r+') as f_lis_input:
        # Create soil_depth dimension if not present
        if 'soil_depth' not in f_lis_input.dimensions:
            f_lis_input.createDimension('soil_depth', soil_depth_data.shape[0])

        # Add variable
        var_out = f_lis_input.createVariable(lis_var_name, 'f8', ('soil_depth', 'north_south', 'east_west'), fill_value=-9999)
        var_out[:] = soil_depth_data

        # copy the attributes
        original_var = f_lis_input.variables[lis_none_vec[lis_var_name]]
        for attr_name in original_var.ncattrs():
            var_out.setncattr(attr_name, original_var.getncattr(attr_name))
    gc.collect()
    end = time.time()
    print(f"Execution time of finalize_lis_file: {end - start} seconds")

    return

def regrid_vec_to_wrf_domain_parallel(n, case_name, lis_input_path, gridinfo_file, gridinfo_var_names, lis_var_names):

    lis_input_file = f'{lis_input_path}/LIS_offline/lis_input.d0{n}.nc'

    # Read LIS grid information
    with nc.Dataset(lis_input_file, 'r') as f_lis_input:
        lat_lis  = f_lis_input.variables['lat'][:]
        lon_lis  = f_lis_input.variables['lon'][:]
        mask_lis = f_lis_input.variables['LANDMASK'][:]

    # Prepare arguments for each variable and soil layer
    args_list = []
    for i in range(len(gridinfo_var_names)):
        for layer in range(6):  # Assuming 6 soil layers
            args_list.append((layer, 
                              gridinfo_file, 
                              gridinfo_var_names[i], 
                              lis_input_file, 
                              lis_var_names[i], 
                              mask_lis, 
                              lat_lis, 
                              lon_lis))

    # Parallelize with a limited number of processes
    with Pool(processes=1) as pool:  # You can experiment with different numbers of processes
        results = pool.starmap(process_soil_layer, args_list)

    # Organize regridded data and finalize LIS file
    for i, var_name in enumerate(lis_var_names):
        soil_depth_data = np.zeros((6, lat_lis.shape[0], lon_lis.shape[1]))
        for layer, regridded_layer_data in results: 
            if var_name in args_list[layer]: 
                soil_depth_data[layer, :, :] = regridded_layer_data
        finalize_lis_file(lis_input_file, var_name, soil_depth_data)

    return 

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Interpolate _vec to lis_input.d0X.nc.')
    parser.add_argument('--case_name', type=str, default=None, help='What is the name of the case study')
    parser.add_argument('--n', type=int, default=None, help='The number of domain')

    # args = parser.parse_args()
    # n = args.n
    # case_name = args.case_name

    n                  = 1
    case_name          = 'test_runs'

    lis_var_names      = ['SAND_VEC', 'CLAY_VEC', 'SILT_VEC', 'OCSOIL_VEC', 'BULKDENSITY_VEC']
    gridinfo_var_names = ['SAND_VEC', 'CLAY_VEC', 'SILT_VEC', 'OC_VEC', 'BULK_DEN_VEC']

    lis_input_path = f'/scratch/w97/mm3972/model/NUWRF/Drought_breaking_rainfall/{case_name}'
    gridinfo_file  = '/g/data/w97/mm3972/scripts/wrf_scripts/make_LIS_landinfo/nc_file/Openlandmap_soilcomposition_CORDEX_180E_depth_varying.nc'

    regrid_vec_to_wrf_domain_parallel(n, case_name, lis_input_path, gridinfo_file, gridinfo_var_names, lis_var_names)
