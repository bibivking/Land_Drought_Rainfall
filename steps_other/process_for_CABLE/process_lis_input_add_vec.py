import gc
import time
import os
import argparse
import logging
import netCDF4 as nc
import numpy as np
import xesmf as xe
import xarray as xr
import shutil
from multiprocessing import Pool
from scipy.interpolate import griddata

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger()

def regrid_xesmf(lis_input_file, gridinfo_file, regrided_tmp_file):
    
    print(f'Processing {lis_input_file}')
    # Read LIS grid information
    with nc.Dataset(lis_input_file, 'r') as f_lis_input:
        lat = f_lis_input.variables['lat'][:]
        lon = f_lis_input.variables['lon'][:]

    # Load your regular lat-lon grid
    ds_in = xr.open_dataset(gridinfo_file)

    # Define target Lambert grid
    ds_out = xr.Dataset(
        { "lat": (["y", "x"], lat), 
          "lon": (["y", "x"], lon), })

    # Regrid
    regridder    = xe.Regridder(ds_in, ds_out, method="bilinear")
    ds_regridded = regridder(ds_in)
    ds_regridded.to_netcdf(regrided_tmp_file)

    return

def regrid_vec_to_wrf_domain(n, case_name, lis_input_path, lis_out_path, gridinfo_file):

    start              = time.time()

    # define files
    lis_input_file     = f'{lis_input_path}/LIS_offline/lis_input.d0{n}.nc'
    lis_out_file       = f'{lis_out_path}/{case_name}/lis_input.d0{n}.nc'
    regrided_tmp_file  = f'{lis_out_path}/{case_name}/Openlandmap_soilcomposition_CORDEX_180E_depth_varying_tmp.d0{n}.nc'

    # Copy file
    if not os.path.exists(f'{lis_out_path}/{case_name}'):
        os.makedirs(f'{lis_out_path}/{case_name}')
    shutil.copy(lis_input_file, lis_out_file)

    # interpolate lis_input file
    regrid_xesmf(lis_input_file, gridinfo_file, regrided_tmp_file)
        
    # Get the regrided data 
    with nc.Dataset(regrided_tmp_file, 'r') as f_regrid:
        sand_vec     = f_regrid.variables['SAND_VEC'][:]
        clay_vec     = f_regrid.variables['CLAY_VEC'][:]
        silt_vec     = f_regrid.variables['SILT_VEC'][:]
        oc_vec       = f_regrid.variables['OC_VEC'][:]
        bulk_den_vec = f_regrid.variables['BULK_DEN_VEC'][:]
    
    # delete the median file
    os.remove(regrided_tmp_file)
    
    # Delete the pervious file 
    if not os.path.exists(lis_out_file):
        os.remove(lis_out_file)

    # Read LIS grid information
    with nc.Dataset(lis_input_file, 'r+') as f_lis_input:
        print('Mask out ocean')
        mask_lis     = f_lis_input.variables['LANDMASK'][:]
        
        for i in np.arange(6):
            
            sand_vec[i,:,:]     = np.where(mask_lis, sand_vec[i,:,:], -9999.)
            clay_vec[i,:,:]     = np.where(mask_lis, clay_vec[i,:,:], -9999.)
            silt_vec[i,:,:]     = np.where(mask_lis, silt_vec[i,:,:], -9999.)
            oc_vec[i,:,:]       = np.where(mask_lis, oc_vec[i,:,:], -9999.)
            bulk_den_vec[i,:,:] = np.where(mask_lis, bulk_den_vec[i,:,:], -9999.)

        # Create soil_depth dimension if not present
        if 'soil_depth' not in f_lis_input.dimensions:
            f_lis_input.createDimension('soil_depth', 6)

        print('Outputing data')
        # Add variable
        sand_vec_out       = f_lis_input.createVariable('SAND_VEC', 'f8', ('soil_depth', 'north_south', 'east_west'), fill_value=-9999)
        original_sand      = f_lis_input.variables['SAND']
        for attr_name in original_sand.ncattrs():
            sand_vec_out.setncattr(attr_name, original_sand.getncattr(attr_name))
        sand_vec_out[:]    = sand_vec

        clay_vec_out       = f_lis_input.createVariable('CLAY_VEC', 'f8', ('soil_depth', 'north_south', 'east_west'), fill_value=-9999)
        original_clay      = f_lis_input.variables['CLAY']
        for attr_name in original_clay.ncattrs():
            sand_vec_out.setncattr(attr_name, original_clay.getncattr(attr_name))
        clay_vec_out[:]    = clay_vec

        silt_vec_out       = f_lis_input.createVariable('SILT_VEC', 'f8', ('soil_depth', 'north_south', 'east_west'), fill_value=-9999)
        original_silt      = f_lis_input.variables['SILT']
        for attr_name in original_silt.ncattrs():
            sand_vec_out.setncattr(attr_name, original_silt.getncattr(attr_name))
        silt_vec_out[:]    = silt_vec

        ocsoil_vec_out     = f_lis_input.createVariable('OCSOIL_VEC', 'f8', ('soil_depth', 'north_south', 'east_west'), fill_value=-9999)
        original_ocsoil    = f_lis_input.variables['OCSOIL']
        for attr_name in original_ocsoil.ncattrs():
            sand_vec_out.setncattr(attr_name, original_ocsoil.getncattr(attr_name))
        ocsoil_vec_out[:]  = oc_vec

        blukdensity_out    = f_lis_input.createVariable('BULKDENSITY_VEC', 'f8', ('soil_depth', 'north_south', 'east_west'), fill_value=-9999)
        original_bulkden   = f_lis_input.variables['BULKDENSITY']
        for attr_name in original_bulkden.ncattrs():
            sand_vec_out.setncattr(attr_name, original_bulkden.getncattr(attr_name))
        blukdensity_out[:] = bulk_den_vec

    gc.collect()
    end = time.time()
    print(f"Execution time of finalize_lis_file: {end - start} seconds")

    return

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Interpolate _vec to lis_input.d0X.nc.')
    parser.add_argument('--case_name', type=str, default=None, help='What is the name of the case study')
    parser.add_argument('--n', type=int, default=None, help='The number of domain')

    args      = parser.parse_args()
    n         = args.n
    case_name = args.case_name

    lis_input_path  = f'/scratch/w97/mm3972/model/NUWRF/Drought_breaking_rainfall/{case_name}'
    lis_out_path    = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/for_LIS_WRF/lis_input'
    gridinfo_file   = '/g/data/w97/mm3972/scripts/wrf_scripts/make_LIS_landinfo/nc_file/Openlandmap_soilcomposition_CORDEX_180E_depth_varying.nc'
    
    print('n, case_name', n, case_name)
    regrid_vec_to_wrf_domain(n, case_name, lis_input_path, lis_out_path, gridinfo_file)
