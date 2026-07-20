import gc
import os
import time
import argparse
import logging
import numpy as np
import xesmf as xe
import xarray as xr
import netCDF4 as nc
from datetime import datetime, timedelta

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger()

def regrid_xesmf(lis_inout_file, gridinfo_file, regrided_tmp_file):

    """
    Regrid the gridinfo file to match the LIS input domain using xESMF.
    """
    
    print(f'Processing regrid_xesmf')
    
    # Read LIS grid information
    with nc.Dataset(lis_inout_file, 'r') as f_lis_input:
        lat = f_lis_input.variables['lat'][:]
        lon = f_lis_input.variables['lon'][:]

    # Load your regular lat-lon grid
    ds_in = xr.open_dataset(gridinfo_file)

    # Define target Lambert grid
    ds_out = xr.Dataset(
        {"lat": (["y", "x"], lat),
         "lon": (["y", "x"], lon), })

    # Regrid
    regridder    = xe.Regridder(ds_in, ds_out, method="bilinear")
    ds_regridded = regridder(ds_in)
    ds_regridded.to_netcdf(regrided_tmp_file)

    return

def regrid_albedo_lai_to_wrf_domain(n, case_name, is_clim, start_date, end_date, lis_inout_path, gridinfo_path):

    """
    Main function to regrid albedo and LAI to the WRF domain.
    """
    
    start = time.time()

    # Define input LIS file path
    lis_inout_file = f'{lis_inout_path}/{case_name}/lis_input.d0{n}.nc'

    date_s  = datetime.strptime(start_date, '%Y%m%d')
    date_e  = datetime.strptime(end_date,   '%Y%m%d')
    one_day = timedelta(days=1)

    # Generate all required dates
    dates   = [date_s + i * one_day for i in range((date_e - date_s).days + 1)]

    # Set previous year
    previous_year = date_s.year

    # Open the LIS input file
    f_lis_input   = nc.Dataset(lis_inout_file, 'r+')

    # Read in LANDMASK
    mask_lis      = f_lis_input.variables['LANDMASK'][:]

    for i, date in enumerate(dates):
        current_year = date.year
        doy          = int(date.strftime('%j'))

        print('Processing ',current_year, doy )
        
        # Regrid the gridinfo file if year changes
        if i == 0 or previous_year != current_year:
            
            # Set input gridinfo file name 
            if is_clim == 'True':
                # Handle leap years
                if current_year % 4 == 0:
                    gridinfo_file = f'{gridinfo_path}/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_clim_leap.nc'
                else:
                    gridinfo_file = f'{gridinfo_path}/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_clim_common.nc'
            else:
                gridinfo_file = f'{gridinfo_path}/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_time_varying_{current_year}.nc'

            # Set regrided gridinfo file name
            regrided_tmp_file = f'{lis_inout_path}/gridinfo_tmp.d0{n}.nc'
            
            # Delete pervious regrided file if it exists
            if os.path.exists(regrided_tmp_file):
                os.remove(regrided_tmp_file)

            # Regridding
            regrid_xesmf(lis_inout_file, gridinfo_file, regrided_tmp_file)

        # Read the regridded data
        with nc.Dataset(regrided_tmp_file, 'r') as f_regrid:
            albedo_modis = f_regrid.variables['Albedo_MODIS'][:, doy - 1, :, :]
            lai          = f_regrid.variables['LAI'][doy - 1, :, :]

        # Mask out ocean areas
        for rad in np.arange(4):
            albedo_modis[rad, :, :] = np.where(mask_lis == 1, albedo_modis[rad, :, :], -9999.)
        lai = np.where(mask_lis == 1, lai, -9999.)

        # Create and save variables in LIS input file
        print('Outputting data')
        
        lai_out            = f_lis_input.createVariable(f'LAI_{date.strftime("%Y%m%d")}', 'f8', ('north_south', 'east_west'), fill_value=-9999)
        lai_out[:]         = lai

        alb_BSA_vis_out    = f_lis_input.createVariable(f'ALBEDO_BSA_vis_{date.strftime("%Y%m%d")}', 'f8', ('north_south', 'east_west'), fill_value=-9999)
        alb_BSA_vis_out[:] = albedo_modis[0, :, :]
        
        alb_BSA_nir_out    = f_lis_input.createVariable(f'ALBEDO_BSA_nir_{date.strftime("%Y%m%d")}', 'f8', ('north_south', 'east_west'), fill_value=-9999)
        alb_BSA_nir_out[:] = albedo_modis[1, :, :]
        
        alb_WSA_vis_out    = f_lis_input.createVariable(f'ALBEDO_WSA_vis_{date.strftime("%Y%m%d")}', 'f8', ('north_south', 'east_west'), fill_value=-9999)
        alb_WSA_vis_out[:] = albedo_modis[2, :, :]
        
        alb_WSA_nir_out    = f_lis_input.createVariable(f'ALBEDO_WSA_nir_{date.strftime("%Y%m%d")}', 'f8', ('north_south', 'east_west'), fill_value=-9999)
        alb_WSA_nir_out[:] = albedo_modis[3, :, :]

        gc.collect()

        previous_year = current_year

    f_lis_input.close()
    end = time.time()
    print(f"Execution time of regrid_albedo_lai_to_wrf_domain: {end - start} seconds")

    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Interpolate MODIS albedo and LAI to lis_input.d0X.nc.')
    parser.add_argument('--n', type=int, help='The number of the domain')
    parser.add_argument('--is_clim', type=str, choices=['True', 'False'], help='Whether using climatology albedo and LAI')
    parser.add_argument('--case_name', type=str, help='Case study name')
    parser.add_argument('--start_date', type=str, help="Start date, format: 'YYYYMMDD'")
    parser.add_argument('--end_date', type=str, help="End date, format: 'YYYYMMDD'")

    args       = parser.parse_args()
    
    print(args)

    n          = args.n
    is_clim    = args.is_clim
    case_name  = args.case_name
    start_date = args.start_date
    end_date   = args.end_date

    lis_inout_path = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/for_LIS_WRF/lis_input'
    gridinfo_path  = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/SG_smooth_anomaly'

    print('n, is_clim, case_name, start_date & end_date:', n, is_clim, case_name, start_date, end_date)
    regrid_albedo_lai_to_wrf_domain(n, case_name, is_clim, start_date, end_date, lis_inout_path, gridinfo_path)
