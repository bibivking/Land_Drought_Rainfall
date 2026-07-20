import os
import argparse
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta

def regrid_C_type_to_5km(var_name, netcdf_path, output_path):

    # File paths
    input_file       = netcdf_path
    output_file      = output_path
    target_grid_file = '/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc'

    # Open input file and read landcover, lat and lon information
    f_in             = nc.Dataset(input_file,'r')
    Ctype_fine       = f_in.variables[var_name][:,:]
    Ctype_fine       = np.where(Ctype_fine>1., np.nan, Ctype_fine)
    Ctype_fine       = np.where(Ctype_fine<0., np.nan, Ctype_fine)
    lat_fine         = f_in.variables['lat'][:]
    lon_fine         = f_in.variables['lon'][:]

    # Read target lat and lon information
    f_grid           = nc.Dataset(target_grid_file,'r')
    lat_grid         = f_grid.variables['latitude'][:]
    lon_grid         = f_grid.variables['longitude'][:]
    nlat_grid        = len(lat_grid)
    nlon_grid        = len(lon_grid)
    Ctype_coarse     = np.zeros((nlat_grid, nlon_grid))

    # Loop input lat and lon
    for i, lat in enumerate(lat_grid):
        for j, lon in enumerate(lon_grid):

            # fine-res pixels in the coarse grid cell
            lat_indices = np.where((lat_fine > lat-0.025) & (lat_fine < lat+0.025))[0]
            lon_indices = np.where((lon_fine > lon-0.025) & (lon_fine < lon+0.025))[0]

            lat_index_s = lat_indices[0]
            lat_index_e = lat_indices[-1]

            lon_index_s = lon_indices[0]
            lon_index_e = lon_indices[-1]

            Ctype_coarse[i, j] = np.nanmean(Ctype_fine[lat_index_s:lat_index_e+1,lon_index_s:lon_index_e+1])

    Ctype_coarse = np.where(np.isnan(Ctype_coarse), -9999., Ctype_coarse)

    # =========== Create nc file ===========
    f_out               = nc.Dataset(output_file, 'w', format='NETCDF4')

    f_out.creation_date = "%s" % (datetime.now())
    f_out.Conventions   = "CF-1.0"
    f_out.CDI           = "Climate Data Interface version 2.1.0 (https://mpimet.mpg.de/cdi)"
    f_out.source        = "/g/data/w97/mm3972/data/AU_C3_C4_plant_map/Total_C3_cover.tif"
    f_out.description   = "NetCDF file created from tiff file"
    f_out.crs           = "EPSG:4326"
    f_out.history       = "Created by: %s" % (os.path.basename(__file__))

    # set dimensions
    f_out.createDimension('lat', nlat_grid)
    f_out.createDimension('lon', nlon_grid)

    lat_out                = f_out.createVariable('lat', 'f4', ('lat'))
    lat_out.long_name      = "Latitude"
    lat_out.standard_name  = "latitude"
    lat_out.axis           = "Y"
    lat_out.units          = "degrees_north"
    lat_out[:]             = lat_grid

    lon_out                = f_out.createVariable('lon', 'f4', ('lon'))
    lon_out.long_name      = "Longitude"
    lon_out.standard_name  = "longitude"
    lon_out.axis           = "X"
    lon_out.units          = "degrees_east"
    lon_out[:]             = lon_grid

    var_out                = f_out.createVariable(var_name, 'f4', ('lat', 'lon'), fill_value=-9999.)
    var_out.long_name      = var_name
    var_out.units          = "fraction"
    var_out.missing_value  = -9999.
    var_out[:]             = Ctype_coarse
    f_out.close()

    lat_out = None
    lon_out = None
    var_out = None

    return

def regrid_C_type_to_MODIS_landcover_res(var_name, netcdf_path, output_path):

    # File paths
    input_file       = netcdf_path
    output_file      = output_path
    target_grid_file = '/g/data/w97/mm3972/data/MODIS/MODIS_landcover/var_landcover_only/MCD12Q1.061_500m_aid0001_2001-2022_PFT.nc'

    # Open input file and read landcover, lat and lon information
    f_in             = nc.Dataset(input_file,'r')
    Ctype_fine       = f_in.variables[var_name][:,:]
    Ctype_fine       = np.where(Ctype_fine>1., np.nan, Ctype_fine)
    Ctype_fine       = np.where(Ctype_fine<0., np.nan, Ctype_fine)
    lat_fine         = f_in.variables['lat'][:]
    lon_fine         = f_in.variables['lon'][:]

    # Read target lat and lon information
    f_grid           = nc.Dataset(target_grid_file,'r')
    lat_grid         = f_grid.variables['lat'][:]
    lon_grid         = f_grid.variables['lon'][:]
    nlat_grid        = len(lat_grid)
    nlon_grid        = len(lon_grid)
    Ctype_coarse     = np.zeros((nlat_grid, nlon_grid))
    interval         = 0.00416667

    # Loop input lat and lon
    for i, lat in enumerate(lat_grid):
        for j, lon in enumerate(lon_grid):

            # fine-res pixels in the coarse grid cell
            lat_indices = np.where((lat_fine > lat-interval/2.) & (lat_fine < lat+interval/2.))[0]
            lon_indices = np.where((lon_fine > lon-interval/2.) & (lon_fine < lon+interval/2.))[0]
            try:
                lat_index_s = lat_indices[0]
                lat_index_e = lat_indices[-1]

                lon_index_s = lon_indices[0]
                lon_index_e = lon_indices[-1]

                Ctype_coarse[i, j] = np.nanmean(Ctype_fine[lat_index_s:lat_index_e+1,lon_index_s:lon_index_e+1])
            except:
                # if no fine resolution pixel in the coarse resolution grid cell, set proportion as 0,
                # the C3 proportion = C4 proportion, define as C3 type
                Ctype_coarse[i, j] = 0.

    Ctype_coarse = np.where(np.isnan(Ctype_coarse), -9999., Ctype_coarse)

    # =========== Create nc file ===========
    f_out               = nc.Dataset(output_file, 'w', format='NETCDF4')

    f_out.creation_date = "%s" % (datetime.now())
    f_out.Conventions   = "CF-1.0"
    f_out.CDI           = "Climate Data Interface version 2.1.0 (https://mpimet.mpg.de/cdi)"
    f_out.source        = "/g/data/w97/mm3972/data/AU_C3_C4_plant_map/Total_C3_cover.tif"
    f_out.description   = "NetCDF file created from tiff file"
    f_out.crs           = "EPSG:4326"
    f_out.history       = "Created by: %s" % (os.path.basename(__file__))

    # set dimensions
    f_out.createDimension('lat', nlat_grid)
    f_out.createDimension('lon', nlon_grid)

    lat_out                = f_out.createVariable('lat', 'f4', ('lat'))
    lat_out.long_name      = "Latitude"
    lat_out.standard_name  = "latitude"
    lat_out.axis           = "Y"
    lat_out.units          = "degrees_north"
    lat_out[:]             = lat_grid

    lon_out                = f_out.createVariable('lon', 'f4', ('lon'))
    lon_out.long_name      = "Longitude"
    lon_out.standard_name  = "longitude"
    lon_out.axis           = "X"
    lon_out.units          = "degrees_east"
    lon_out[:]             = lon_grid

    var_out                = f_out.createVariable(var_name, 'f4', ('lat', 'lon'), fill_value=-9999.)
    var_out.long_name      = var_name
    var_out.units          = "fraction"
    var_out.missing_value  = -9999.
    var_out[:]             = Ctype_coarse
    f_out.close()

    lat_out = None
    lon_out = None
    var_out = None

    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Regrid to coarse resolution.')
    parser.add_argument('--var_name',    type=str, default=None, help='The name of output variable ')
    parser.add_argument('--netcdf_path', type=str, default=None, help='Input fine res file path.')
    parser.add_argument('--output_path', type=str, default=None, help='Output file path.')
    parser.add_argument('--resolution', type=str, default=None, help='MODIS_res or 5km.')

    args = parser.parse_args()
    if args.resolution == '5km':
        regrid_C_type_to_5km(args.var_name, args.netcdf_path, args.output_path)
    elif args.resolution == 'MODIS_res':
        # 10 hour cannot finish running regrid_C_type_to_MODIS_landcover_res, so it is bad function 
        regrid_C_type_to_MODIS_landcover_res(args.var_name, args.netcdf_path, args.output_path)
