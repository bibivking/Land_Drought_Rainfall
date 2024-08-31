import os
import rasterio
import argparse
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta

def generate_dominant_type_map(resolution):
    if resolution == "5km":
        C3_path     = "/g/data/w97/mm3972/data/AU_C3_C4_plant_map/nc_file/Proportional_C3_vegation_5km.nc"
        C4_path     = "/g/data/w97/mm3972/data/AU_C3_C4_plant_map/nc_file/Proportional_C4_vegation_5km.nc"
        output_path = "/g/data/w97/mm3972/data/AU_C3_C4_plant_map/nc_file/C3_or_C4_dominate_map.nc"
    elif resolution == "MODIS_res":
        C3_path     = "/g/data/w97/mm3972/data/AU_C3_C4_plant_map/nc_file/Proportional_C3_vegation_MODIS_res.nc"
        C4_path     = "/g/data/w97/mm3972/data/AU_C3_C4_plant_map/nc_file/Proportional_C4_vegation_MODIS_res.nc"
        output_path = "/g/data/w97/mm3972/data/AU_C3_C4_plant_map/nc_file/C3_or_C4_dominate_map_MODIS_res.nc"
    elif resolution == "NVIS_res":
        C3_path     = "/g/data/w97/mm3972/data/AU_C3_C4_plant_map/nc_file/Proportional_C3_vegation_NVIS_res.nc"
        C4_path     = "/g/data/w97/mm3972/data/AU_C3_C4_plant_map/nc_file/Proportional_C4_vegation_NVIS_res.nc"
        output_path = "/g/data/w97/mm3972/data/AU_C3_C4_plant_map/nc_file/C3_or_C4_dominate_map_NVIS_res.nc"

    # Read C3 and C4 fraction
    f_C3        = nc.Dataset(C3_path, 'r', format='NETCDF4')
    C3_fraction = f_C3.variables['C3_proportion'][:,:]
    lat_in      = f_C3.variables['lat'][:]
    lon_in      = f_C3.variables['lon'][:]
    f_C3.close()

    f_C4        = nc.Dataset(C4_path, 'r', format='NETCDF4')
    C4_fraction = f_C4.variables['C4_proportion'][:,:]
    f_C4.close()

    var_output  = np.where(C3_fraction >= C4_fraction, 3, 4)

    # Output C3 C4 dominant

    f_out        = nc.Dataset(output_path, 'w', format='NETCDF4')

    f_out.creation_date = "%s" % (datetime.now())
    f_out.Conventions   = "CF-1.0"
    f_out.CDI           = "Climate Data Interface version 2.1.0 (https://mpimet.mpg.de/cdi)"
    f_out.source        = C3_path+ " and "+C4_path
    f_out.crs           = "EPSG:4326"
    f_out.history       = "Created by: %s" % (os.path.basename(__file__))

    # set dimensions
    f_out.createDimension('lat', len(lat_in))
    f_out.createDimension('lon', len(lon_in))

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

    var_out                = f_out.createVariable('Ctype', 'f4', ('lat', 'lon'), fill_value=-9999.)
    var_out.long_name      = 'dominant C type, 3: C3, 4: C4'
    var_out.missing_value  = -9999.
    var_out[:]             = var_output
    f_out.close()

    lat_out = None
    lon_out = None
    var_out = None

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Regrid to coarse resolution.')
    parser.add_argument('--resolution', type=str, default=None, help='MODIS_res, NVIS_res or 5km.')

    args = parser.parse_args()

    generate_dominant_type_map(args.resolution)
