import rasterio
import argparse
import netCDF4 as nc
import numpy as np


def tiff_to_nc(var_name, tiff_path, netcdf_path):

    # Open the TIFF file
    with rasterio.open(tiff_path) as src:
        # Read the data from the TIFF file
        tiff_data = src.read(1)

        if 'cover' in var_name:
            tiff_data = np.where(tiff_data<0., -9999., tiff_data)
        elif 'proportion' in var_name:
            tiff_data = np.where(tiff_data>1., -9999., tiff_data)
            tiff_data = np.where(tiff_data<0., -9999., tiff_data)

        # Extract the geospatial metadata
        transform = src.transform
        width     = src.width
        height    = src.height
        crs       = src.crs
        
        print('transform, width, height, crs',transform, width, height, crs)
        # Define latitude and longitude arrays
        lon = np.linspace(108.6900885, 108.6900885 + 0.000966697366395897 * width, width)
        lat = np.linspace(-6.4994993, -6.4994993 - 0.000966697366395897 * height, height)
        
        # Reverse latitude to ensure north to south orientation
        lat = lat[::-1]

    # Create a NetCDF file
    with nc.Dataset(netcdf_path, 'w', format='NETCDF4') as dst:
        # Create dimensions
        dst.createDimension('lat', len(lat))
        dst.createDimension('lon', len(lon))
        
        # Create coordinate variables
        latitudes  = dst.createVariable('lat', 'f4', ('lat',))
        longitudes = dst.createVariable('lon', 'f4', ('lon',))
        
        # Assign data to coordinate variables
        latitudes[:]  = lat
        longitudes[:] = lon
        
        # Assign metadata to coordinate variables
        latitudes.units    = 'degrees_north'
        longitudes.units   = 'degrees_east'
        
        # Create the main variable to store the data
        data_var = dst.createVariable(var_name, 'f4', ('lat', 'lon'), zlib=True, fill_value=-9999) 
        # the zlib=True argument is used to enable compression for the variable
        if 'cover' in var_name:
            data_var.units     = '%'  # Change as necessary
            data_var.long_name = f'{var_name}age, can >100%, do not know why'
        elif 'proportion' in var_name:
            data_var.units     = 'fraction'  # Change as necessary
            data_var.long_name = f'{var_name}, 0-1'
        
        # Assign the data to the NetCDF variable
        data_var[:, :]     = tiff_data[::-1, :]  # Flip the data along the latitude axis

        # Add global attributes
        dst.description    = f'NetCDF file created from tiff file' 
        dst.source         = tiff_path
        dst.crs            = str(crs)
        
    print(f"Conversion complete: {netcdf_path}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Translate tiff file to nc file.')
    parser.add_argument('--var_name',    type=str, default=None, help='The name of output variable ')
    parser.add_argument('--tiff_path',   type=str, default=None, help='Input tiff file path.')
    parser.add_argument('--netcdf_path', type=str, default=None, help='Output nc file path.')

    args = parser.parse_args()

    tiff_to_nc(args.var_name, args.tiff_path, args.netcdf_path)
