from netCDF4 import Dataset
import numpy as np

def create_and_copy_agcd_rainfall(input_file, output_file):

    # Open the original dataset
    src_dataset = Dataset(input_file, 'r')

    # Create a new NetCDF file
    dst_dataset = Dataset(output_file, 'w', format='NETCDF4')

    # Copy dimensions
    for name, dimension in src_dataset.dimensions.items():
        dst_dataset.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

    # Copy variables
    for name, variable in src_dataset.variables.items():
        dst_var = dst_dataset.createVariable(name, variable.datatype, variable.dimensions)
        dst_var.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
        dst_var[:] = variable[:]

    # Copy global attributes
    dst_dataset.setncatts({k: src_dataset.getncattr(k) for k in src_dataset.ncattrs()})

    # Close datasets
    src_dataset.close()
    dst_dataset.close()

    # Verify the new file
    new_dataset = Dataset(output_file, 'r')
    print("Dimensions:")
    for name, dimension in new_dataset.dimensions.items():
        print(f"{name}: {len(dimension)}")

    print("\nVariables:")
    for name, variable in new_dataset.variables.items():
        print(f"{name}: {variable.shape}")

    print("\nGlobal Attributes:")
    for attr in new_dataset.ncattrs():
        print(f"{attr}: {new_dataset.getncattr(attr)}")
    new_dataset.close()
    return

if __name__ == "__main__":

    # Paths
    for year in np.arange(1971,2024,1):
        input_file = '/g/data/zv2/agcd/v1-0-2/precip/total/r005/01day/agcd_v1_precip_total_r005_daily_'+str(year)+'.nc'
        output_file = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/calc_AUS_SPI/nc_files/AGCD_rain_daily_1970_2023/agcd_v1_precip_total_r005_daily_'+str(year)+'.nc'
        create_and_copy_agcd_rainfall(input_file, output_file)