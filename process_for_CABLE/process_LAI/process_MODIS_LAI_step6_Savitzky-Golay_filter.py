import os
import glob
import xarray as xr
import numpy as np
import argparse
from scipy.signal import savgol_filter
import pandas as pd
import netCDF4 as nc

def run_savizky_golay_filter(file_in, file_out, window_length=15, polyorder=2):

    # Load the NetCDF file into an xarray Dataset
    ds = xr.open_dataset(file_in, engine="netcdf4")

    # Chunk the dataset along 'lat' and 'lon' dimensions
    # Example chunk sizes: 50 for lat and lon, time remains unchunked
    ds_chunked = ds.chunk({"time": -1, "latitude": 50, "longitude": 50})

    # Extract the LAI array
    lai        = ds_chunked["Lai_500m"]
    print(lai.chunks)  # Verify chunking

    # Apply the filter along the time dimension
    smoothed_lai = xr.apply_ufunc(
        savgol_filter,
        lai,                # First argument: has a "time" dimension
        window_length,      # Second argument: scalar, no dimensions
        polyorder,          # Third argument: scalar, no dimensions
        kwargs={"axis": 0},  # Smooth along the time dimension
        # input_core_dims=[["time"], [], []],  # One entry per argument, for my case savgol_filter(lai, window_length,polyorder), ["time"] is for lai, the two [] are for window_length and polyorder
        # output_core_dims=[["time"]],         # Output has a "time" dimension
        input_core_dims=[["time"], [], []],  # One entry per argument, for my case savgol_filter(lai, window_length,polyorder), ["time"] is for lai, the two [] are for window_length and polyorder
        output_core_dims=[["time"]],         # Output has a "time" dimension
        vectorize=True,                      # Enable vectorization
        dask="parallelized",                 # Use Dask if data is chunked
        output_dtypes=[lai.dtype],           # Keep the same data type
    )

    # Copy the attributes from Lai_500m to the smoothed version
    smoothed_lai.attrs = lai.attrs.copy()

    # Remove the 'Lai_500m' variable from the dataset
    ds_chunked         = ds_chunked.drop("Lai_500m")

    # Rename the smoothed variable to "Lai_500m"
    smoothed_lai.name  = "Lai_500m"

    # Reorder the dimensions to (time, lat, lon)
    smoothed_lai       = smoothed_lai.transpose("time", "latitude", "longitude")

    # Add the smoothed data back to the Dataset
    ds_chunked["Lai_500m"] = smoothed_lai

    # Check if the file exists
    if os.path.exists(file_out):
         os.remove(file_out)

    # Save to a new NetCDF file
    ds_chunked.to_netcdf(file_out)

    print(f"Smoothed LAI saved to {file_out}")

    return

if __name__ == "__main__":

    # Define Savitzky-Golay filter parameters
    parser = argparse.ArgumentParser(description='set parameters for Savitzky-Golay filter')
    parser.add_argument('--window', type=int,default=15, help='window_length')
    parser.add_argument('--order', type=int, default=2, help='polyorder')

    args = parser.parse_args()

    file_in  = "/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km/MCD15A3H.061_500m_aid0001_LAI_regridded_2002-2024.nc"
    file_out = f"/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km/MCD15A3H.061_500m_aid0001_LAI_regridded_2002-2024_SG_filter_window={args.window}timestep_order={args.order}.nc"

    run_savizky_golay_filter(file_in, file_out, args.window, args.order)
