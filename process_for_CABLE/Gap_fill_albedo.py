import os
import gc
import argparse
import copy
import logging
import netCDF4 as nc
import numpy as np
from multiprocessing import Pool, cpu_count
from common_utils import *

def fill_gap_for_point(args):

    var, lat, lon, window = args

    time_series = var[:, lat, lon]

    nan_indices   = np.where(np.isnan(time_series))[0]
    valid_indices = np.where(~np.isnan(time_series))[0]

    if np.all(np.isnan(time_series)):
        return time_series

    if len(nan_indices) == 0:
        return time_series

    if len(valid_indices) <= 2:
        return time_series

    # logging.debug(f"Original time series at (lat={lat}, lon={lon}): {time_series}")
    # logging.debug(f"Valid indices at (lat={lat}, lon={lon}): {valid_indices}")
    # logging.debug(f"NaN indices at (lat={lat}, lon={lon}): {nan_indices}")

    for nan_index in nan_indices:
        distances = np.abs(valid_indices - nan_index)
        sorted_indices = np.argsort(distances)
        # logging.debug(f"Distances at (lat={lat}, lon={lon}): {distances}")
        # logging.debug(f"Sorted indices at (lat={lat}, lon={lon}): {sorted_indices}")

        closest_first_date = valid_indices[sorted_indices[0]]
        closest_second_date = valid_indices[sorted_indices[1]]

        if distances[sorted_indices[1]] == distances[sorted_indices[0]]:
            if_same_side = (closest_first_date - nan_index) * (closest_second_date - nan_index) > 0
            if if_same_side:
                closest_second_date = valid_indices[sorted_indices[2]]

        # logging.debug(f"First date at (lat={lat}, lon={lon}): {closest_first_date}, Second date at (lat={lat}, lon={lon}): {closest_second_date}")

        if (window > 0) and (distances[sorted_indices[1]] <= window):
            # and abs(closest_second_date - closest_first_date) <= window
            delta = (time_series[closest_second_date] - time_series[closest_first_date]) / \
                    (closest_second_date - closest_first_date)
            time_series[nan_index] = time_series[closest_first_date] + delta * (nan_index - closest_first_date)
        else:
            delta = (time_series[closest_second_date] - time_series[closest_first_date]) / \
                    (closest_second_date - closest_first_date)
            time_series[nan_index] = time_series[closest_first_date] + delta * (nan_index - closest_first_date)

        # logging.debug(f"Filled time series at (lat={lat}, lon={lon}): {time_series}\n")

    # var[:, lat, lon] = time_series ??? Why adding this sentence will fail pass the filled data to nc file
    logging.debug(f"Filled time_series at (lat={lat}, lon={lon}): {time_series}\n")

    return time_series

def gap_fill(var, window=30):

    nlat, nlon = var.shape[1], var.shape[2]
    var_fill   = copy.deepcopy(var)

    # Create the argument list for each process
    args = [(var, lat, lon, window) for lat in range(nlat) for lon in range(nlon)]

    # Use multiprocessing to parallelize the task
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(fill_gap_for_point, args)

    gc.collect()

    # Reassemble the results back into the original array shape
    for i, (lat, lon) in enumerate([(lat, lon) for lat in range(nlat) for lon in range(nlon)]):
        var_fill[:, lat, lon] = results[i]
        if ~np.all(np.isnan(results[i])):
            logging.debug(f"Original var[:, lat, lon] at (lat={lat}, lon={lon}): {var[:, lat, lon]}\n")
            logging.debug(f"Results at (lat={lat}, lon={lon}): {results[i]}\n")
        # logging.debug(f"Filled var_fill[:, lat, lon] at (lat={lat}, lon={lon}): {var_fill[:, lat, lon]}\n")

    return var_fill #var

if __name__ == "__main__":

    # Configure logging
    logging.basicConfig(filename='/g/data/w97/mm3972/scripts/Drought/Post_drought_rainfall/process_for_CABLE/process_albedo/gap_fill.log', level=logging.DEBUG)

    parser = argparse.ArgumentParser(description='Gap-fill MODIS albedo data.')
    parser.add_argument('file_path', type=str, help='Path to the netCDF file')
    parser.add_argument('--window', type=int, default=0, help='Window size for gap filling (default: 90 days)')
    parser.add_argument('--albedo_band', type=str, default="BSA_nir", help='Albedo band to process, BSA_nir, WSA_nir, BSA_vis, WSA_vis')

    args = parser.parse_args()

    scale_factor   = 0.001
    missing_value  = 32767
    var_name       = "Albedo_"+args.albedo_band

    f_gap_fill     = nc.Dataset(args.file_path, 'r+')
    var            = f_gap_fill.variables[var_name][:,:,:] # do not need to time scale_factor, since f_gap_fill.variables by default times scale_factor

    gc.collect()

    var            = np.where(var > 1, np.nan, var)
    var_gap_filled = gap_fill(var, window=args.window)
    
    var_gap_filled = np.where(np.isnan(var_gap_filled), missing_value*scale_factor, var_gap_filled)
    var_gap_filled = np.where(var_gap_filled>1, missing_value*scale_factor, var_gap_filled)
    var_gap_filled = np.where(var_gap_filled<0, missing_value*scale_factor, var_gap_filled)

    f_gap_fill.variables[var_name][:,:,:] = var_gap_filled
    f_gap_fill.close()
