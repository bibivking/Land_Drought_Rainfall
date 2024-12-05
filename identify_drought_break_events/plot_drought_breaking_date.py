import os
import sys
import glob
import numpy as np
import pandas as pd
import netCDF4 as nc
import cartopy.crs as ccrs
from datetime import datetime, timedelta
from scipy.ndimage import uniform_filter
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import matplotlib.ticker as mticker
from scipy.ndimage import label, sum as nd_sum
from cartopy.feature import NaturalEarthFeature, OCEAN
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from multiprocessing import Pool
from common_utils import *

def plot_drought_breaking_date(i, time, drought_break_day, precip_day, lat, lon):

    loc_lat  = [-45,-5]
    loc_lon  = [110,156]
    year     = time.year
    print(time)

    # print('drought_break_day', drought_break_day)

    if np.nansum(drought_break_day)>10:

        # ================== Start Plotting =================
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[5,5],sharex=False,
                                sharey=False, squeeze=True, subplot_kw={'projection': ccrs.PlateCarree()})

        plt.subplots_adjust(wspace=-0.15, hspace=0.105)

        plt.rcParams['text.usetex']     = False
        plt.rcParams['font.family']     = "sans-serif"
        plt.rcParams['font.serif']      = "Helvetica"
        plt.rcParams['axes.linewidth']  = 1.5
        plt.rcParams['axes.labelsize']  = 12
        plt.rcParams['font.size']       = 12
        plt.rcParams['legend.fontsize'] = 12
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12

        almost_black                    = '#262626'
        # change the tick colors also to the almost black
        plt.rcParams['ytick.color']     = almost_black
        plt.rcParams['xtick.color']     = almost_black

        # change the text colors also to the almost black
        plt.rcParams['text.color']      = almost_black

        # Change the default axis colors from black to a slightly lighter black,
        # and a little thinner (0.5 instead of 1)
        plt.rcParams['axes.edgecolor']  = almost_black
        plt.rcParams['axes.labelcolor'] = almost_black

        # set the box type of sequence number
        props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')

        states= NaturalEarthFeature(category="cultural", scale="50m",
                                            facecolor="none",
                                            name="admin_1_states_provinces_shp")

        # Download and add the states and coastlines
        states = NaturalEarthFeature(category="cultural", scale="50m",
                                            facecolor="none",
                                            name="admin_1_states_provinces_shp")

        ax.add_feature(states, linewidth=.5, edgecolor="black")
        ax.coastlines('50m', linewidth=0.8)

        # Set map extent for specific region
        ax.set_extent([110, 156, -45, -5])
        ax.coastlines(resolution="50m", linewidth=1)

        # Plot the original accumulated day data
        precip_day = np.where(drought_break_day < 1, np.nan, precip_day)
        im1        = ax.imshow(precip_day[:, :], vmin=0, vmax=30, cmap='YlGnBu', extent=[lon.min(), lon.max(), lat.min(), lat.max()])
        cbar1      = plt.colorbar(im1, ax=ax, orientation="horizontal", pad=0.02, aspect=20, shrink=0.7)
        cbar1.set_label('Rainfall (mm/day)')
        # ax.text(0.02, 0.95, texts[i], transform=ax.transAxes, fontsize=12, bbox=props)

        output_path = f"/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/plots/drought_breaking_date/{year}"

        print(f'output_path is {output_path}')

        if not os.path.exists(output_path):  # Check if the path does NOT exist
            os.mkdir(output_path)  # Create the directory
        plt.tight_layout()

        plt.savefig(f"{output_path}/Drought_breaking_{time}" , dpi=300, bbox_inches="tight")

    return

def plot_drought_breaking_multi_date():

    input_file = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_90_reorder_nan_filled_drought_breaking_date.nc'
    mask_file  = '/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_DLCM_mask.nc'
    rain_file  = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/agcd_v1_precip_total_r005_daily_1950_2023.nc'

    with nc.Dataset(input_file, mode='r') as f_in:
        lat_out            = f_in.variables['lat'][:]
        lon_out            = f_in.variables['lon'][:]
        drought_break_days = f_in.variables['drought_break_days'][:,:,:]
        time               = nc.num2date(f_in.variables['time'][:],f_in.variables['time'].units,
                             only_use_cftime_datetimes=False,only_use_python_datetimes=True)
        ntime              = len(time)

    with nc.Dataset(rain_file, mode='r') as f_rain:
        lat_rain         = f_rain.variables['lat'][10:]
        lon_rain         = f_rain.variables['lon'][:-45]
        precip           = f_rain.variables['precip'][:,10:,:-45]
        time_rain        = nc.num2date(f_rain.variables['time'][:],f_rain.variables['time'].units,
                            only_use_cftime_datetimes=False,only_use_python_datetimes=True)
        ntime_rain       = len(time)

    with nc.Dataset(mask_file, mode='r') as f_mask:
        landsea          = f_mask.variables['landsea'][:,:]

    landsea_3d           = np.repeat(landsea[np.newaxis, :, :], ntime, axis=0)
    drought_break_days   = np.where(landsea_3d==0, drought_break_days, np.nan)
    precip               = np.where(landsea_3d==0, precip, np.nan)
    print('finish reading')
    # Prepare arguments for parallel processing
    args_list = [ (i, time[i], drought_break_days[i,:,:], precip[i,:,:],
                    lat_out, lon_out) for i in np.arange(18263, ntime)]

    # Use multiprocessing to calculate accumulated days for each grid point
    with Pool() as pool:
        pool.starmap(plot_drought_breaking_date, args_list)

    return

# Main script
if __name__ == "__main__":

    plot_drought_breaking_multi_date()
