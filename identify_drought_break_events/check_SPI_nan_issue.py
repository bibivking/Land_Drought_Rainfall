import os
import gc
import sys
import glob
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from common_utils import *

# Paths
file_spi_30 = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_30_reorder.nc'
file_spi_90 = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_90_reorder.nc'
file_mask   = '/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_CSIRO_AU_NAT_ELEV_DLCM_mask.nc'
file_rain   = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/agcd_v1_precip_total_r005_daily_1970_2023.nc'

f_spi       = nc.Dataset(file_spi_30, 'r')
spi         = f_spi.variables['spi_pearson_30'][:,10:,:-45]

f_rain      = nc.Dataset(file_rain, 'r')
rain        = f_rain.variables['precip'][:,10:,:-45]
lat_rian    = f_rain.variables['lat'][10:]
lon_rain    = f_rain.variables['lon'][:-45]
lons, lats  = np.meshgrid(lon_rain, lat_rian)

f_mask      = nc.Dataset(file_mask, 'r')
landsea     = f_mask.variables['landsea']
lat_mask    = f_rain.variables['lat']
lon_mask    = f_rain.variables['lon']

ntime       = len(spi[:,0,0])

for i in np.arange(29,ntime):

    spi_tmp  = np.where(landsea == 1, -9999., spi[i,:,:].data)
    rain_tmp = np.where(landsea == 1, np.nan, rain[i,:,:].data  )
    rain_tmp = np.where(np.isnan(spi_tmp), rain_tmp, np.nan     )

    if np.sum(np.isnan(spi_tmp)) > 0:
        print(i,' ', np.sum(np.isnan(spi_tmp)))

        # Create a figure with 3 subplots in 1 row
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})

        # Adjust spacing between subplots
        plt.subplots_adjust(wspace=0.3)

        # Loop through each axis
        for ax in axs:

            # Set extent based on loc_lat and loc_lon
            ax.set_extent([112, 156, -45, -10])  # Example extent, adjust as needed

            ax.coastlines(resolution="50m", linewidth=1)

            # Add gridlines
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='black', linestyle='--')
            gl.xlabels_top = False
            gl.ylabels_right = False
            gl.xlines = True

            gl.xlocator = mticker.FixedLocator([140, 145, 150])
            gl.ylocator = mticker.FixedLocator([-40, -35, -30])

            gl.xformatter = LongitudeFormatter()
            gl.yformatter = LatitudeFormatter()
            gl.xlabel_style = {'size': 10, 'color': 'black'}
            gl.ylabel_style = {'size': 10, 'color': 'black'}

        plot1 = axs[0].contourf(lons, lats, np.where(spi_tmp == -9999., np.nan, spi_tmp), transform=ccrs.PlateCarree(), extend='both', cmap=plt.cm.BrBG) # levels=clevs, 
        plot2 = axs[1].contourf(lons, lats, rain_tmp, transform=ccrs.PlateCarree(), extend='both', cmap=plt.cm.BrBG) # levels=clevs, 

        cb = plt.colorbar(plot1, ax=axs[0], orientation="horizontal", pad=0.02, aspect=16, shrink=0.8)
        cb = plt.colorbar(plot2, ax=axs[1], orientation="horizontal", pad=0.02, aspect=16, shrink=0.8)
        cb.ax.tick_params(labelsize=10)

        plt.savefig(f'./plots/check_SPI_issue/SPI_nan_day{i}.png',dpi=300)
