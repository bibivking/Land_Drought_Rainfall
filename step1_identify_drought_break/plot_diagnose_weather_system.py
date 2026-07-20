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

def plot_diagnose_weather_system(time, lat_AGCD, lon_AGCD, accumulated_day, accumulated_day_large_drought, precip_day):
    
    def is_leap(year):
        if year%4 == 0:
            return True
        else:
            return False
    
    # time     = datetime(2017,1,2,0,0,0,0)
    year     = time.year
    month    = time.month
    day      = time.day

    dom_leap = {1:31, 2:29, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 9:30, 10:31, 11:30, 12:31}
    dom_com  = {1:31, 2:28, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 9:30, 10:31, 11:30, 12:31}

    if is_leap(year):
        dom = dom_leap
    else:
        dom = dom_com

    # ERA5 MSLP
    ERA5_MSLP_path    = "/g/data/rt52/era5/single-levels/reanalysis"
    ERA5_MSLP_T_file  = f'{ERA5_MSLP_path}/2t/{year}/2t_era5_oper_sfc_{year}{month:02d}01-{year}{month:02d}{dom[month]:02d}.nc'   # air temperature
    ERA5_MSLP_P_file  = f'{ERA5_MSLP_path}/sp/{year}/sp_era5_oper_sfc_{year}{month:02d}01-{year}{month:02d}{dom[month]:02d}.nc'   # surface pressure
    ERA5_MSLP_U_file  = f'{ERA5_MSLP_path}/10u/{year}/10u_era5_oper_sfc_{year}{month:02d}01-{year}{month:02d}{dom[month]:02d}.nc' # 10 m wind speed
    ERA5_MSLP_V_file  = f'{ERA5_MSLP_path}/10v/{year}/10v_era5_oper_sfc_{year}{month:02d}01-{year}{month:02d}{dom[month]:02d}.nc' # 10 m wind speed
    ERA5_MSLP_R_file  = f'{ERA5_MSLP_path}/tp/{year}/tp_era5_oper_sfc_{year}{month:02d}01-{year}{month:02d}{dom[month]:02d}.nc'   # Total rainfall

    # ERA5 Pressure level
    ERA5_Press_path    = "/g/data/rt52/era5/pressure-levels/reanalysis"
    ERA5_Press_T_file  = f'{ERA5_Press_path}/t/{year}/t_era5_oper_pl_{year}{month:02d}01-{year}{month:02d}{dom[month]:02d}.nc'   # air temperature
    ERA5_Press_Z_file  = f'{ERA5_Press_path}/z/{year}/z_era5_oper_pl_{year}{month:02d}01-{year}{month:02d}{dom[month]:02d}.nc'   # Geopotential
    ERA5_Press_U_file  = f'{ERA5_Press_path}/u/{year}/u_era5_oper_pl_{year}{month:02d}01-{year}{month:02d}{dom[month]:02d}.nc'   # wind speed
    ERA5_Press_V_file  = f'{ERA5_Press_path}/v/{year}/v_era5_oper_pl_{year}{month:02d}01-{year}{month:02d}{dom[month]:02d}.nc'   # wind speed
    ERA5_Press_PV_file = f'{ERA5_Press_path}/pv/{year}/pv_era5_oper_pl_{year}{month:02d}01-{year}{month:02d}{dom[month]:02d}.nc' # Potential vorticity
    ERA5_Press_w_file  = f'{ERA5_Press_path}/w/{year}/w_era5_oper_pl_{year}{month:02d}01-{year}{month:02d}{dom[month]:02d}.nc'   # Vertical velocity

    # Open the NetCDF file
    # -10: convert UTC to AEDT 
    t_s        = max(0,(day-1)*24-10)
    t_e        = day*24-10

    loc_lat    = [-45,-5]
    loc_lon    = [110,156]
    
    with nc.Dataset(ERA5_MSLP_T_file) as f_T_mslp:
        lat       = f_T_mslp.variables['latitude'][:]
        lon       = f_T_mslp.variables['longitude'][:]
        lat_idx_s = np.argmin(np.abs(lat - loc_lat[0]))
        lat_idx_e = np.argmin(np.abs(lat - loc_lat[1]))
        lon_idx_s = np.argmin(np.abs(lon - loc_lon[0]))
        lon_idx_e = np.argmin(np.abs(lon - loc_lon[1]))
        T_mslp    = np.nanmean(f_T_mslp.variables['t2m'][t_s:t_e,lat_idx_s:lat_idx_e:-1,lon_idx_s:lon_idx_e], axis=0)
        T_mslp    = T_mslp - 273.15
        lat_AU    = f_T_mslp.variables['latitude'][lat_idx_s:lat_idx_e:-1]
        lon_AU    = f_T_mslp.variables['longitude'][lon_idx_s:lon_idx_e]

    with nc.Dataset(ERA5_MSLP_P_file) as f_P_mslp:
        P_mslp = np.nanmean(f_P_mslp.variables['sp'][t_s:t_e,lat_idx_s:lat_idx_e:-1,lon_idx_s:lon_idx_e], axis=0)

    with nc.Dataset(ERA5_MSLP_U_file) as f_U_mslp:
        U_mslp = np.nanmean(f_U_mslp.variables['u10'][t_s:t_e,lat_idx_s:lat_idx_e:-1,lon_idx_s:lon_idx_e], axis=0)

    with nc.Dataset(ERA5_MSLP_V_file) as f_V_mslp:
        V_mslp = np.nanmean(f_V_mslp.variables['v10'][t_s:t_e,lat_idx_s:lat_idx_e:-1,lon_idx_s:lon_idx_e], axis=0)

    with nc.Dataset(ERA5_MSLP_R_file) as f_R_mslp:
        R_mslp = np.nansum(f_R_mslp.variables['tp'][t_s:t_e,lat_idx_s:lat_idx_e:-1,lon_idx_s:lon_idx_e], axis=0)

    with nc.Dataset(ERA5_Press_T_file) as f_T_press:
        T_press = np.nanmean(f_T_press.variables['t'][t_s:t_e,:,lat_idx_s:lat_idx_e:-1,lon_idx_s:lon_idx_e], axis=0)
        T_press = T_press - 273.15

    with nc.Dataset(ERA5_Press_Z_file) as f_Z_press:
        Z_press = np.nanmean(f_Z_press.variables['z'][t_s:t_e,:,lat_idx_s:lat_idx_e:-1,lon_idx_s:lon_idx_e], axis=0)
        #geopotential to geopotential height
        Z_press = Z_press/9.81

    with nc.Dataset(ERA5_Press_U_file) as f_U_press:
        U_press = np.nanmean(f_U_press.variables['u'][t_s:t_e,:,lat_idx_s:lat_idx_e:-1,lon_idx_s:lon_idx_e], axis=0)

    with nc.Dataset(ERA5_Press_V_file) as f_V_press:
        V_press = np.nanmean(f_V_press.variables['v'][t_s:t_e,:,lat_idx_s:lat_idx_e:-1,lon_idx_s:lon_idx_e], axis=0)

    with nc.Dataset(ERA5_Press_PV_file) as f_PV_press:
        PV_press = np.nanmean(f_PV_press.variables['pv'][t_s:t_e,:,lat_idx_s:lat_idx_e:-1,lon_idx_s:lon_idx_e], axis=0)

    with nc.Dataset(ERA5_Press_w_file) as f_w_press:
        w_press = np.nanmean(f_w_press.variables['w'][t_s:t_e,:,lat_idx_s:lat_idx_e:-1,lon_idx_s:lon_idx_e], axis=0)

    # ================== Start Plotting =================
    fig, axs = plt.subplots(nrows=2, ncols=4, figsize=[16,8],sharex=False,
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
    texts = ["(a)","(b)","(c)",
             "(d)","(e)","(f)",
             "(g)","(h)","(i)"]

    # Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                         facecolor="none",
                                         name="admin_1_states_provinces_shp")
                    #  925 850 500 300 hPa
    pressure_levels = [33, 30, 21, 17]  # Pressure levels for the plots
    variables = {"geopotential": Z_press, "wind_u": U_press, "wind_v": V_press, "w": w_press, "temperature": T_press}

    for i, ax in enumerate(axs.flat):

        ax.add_feature(states, linewidth=.5, edgecolor="black")
        ax.coastlines('50m', linewidth=0.8)
    
        # Set map extent for specific region
        ax.set_extent([110, 156, -45, -5])
        ax.coastlines(resolution="50m", linewidth=1)

        if i == 0:  # SPI 
            # Plot the original accumulated day data
            accumulated_day               = np.where(accumulated_day<0, np.nan, accumulated_day)
            accumulated_day_large_drought = np.where(accumulated_day_large_drought<0, np.nan, accumulated_day_large_drought)
            im1                           = ax.imshow(accumulated_day_large_drought[:, :], vmin=0, vmax=200, cmap='hot_r', 
                                                      extent=[lon_AGCD.min(), lon_AGCD.max(), lat_AGCD.min(), lat_AGCD.max()])
            cbar1                         = plt.colorbar(im1, ax=ax, orientation="horizontal", pad=0.02, aspect=20, shrink=0.7)
            cbar1.set_label('days')
            ax.text(0.02, 0.95, texts[i], transform=ax.transAxes, fontsize=12, bbox=props)
        
        elif i == 1: # AGCD rainfall
            precip_day                    = np.where(precip_day<0.0001, np.nan, precip_day)
            im2                           = ax.imshow(precip_day[:, :], vmin=0, vmax=20, cmap='YlGnBu', 
                                                      extent=[lon_AGCD.min(), lon_AGCD.max(), lat_AGCD.min(), lat_AGCD.max()])
            cbar2                         = plt.colorbar(im2, ax=ax, orientation="horizontal", pad=0.02, aspect=20, shrink=0.7)
            cbar2.set_label('rainfall')
            ax.text(0.02, 0.95, texts[i], transform=ax.transAxes, fontsize=12, bbox=props)
        
        elif i == 2: # 2m temperature
            # Temperature as contourf
            level_temp = [-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55]
            pcm  = ax.contourf(lon_AU, lat_AU, T_mslp, levels=level_temp, cmap="Spectral_r", extend="max")
            cbar = fig.colorbar(pcm, ax=ax, orientation="horizontal", pad=0.05, aspect=20, shrink=0.7)
            cbar.set_label("Temperature (C)")        

            # MSLP as contour lines
            c = ax.contour(lon_AU, lat_AU, P_mslp / 100, levels=np.arange(980, 1050, 10), colors="black", linewidths=0.8) # 
            ax.clabel(c, inline=True, fontsize=8, fmt='%d hPa')

            # 10m Wind as arrows
            ax.quiver(lon_AU[::5], lat_AU[::5], U_mslp[::5, ::5], V_mslp[::5, ::5],  scale=500, width=0.002, color="black")
            ax.text(0.02, 0.95, texts[i], transform=ax.transAxes, fontsize=12, bbox=props)

        elif i == 3: # MSLP
            # Precipitation as contourf
            R_mslp = np.where(R_mslp*1000<0.01, np.nan,R_mslp)
            pcm  = ax.contourf(lon_AU, lat_AU, R_mslp * 1000, levels=np.arange(0, 30, 2), cmap="YlGnBu", extend="max")
            cbar = fig.colorbar(pcm, ax=ax, orientation="horizontal", pad=0.05, aspect=20, shrink=0.7)
            cbar.set_label("Precipitation (mm/day)")     

            # Smooth the MSLP field using a moving average (uniform filter)
            P_mslp_smoothed = uniform_filter(P_mslp, size=5)  # Adjust size for the desired smoothness

            # MSLP as contour lines
            c = ax.contour(lon_AU, lat_AU, P_mslp_smoothed / 100, levels=np.arange(980, 1050, 10), colors="black", linewidths=0.8) # 
            ax.clabel(c, inline=True, fontsize=8, fmt='%d hPa')

            # 10m Wind as arrows
            ax.quiver(lon_AU[::5], lat_AU[::5], U_mslp[::5, ::5], V_mslp[::5, ::5],  scale=500, width=0.002, color="black")
            ax.text(0.02, 0.95, texts[i], transform=ax.transAxes, fontsize=12, bbox=props)

        else:  # Pressure levels
            level_idx = pressure_levels[i - 4]  # Select pressure level (925, 850, or 500 hPa)

            # # Temperature as contourf
            # level_temp = [-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55]
            # pcm = ax.contourf(lon_AU, lat_AU, variables["temperature"][level_idx], levels=level_temp, cmap="hot_r", extend="both")
            # cbar = fig.colorbar(pcm, ax=ax, orientation="horizontal", pad=0.05, aspect=20, shrink=0.7)
            # cbar.set_label("Temperature (K)")

            # Vertical velocity as red contours
            level_velo = [-1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.,1.2]
            pcm = ax.contourf(lon_AU, lat_AU, variables["w"][level_idx], levels=level_velo, cmap="bwr", extend="both")
            cbar = fig.colorbar(pcm, ax=ax, orientation="horizontal", pad=0.05, aspect=20, shrink=0.7)
            cbar.set_label("Vertical Velocity (Pa/s)")

            # Geopotential as contour lines
            c = ax.contour(lon_AU, lat_AU, variables["geopotential"][level_idx]/10., levels=np.arange(0, 600, 5), colors="black", linewidths=0.8)
            ax.clabel(c, inline=True, fontsize=8, fmt='%d x 10 gpm')

            # Wind as arrows
            ax.quiver(lon_AU[::5], lat_AU[::5], variables["wind_u"][level_idx, ::5, ::5], variables["wind_v"][level_idx, ::5, ::5], scale=300, width=0.002, color="black")

            # # Vertical velocity as red contours
            # c = ax.contour(lon_AU, lat_AU, variables["w"][level_idx], levels=np.arange(-1, 1.1, 0.4), colors="white", linewidths=0.8)
            # ax.clabel(c, inline=True, fontsize=8, fmt='%0.1f Pa/s')
            ax.text(0.02, 0.95, texts[i], transform=ax.transAxes, fontsize=12, bbox=props)

    # Add titles and save the figure    
    axs[0, 0].set_title("Accumulated Days in Drought")
    axs[0, 1].set_title("Rainfall")
    axs[0, 2].set_title("Temperature")
    axs[0, 3].set_title("MSLP & Precipitation")
    axs[1, 0].set_title("925 hPa Level")
    axs[1, 1].set_title("850 hPa Level")
    axs[1, 2].set_title("500 hPa Level")
    axs[1, 3].set_title("300 hPa Level")
    
    output_path = f"/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/plots/{year}"
    if not os.path.exists(output_path):  # Check if the path does NOT exist
        os.mkdir(output_path)  # Create the directory
    plt.tight_layout()
    plt.savefig(f"{output_path}/Diagnose_weather_system_{time}" , dpi=300, bbox_inches="tight")

    return

def keep_only_large_drought_region(i, time, accumulated_day, precip_day, lat, lon, min_drought_area_km2=10000):
    
    # Define the drought threshold
    drought_binary = accumulated_day > 0

    # Label connected regions in the binary drought mask
    labeled_droughts, num_features = label(drought_binary)

    # Calculate grid cell size (resolution) based on lat/lon coordinates
    lat_resolution = abs(lat[1] - lat[0])  # Latitude resolution in degrees
    lon_resolution = abs(lon[1] - lon[0])  # Longitude resolution in degrees

    # Approximate cell area in km^2
    earth_radius_km = 6371.0
    cell_area_km2 = (
        (lat_resolution * np.pi / 180) *
        (lon_resolution * np.pi / 180) *
        (earth_radius_km ** 2) )

    # Calculate area in km^2 for each labeled region
    drought_areas     = nd_sum(drought_binary, labeled_droughts, index=range(1, num_features + 1))
    drought_areas_km2 = drought_areas * cell_area_km2

    # Filter regions that meet the minimum area threshold
    large_drought_labels = [i + 1 for i, area in enumerate(drought_areas_km2) if area >= min_drought_area_km2]
    large_drought_map    = np.isin(labeled_droughts, large_drought_labels)
    
    accumulated_day_large_drought = np.where(~np.isnan(large_drought_map), accumulated_day, np.nan)
    
    if large_drought_labels != None:
        
        # ============ Plotting weather system diagnose =============
        plot_diagnose_weather_system(time, lat, lon, accumulated_day, accumulated_day_large_drought, precip_day)
    
    large_drought_labels= None
    large_drought_map   = None
    precip_day          = None
    
    return 

def plot_disgnose_weather_system_multi_days():

    input_file = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/spi_pearson_90_reorder_nan_filled_drought_periods.nc'
    mask_file  = '/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_DLCM_mask.nc'
    rain_file  = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files/agcd_v1_precip_total_r005_daily_1950_2023.nc'

    with nc.Dataset(input_file, mode='r') as f_in:
        lat_out          = f_in.variables['lat'][:]
        lon_out          = f_in.variables['lon'][:]
        accumulated_days = f_in.variables['accumulated_days'][:,:,:]
        time             = nc.num2date(f_in.variables['time'][:],f_in.variables['time'].units,
                            only_use_cftime_datetimes=False,only_use_python_datetimes=True)
        ntime            = len(time)
        
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
    accumulated_days     = np.where(landsea_3d==0, accumulated_days, np.nan)
    precip               = np.where(landsea_3d==0, precip, np.nan)

    min_drought_area_km2 = 10000  # Minimum size in square km

    # Prepare arguments for parallel processing
    args_list = [ (i, time[i], accumulated_days[i,::-1,:], precip[i,::-1,:], 
                    lat_out, lon_out, min_drought_area_km2) 
                    for i in np.arange(18263, ntime)]

    # Use multiprocessing to calculate accumulated days for each grid point
    with Pool() as pool:
        pool.starmap(keep_only_large_drought_region, args_list)

    return

# Main script
if __name__ == "__main__":

    plot_disgnose_weather_system_multi_days()