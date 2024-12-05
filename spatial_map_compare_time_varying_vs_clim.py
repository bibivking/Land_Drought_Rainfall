#!/usr/bin/python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"

'''
Functions:
1. Compare LIS-CABLE with GRACE, GLEAM, & DOLCE
2. GW vs FD
3. plot time-series and spitial (difference) map
'''
import glob
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.interpolate import griddata
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from common_utils import *

def calc_cable_tws(file_path, loc_lat, loc_lon, lat_name, lon_name):

    # calculate TWS for CABLE simulations

    print("calc_cable_tws")

    Soil_thickness = [0.022, 0.058, 0.154, 0.409, 1.085, 2.872]
    off_file       = Dataset(file_path, mode='r')
    GWdz           = off_file.variables['GWdz'][:]
    SoilMoist      = off_file.variables['SoilMoist'][:]

    time1, CanopInt= read_var(file_path, "CanopInt", loc_lat, loc_lon, lat_name, lon_name)

    time2, GWMoist = read_var(file_path, "GWMoist", loc_lat, loc_lon, lat_name, lon_name)
    GWMoist        = GWMoist*GWdz*1000.

    time3, SWE     = read_var(file_path, "SWE", loc_lat, loc_lon, lat_name, lon_name)
    print('np.shape(SWE)')
    print(np.shape(SWE))

    TWS            = GWMoist + CanopInt + SWE
    print('np.shape(TWS)')
    print(np.shape(TWS))

    for i in np.arange(6):
        TWS = TWS + SoilMoist[:,i,:,:]*Soil_thickness[i]*1000

    year_2004       = datetime(2004,1,1)
    year_2009       = datetime(2009,12,31)

    TWS_0409 = spatial_var(time1,TWS,year_2004,year_2009)
    print('np.shape(TWS_0409)')
    print(np.shape(TWS_0409))
    TWS      = TWS - TWS_0409
    print('np.shape(TWS)')
    print(np.shape(TWS))
    return time1,TWS

def plot_spital_map( yr, file_paths, var_names, time_s, time_e, file_paths2=None, loc_lat=None, loc_lon=None,
                    lat_names=None, lon_names=None, message=None,soil_layer=None):

    '''
    Plot either value or difference
    '''
    print("======== In plot_spital_map =========")

    # Open the NetCDF4 file (add a directory path if necessary) for reading:

    time1, Var1     = read_var_multi_file(file_paths, var_names[0], loc_lat, loc_lon, lat_name=lat_names[0], lon_name=lon_names[0])
    time_tmp, lats1 = read_var(file_paths[0], lat_names[0], loc_lat, loc_lon, lat_name=lat_names[0], lon_name=lon_names[0])
    time_tmp, lons1 = read_var(file_paths[0], lon_names[0], loc_lat, loc_lon, lat_name=lat_names[0], lon_name=lon_names[0])

    scale           = get_scale(var_names[0])
    var1            = spatial_var_mean(time1,Var1,time_s,time_e)*scale

    if file_paths2 is not None:

        time2, Var2     = read_var_multi_file(file_paths2, var_names[1], loc_lat=None, loc_lon=None, lat_name=lat_names[1], lon_name=lon_names[1])
        time_tmp, lats2 = read_var(file_paths2[0], lat_names[1], loc_lat=None, loc_lon=None, lat_name=lat_names[1], lon_name=lon_names[1])
        time_tmp, lons2 = read_var(file_paths2[0], lon_names[1], loc_lat=None, loc_lon=None, lat_name=lat_names[1], lon_name=lon_names[1])

        scale           = get_scale(var_names[1])

        if var_names[1] in ['E','Et','Ei','Es','etot']:
            var2        = np.nansum(Var2,axis=0)
        else:
            var2        = np.nanmean(Var2,axis=0)*scale

        if soil_layer != None:
            if np.all(lats1 == lats2) and np.all(lons1 == lons2):
                var          = var1[soil_layer,:,:] - var2[soil_layer,:,:]
            else:
                var2_regrid  = regrid_data(lats2, lons2, lats1, lons1, var2[soil_layer,:,:], method='nearest',threshold=None)
                var          = var1[soil_layer,:,:] - var2_regrid
        else:
            if np.all(lats1 == lats2) and np.all(lons1 == lons2):
                var             = var1 - var2
            else:
                var2_regrid     = regrid_data(lats2, lons2, lats1, lons1, var2, method='nearest',threshold=None)
                var             = var1 - var2_regrid
    else:
        if soil_layer != None:
            var             = var1[soil_layer,:,:]
        else:
            var             = var1

    # Create a figure with 3 subplots in 1 row
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(18, 5), subplot_kw={'projection': ccrs.PlateCarree()})

    # Adjust spacing between subplots
    plt.subplots_adjust(wspace=0.3)

    # Loop through each axis
    for ax in axs:

        # Set extent based on loc_lat and loc_lon
        if loc_lat is None:
            ax.set_extent([140, 154, -40, -28])  # Example extent, adjust as needed
        else:
            ax.set_extent([loc_lon[0], loc_lon[1], loc_lat[0], loc_lat[1]])

        ax.coastlines(resolution="50m", linewidth=1)

        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='black', linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = True

        if loc_lat is None:
            gl.xlocator = mticker.FixedLocator([140, 145, 150])
            gl.ylocator = mticker.FixedLocator([-40, -35, -30])
        else:
            gl.xlocator = mticker.FixedLocator(loc_lon)
            gl.ylocator = mticker.FixedLocator(loc_lat)

        gl.xformatter = LongitudeFormatter()
        gl.yformatter = LatitudeFormatter()
        gl.xlabel_style = {'size': 10, 'color': 'black'}
        gl.ylabel_style = {'size': 10, 'color': 'black'}

    # Plot windspeed based on variable names
    if var_names[0] in ['SoilMoist', 'GWMoist']:
        clevs      = np.linspace(0, 0.400, num=21)
        clevs_diff = np.linspace(-0.200, 0.200, num=21)
    elif var_names[0] == 'Evap':
        clevs      = np.linspace(0., 1000., num=21)
        clevs_diff = np.linspace(-200., 200., num=21)
    elif var_names[0] == 'Rainf':
        var1       = var1 * 365
        var2       = var2 * 365
        var        = var  * 365
        clevs      = np.linspace(0., 2000., num=21)
        clevs_diff = np.linspace(-200., 200., num=21)
    elif var_names[0] in ['Qh', 'Qle']:
        clevs      = np.linspace(-100., 800., num=19)
        clevs_diff = np.linspace(-200., 200., num=21)
    elif var_names[0] == 'LAI':
        clevs      = np.linspace(0., 7., num=21)
        clevs_diff = np.linspace(-2., 2., num=21)
    elif var_names[0] == 'Albedo':
        clevs      = np.linspace(0., 0.5, num=21)
        clevs_diff = np.linspace(-0.02, 0.02, num=21)
        
    if soil_layer != None:
        plot1 = axs[0].contourf(lons1, lats1, var1[soil_layer,:,:], levels=clevs, transform=ccrs.PlateCarree(), extend='both', cmap=plt.cm.BrBG)
        plot2 = axs[1].contourf(lons2, lats2, var2[soil_layer,:,:], levels=clevs, transform=ccrs.PlateCarree(), extend='both', cmap=plt.cm.BrBG)
        plot3 = axs[2].contourf(lons1, lats1, var,  levels=clevs_diff, transform=ccrs.PlateCarree(), extend='both', cmap=plt.cm.BrBG)
    else:
        plot1 = axs[0].contourf(lons1, lats1, var1, levels=clevs, transform=ccrs.PlateCarree(), extend='both', cmap=plt.cm.BrBG)
        plot2 = axs[1].contourf(lons2, lats2, var2, levels=clevs, transform=ccrs.PlateCarree(), extend='both', cmap=plt.cm.BrBG)
        plot3 = axs[2].contourf(lons1, lats1, var,  levels=clevs_diff, transform=ccrs.PlateCarree(), extend='both', cmap=plt.cm.BrBG)

    axs[0].set_title(var_names[0], size=16)
    axs[1].set_title(var_names[0], size=16)
    axs[2].set_title(var_names[0], size=16)

    cb = plt.colorbar(plot1, ax=axs[0], orientation="horizontal", pad=0.02, aspect=16, shrink=0.8)
    cb = plt.colorbar(plot2, ax=axs[1], orientation="horizontal", pad=0.02, aspect=16, shrink=0.8)
    cb = plt.colorbar(plot3, ax=axs[2], orientation="horizontal", pad=0.02, aspect=16, shrink=0.8)
    cb.ax.tick_params(labelsize=10)

    # Construct message for saving
    if message is None:
        message = var_names[0]
    else:
        message = message + "_" + var_names[0]

    plt.savefig('./plots/spatial_map_obs_'+message+'.png',dpi=300)

def plot_time_series( file_paths, var_names, year_s, year_e, loc_lat=None, loc_lon=None,
                      lat_names=None, lon_names=None, message=None ):

    print("======== In plot_time_series =========")

    fig, ax = plt.subplots()

    # plot line 1
    if var_names[0] == "TWS":
        Time1, Var1 = calc_cable_tws(file_paths[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
        time1, var1 = time_series_var(Time1, Var1, year_s, year_e)
    elif var_names[0] == "lwe_thickness":
        Time1, Var1 = read_var(file_paths[0], var_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
        time1, var1 = time_series_var(Time1, Var1*10., year_s, year_e)
    else:
        Time1, Var1 = read_var(file_paths[0], var_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
        Var1        = np.where(Var1<=-9999, np.nan, Var1)
        time1, var1 = time_series_var(Time1, Var1, year_s, year_e) # np.cumsum()

    t1 = []
    for i in np.arange(len(time1)):
        t1.append(time1[i].days)
    if var_names[0] in ['Qs','Qsb','Rainf','Evap','ESoil','ECanop','TVeg']:
        scale = 24.*3600.
    else:
        scale = 1.
    print("var1*scale")
    print(var1*scale)
    ax.plot(t1, var1*scale, c = "blue", label=var_names[0], alpha=0.5)

    # plot line 2
    if len(file_paths) > 1:
        if var_names[1] == "TWS":
            Time2, Var2 = calc_cable_tws(file_paths[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
            time2, var2 = time_series_var(Time2, Var2, year_s, year_e)
        elif var_names[1] == "lwe_thickness":
            Time2, Var2 = read_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
            time2, var2 = time_series_var(Time2, Var2*10., year_s, year_e)
        else:
            Time2, Var2 = read_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
            time2, var2 = time_series_var(Time2, Var2, year_s, year_e) #np.cumsum()
        t2 = []
        for i in np.arange(len(time2)):
            t2.append(time2[i].days)
        if var_names[1] in ['Qs','Qsb','Rainf','Evap','ESoil','ECanop','TVeg']:
            scale = 24.*3600.*31 #MMY monthly
        else:
            scale = 1.

        # plt.contourf(Var2[0,:,:])
        # plt.show()

        print("var2*scale")
        print(var2*scale)
        ax.plot(t2, var2*scale, c = "green", label=var_names[1], alpha=0.5)

    # plot line 3
    if len(file_paths) > 2:
        if var_names[2] == "TWS":
            Time3, Var3 = calc_cable_tws(file_paths[2], loc_lat, loc_lon, lat_names[2], lon_names[2])
            time3, var3 = time_series_var(Time3, Var3, year_s, year_e)
        elif var_names[2] == "lwe_thickness":
            Time3, Var3 = read_var(file_paths[2], var_names[2], loc_lat, loc_lon, lat_names[2], lon_names[2])
            time3, var3 = time_series_var(Time3, Var3*10., year_s, year_e)
        else:
            Time3, Var3 = read_var(file_paths[2], var_names[2], loc_lat, loc_lon, lat_names[2], lon_names[2])
            time3, var3 = time_series_var(Time3, Var3, year_s, year_e) # np.cumsum()
        t3 = []
        for i in np.arange(len(time3)):
            t3.append(time3[i].days)
        if var_names[2] in ['Qs','Qsb','Rainf','Evap','ESoil','ECanop','TVeg']:
            scale = 24.*3600.*31 #MMY monthly
        else:
            scale = 1.

        print("var3*scale")
        print(var3*scale)
        ax.plot(t3, var3*scale, c = "red", label=var_names[2], alpha=0.5)

    # ax.set_xlim([np.min(var1*scale,var2*scale), np.max(var1*scale,var2*scale)])
    # Time2, Var2 = read_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_name[1], lon_name[1])
    # time2, var2 = time_series_var(Time2,Var2,year_s,year_e)
    # var = np.zeros((2,len(var1)))
    # var[0,:] = var1
    # var[1,:] = var2
    # ax.plot(t1, var*scale, alpha=0.5)
    # ax.set_ylabel('mm')
    ax.set_title(var_names[0])
    # ax.set_xticks(x1[::1440])
    # ax.set_xticklabels(np.arange(year_s,year_e,1))
    ax.legend()

    fig.tight_layout()
    if message == None:
        message = var_names[0]
    else:
        message = message + "_" + var_names[0]
    if loc_lat != None:
        message = message + "_lat="+str(loc_lat[0]) +"-"+str(loc_lat[1]) + "_lon="+str(loc_lon[0])+"-"+str(loc_lon[1])

    plt.savefig('./plots/time_series_'+message+'.png',dpi=300)

if __name__ == "__main__":

    # #######################
    #     path setting      #
    # #######################

    DOLCE_path   = "/g/data/w97/mm3972/data/DOLCE/v3/"
    DOLCE_file   = DOLCE_path + "DOLCE_v3_2000-2018.nc"
    GLEAM_path   = "/g/data/w97/Shared_data/Observations/Global_ET_products/GLEAM_v3_5/v3-6a/monthly/"
    GLEAM_file   = GLEAM_path + "E_1980-2021_GLEAM_v3.6a_MO.nc"
    GRACE_path   = "/g/data/w97/mm3972/data/GRACE/GRACE_JPL_RL06/GRACE_JPLRL06M_MASCON/"
    GRACE_file   = GRACE_path + "GRCTellus.JPL.200204_202004.GLO.RL06M.MSCNv02CRI.nc"
    AWRA_ET_path = '/g/data/w97/mm3972/data/Soil_Moisture/AWRA_v6/etot/'

    ERA5_rain_path = '/g/data/rt52/era5/single-levels/monthly-averaged/tp'

    SP_off_path = "/g/data/w97/mm3972/model/cable/runs/Land_drought_rainfall_runs/spinup_run_1970_1999/outputs/"
    SP_off_file = SP_off_path + "cable_out_1970-1999.nc"
    SP_off_clim_file = SP_off_path + "cable_out_1970-1999_mean.nc"

    SP_GW_reduce10_path = "/g/data/w97/mm3972/model/cable/runs/Land_drought_rainfall_runs/spinup_run_1970_1999_GWMoist_reduce_10percent/outputs/"
    SP_GW_reduce30_path = "/g/data/w97/mm3972/model/cable/runs/Land_drought_rainfall_runs/spinup_run_1970_1999_GWMoist_reduce_30percent/outputs/"
    SP_GW_reduce50_path = "/g/data/w97/mm3972/model/cable/runs/Land_drought_rainfall_runs/spinup_run_1970_1999_GWMoist_reduce_50percent/outputs/"
    
    SP_GW_reduce50_file = SP_GW_reduce50_path + "cable_out_1970-1999.nc"
    
    CABLE_off_path = "/g/data/w97/mm3972/model/cable/runs/Land_drought_rainfall_runs/run_2000_2019_daily_rst_control/outputs/"
    CABLE_off_file = CABLE_off_path + "cable_out_2000-2019.nc"


    # #######################
    #   plot_spital_map     #
    # #######################

    if 1:

        ### Compare CABLE 2000-2019 with SP clim
        
        '''
        plot OFF vs SP clim
        '''
        
        for yr in np.arange(2000, 2020):
            year_s      = datetime(yr,1,1)
            year_e      = datetime(yr,12,31)
            loc_lat     = [-45,-5]
            loc_lon     = [112,154]
        
            #### plot Evap new SP vs GLEAM ###
            print("plot OFF vs SP clim")
            file_paths1 = [f"{CABLE_off_path}cable_out_{yr}.nc"]
            file_paths2 = [SP_off_clim_file]
            lat_names   = ["latitude","latitude"]#"lat"
            lon_names   = ["longitude","longitude"]#"lon"
    
            message     = f"OFF_LAI_SG-SP_Clim_{yr}"
            
            # var_name    = ["Evap","Evap"]
            # plot_spital_map( yr, file_paths1, var_name, year_s, year_e, file_paths2=file_paths2, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
            #                 lon_names=lon_names,message=message)
        
            var_name    = ["Qh","Qh"]
            plot_spital_map( yr, file_paths1, var_name, year_s, year_e, file_paths2=file_paths2, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
                            lon_names=lon_names,message=message)
        
            var_name    = ["Qle","Qle"]
            plot_spital_map( yr, file_paths1, var_name, year_s, year_e, file_paths2=file_paths2, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
                            lon_names=lon_names,message=message)

            var_name    = ["LAI","LAI"]
            plot_spital_map( yr, file_paths1, var_name, year_s, year_e, file_paths2=file_paths2, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
                            lon_names=lon_names,message=message)

            var_name    = ["Albedo","Albedo"]
            plot_spital_map( yr, file_paths1, var_name, year_s, year_e, file_paths2=file_paths2, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
                            lon_names=lon_names,message=message)

            var_name    = ["SoilMoist","SoilMoist"]
            for soil_layer in np.arange(6):
                message     = f"OFF_LAI_SG-SP_Clim_{yr}_layer={soil_layer}"
                plot_spital_map( yr, file_paths1, var_name, year_s, year_e, file_paths2=file_paths2, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
                                lon_names=lon_names,message=message,soil_layer=soil_layer)