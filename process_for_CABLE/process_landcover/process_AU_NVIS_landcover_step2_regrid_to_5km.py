import os
import argparse
import netCDF4 as nc
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from datetime import datetime, timedelta
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import matplotlib.ticker as mticker

def find_dominant_lc_type(lc_coarse):

    # Flatten the array
    lc_coarse_flat = lc_coarse.flatten()

    # Find the most frequent value
    dominant_lc_index, dominant_lc_count = stats.mode(lc_coarse_flat)

    # print("Most frequent value:", dominant_lc_index)
    # print("Count of most frequent value:", dominant_lc_count)

    return dominant_lc_index

def save_land_cover_to_ncfile(lc_coarse,lat_in,lon_in,lc_dyn_clim="time_varying",lc_type="PFT",year=None):

    if lc_dyn_clim=="time_varying":
        output_file = f'/g/data/w97/mm3972/data/MODIS/MODIS_landcover/var_landcover_only/MCD12Q1.061_landcover_{year}_{lc_type}_5km.nc'
    elif lc_dyn_clim=="climatology":
        output_file = f'/g/data/w97/mm3972/data/MODIS/MODIS_landcover/var_landcover_only/MCD12Q1.061_landcover_clim_{lc_type}_5km.nc'

    # set dimension
    nlat = len(lat_in)
    nlon = len(lon_in)

    # Create new nc file
    f_out               = nc.Dataset(output_file, 'w', format='NETCDF4')
    f_out.history_mmy   = "Created by: %s" % (os.path.basename(__file__))
    f_out.creation_date = "%s" % (datetime.now())

    # set dimensions
    f_out.createDimension('lat', nlat)
    f_out.createDimension('lon', nlon)
    f_out.Conventions   = "CF-1.0"
    f_out.title         = "MCD12Q1.061 for aid0001"
    f_out.institution   = "Land Processes Distributed Active Archive Center (LP DAAC)"
    f_out.source        = "AppEEARS v3.60"
    f_out.history_old   = "Mon Aug 19 11:01:39 2024: ncatted -O -a calendar,time,m,c,gregorian MCD12Q1.061_500m_aid0001_2021-2022.nc\nSee README.md"
    f_out.NCO           = "netCDF Operators version 5.1.4 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)"

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

    var_out                = f_out.createVariable('LC', 'i2', ('lat', 'lon'), fill_value=255)
    var_out.long_name      = "Land cover type for CABLE"
    var_out.coordinates    = "time lat lon"
    var_out.grid_mapping   = "crs"
    var_out.description    =(   "1: evergreen_needleleaf forest; "
                                "2: evergreen_broadleaf forest; "
                                "3: deciduous_needleleaf forest; "
                                "4: deciduous_broadleaf forest; "
                                "5: shrub; "
                                "6: C3 grassland; "
                                "7: C4 grassland; "
                                "8: Tundra grass; "
                                "9: C3 cropland; "
                                "10: C4 cropland; "
                                "11: wetland; "
                                "12: empty (grass); "
                                "13: empty (not used); "
                                "14: barren no veg; "
                                "15: urban no veg; "
                                "16: lakes no veg; "
                                "17: ice no veg"  )
    var_out[:]             = lc_coarse
    f_out.close()

    lat_out = None
    lon_out = None
    var_out = None

    return

def recast_NVIS_to_CABLE_PFT(lc_fine_res_tmp, lc_type):

    '''
    In  National Vegetation Information System (NVIS) classification:
    1: Rainforests and Vine Thickets, 2:Eucalypt Tall Open Forests, 3: Eucalypt Open Forests
    4: Eucalypt Low Open Forests, 5: Eucalypt Woodlands, 6: Acacia Forests and Woodlands,
    7: Callitris Forests and Woodlands, 8: Casuarina Forests and Woodlands, 9: Melaleuca Forests and Woodlands,
    10: Other Forests and Woodlands, 11: Eucalypt Open Woodlands, 12: Tropical Eucalypt Woodlands/Grasslands,
    13: Acacia Open Woodlands, 14: Mallee Woodlands and Shrublands, 15: Low Closed Forests and Tall Closed Shrublands,
    16: Acacia Shrublands, 17: Other Shrublands, 18: Heathlands, 19: Tussock Grasslands, 20: Hummock Grasslands, 
    21: Other Grasslands, Herblands, Sedgelands and Rushlands, 22: Chenopod Shrublands, Samphire Shrublands and Forblands
    23: Mangroves, 24: Inland Aquatic - freshwater, salt lakes, lagoons, 25: Cleared, non-native vegetation, buildings, 
    26: Unclassified native vegetation, 27: Naturally bare - sand, rock, claypan, mudflat, 28: Sea and estuaries
    29: Regrowth, modified native vegetation, 30: Unclassified forest, 31: Other Open Woodlands, 
    32, Mallee Open Woodlands and Sparse Mallee Shrublands, 99: Unknown/no data

    NVIS => CABLE: no change: 1 = 1;  2 = 2;  3 = 3;  4 = 4; 11 = 11
                    change   : 5 => 2; 6, 7 => 5; 8 => 2; 9=> 6 or 7; 10 => 6 or 7;
                                12 => 9 or 10; 13 => 15;
                                14 => 9 (might be 1,2,3,4,5,6,7,10); 15 => 17
                                16 => 14, 17 => 16
    '''

    "1: evergreen_needleleaf forest; "
                                "2: evergreen_broadleaf forest; "
                                "3: deciduous_needleleaf forest; "
                                "4: deciduous_broadleaf forest; "
                                "5: shrub; "
                                "6: C3 grassland; "
                                "7: C4 grassland; "
                                "8: Tundra grass; "
                                "9: C3 cropland; "
                                "10: C4 cropland; "
                                "11: wetland; "
                                "12: empty (grass); "
                                "13: empty (not used); "
                                "14: barren no veg; "
                                "15: urban no veg; "
                                "16: lakes no veg; "
                                "17: ice no veg"  )

    # Set C4 plants
    Ctype_path = "/g/data/w97/mm3972/data/AU_C3_C4_plant_map/nc_file/C3_or_C4_dominate_map_NVIS_res.nc"
    f_ctype    = nc.Dataset(Ctype_path,'r')
    ctype      = f_ctype.variables['Ctype'][:,:]
    f_ctype.close()

    lc_coarse  = np.where((lc_coarse==6) & (ctype==4), 7, lc_coarse)
    lc_coarse  = np.where((lc_coarse==9) & (ctype==4), 10, lc_coarse)

    lc_fine_res = np.where(lc_fine_res_tmp == 5, 2,  lc_fine_res_tmp)
    lc_fine_res = np.where(lc_fine_res_tmp == 6, 5,  lc_fine_res)
    lc_fine_res = np.where(lc_fine_res_tmp == 7, 5,  lc_fine_res)
    lc_fine_res = np.where(lc_fine_res_tmp == 8, 2,  lc_fine_res)
    lc_fine_res = np.where(lc_fine_res_tmp == 9, 6,  lc_fine_res)
    lc_fine_res = np.where(lc_fine_res_tmp == 10, 6,  lc_fine_res)
    lc_fine_res = np.where(lc_fine_res_tmp == 12, 9,  lc_fine_res)
    lc_fine_res = np.where(lc_fine_res_tmp == 14,  9, lc_fine_res) ## Is it right?
    lc_fine_res = np.where(lc_fine_res_tmp == 16, 14, lc_fine_res)
    lc_fine_res = np.where(lc_fine_res_tmp == 17, 16, lc_fine_res)
    lc_fine_res = np.where(lc_fine_res_tmp == 15, 17, lc_fine_res)
    lc_fine_res = np.where(lc_fine_res_tmp == 13, 15, lc_fine_res)



    return lc_fine_res

def regrid_lc_to_coarse_res(lc_dyn_clim="time_varying"):

    # File paths
    input_file       = f'/g/data/w97/mm3972/data/AU_NVIS_landcover/nc_file/AU_NVIS_landcover.nc'
    target_grid_file = '/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc'

    # Open input file and read landcover, lat and lon information
    f_in             = nc.Dataset(input_file,'r')
    lc_fine_res_tmp  = f_in.variables['LC'][:,:,:] #*scale_factor
    lat_fine         = f_in.variables['lat'][:]
    lon_fine         = f_in.variables['lon'][:]
    ntime            = len(lc_fine_res_tmp[:,0,0])

    # Translate input landcover types to CABLE landcover types
    lc_fine_res      = recast_NVIS_to_CABLE_PFT(lc_fine_res_tmp, lc_type)

    # ================== Start Plotting =================
    fig, axs = plt.subplots(nrows=2, ncols=1, squeeze=True)
    cf = axs[0].imshow(lc_fine_res_tmp[0], origin="lower", interpolation="none", vmin=0, vmax=17)
    cf = axs[1].imshow(lc_fine_res[0], origin="lower", interpolation="none", vmin=0, vmax=17)
    plt.savefig(f'./plots/Check_cast_{lc_type}.png',dpi=300)

    # Read target lat and lon information
    f_grid       = nc.Dataset(target_grid_file,'r')
    lat_grid     = f_grid.variables['latitude'][:]
    lon_grid     = f_grid.variables['longitude'][:]
    nlat_grid    = len(lat_grid)
    nlon_grid    = len(lon_grid)

    if lc_dyn_clim=="time_varying":
        lc_coarse = np.zeros((ntime, nlat_grid, nlon_grid))
    elif lc_dyn_clim=="climatology":
        lc_coarse = np.zeros((nlat_grid, nlon_grid))

    # Loop input lat and lon
    for i, lat in enumerate(lat_grid):
        for j, lon in enumerate(lon_grid):

            # fine-res pixels in the coarse grid cell
            lat_index_s = np.where((lat_fine > lat-0.025) & (lat_fine < lat+0.025))[0][0]
            lat_index_e = np.where((lat_fine > lat-0.025) & (lat_fine < lat+0.025))[0][-1]

            lon_index_s = np.where((lon_fine > lon-0.025) & (lon_fine < lon+0.025))[0][0]
            lon_index_e = np.where((lon_fine > lon-0.025) & (lon_fine < lon+0.025))[0][-1]

            # print(lat,lon,
            #       'lat_index_s',lat_fine[lat_index_s],'lat_index_e',lat_fine[lat_index_e],'lat_index_e + 1',lat_fine[lat_index_e+1],
            #       'lon_index_s',lon_fine[lon_index_s],'lon_index_e',lon_fine[lon_index_e],'lon_index_e + 1',lon_fine[lon_index_e+1])

            # find the domainant land cover type in the coarse grid cell
            if lc_dyn_clim=="time_varying":
                for t in np.arange(ntime):
                    lc_coarse_tmp      = lc_fine_res[t,lat_index_s:lat_index_e+1,lon_index_s:lon_index_e+1]
                    lc_coarse[t, i, j] = find_dominant_lc_type(lc_coarse_tmp)
            elif lc_dyn_clim=="climatology":
                lc_coarse_tmp   = lc_fine_res[:,lat_index_s:lat_index_e+1,lon_index_s:lon_index_e+1]
                lc_coarse[i, j] = find_dominant_lc_type(lc_coarse_tmp)


    # Create output file
    if lc_dyn_clim=="time_varying":
        for i, year in enumerate(np.arange(2001,2023,1)):
            save_land_cover_to_ncfile(lc_coarse[i,:,:], lat_grid, lon_grid, lc_dyn_clim, lc_type, year)
    elif lc_dyn_clim=="climatology":
        save_land_cover_to_ncfile(lc_coarse, lat_grid, lon_grid, lc_dyn_clim, lc_type)

    return

if __name__ == "__main__":

    lc_type     = "PFT" #"PFT"
    lc_dyn_clim = "time_varying" #"time_varying" # "climatology"
    regrid_lc_to_coarse_res(lc_type=lc_type, lc_dyn_clim=lc_dyn_clim)

    lc_type     = "PFT" #"PFT"
    lc_dyn_clim = "climatology" #"time_varying" # "climatology"
    regrid_lc_to_coarse_res(lc_type=lc_type, lc_dyn_clim=lc_dyn_clim)

    lc_type     = "IGBP" #"PFT"
    lc_dyn_clim = "time_varying" #"time_varying" # "climatology"
    regrid_lc_to_coarse_res(lc_type=lc_type, lc_dyn_clim=lc_dyn_clim)

    lc_type     = "IGBP" #"PFT"
    lc_dyn_clim = "climatology" #"time_varying" # "climatology"
    regrid_lc_to_coarse_res(lc_type=lc_type, lc_dyn_clim=lc_dyn_clim)