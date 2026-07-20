import argparse
import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata

def is_leap(year):
    """Check if a year is a leap year."""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)

def regrid_data(lat_in, lon_in, lat_out, lon_out, input_data, method='linear',threshold=None):

    if len(np.shape(lat_in)) == 1:
        lon_in_2D, lat_in_2D = np.meshgrid(lon_in,lat_in)
        lon_in_1D            = np.reshape(lon_in_2D,-1)
        lat_in_1D            = np.reshape(lat_in_2D,-1)
    elif len(np.shape(lat_in)) == 2:
        lon_in_1D            = np.reshape(lon_in,-1)
        lat_in_1D            = np.reshape(lat_in,-1)
    else:
        print("ERROR: lon_in has ", len(np.shape(lat_in)), "dimensions")

    if len(np.shape(lat_out)) == 1:
        lon_out_2D, lat_out_2D = np.meshgrid(lon_out,lat_out)
    elif len(np.shape(lat_out)) == 2:
        lon_out_2D            = lon_out
        lat_out_2D            = lat_out
    else:
        print("ERROR: lon_out has ", len(np.shape(lat_in)), "dimensions")

    value_tmp = np.reshape(input_data,-1)

    # Check NaN - input array shouldn't have NaN
    if threshold is None:
        mask_values     = ~np.isnan(value_tmp)
    else:
        #print("np.all([~np.isnan(value_tmp),value_tmp>threshold],axis=0) ",np.all([~np.isnan(value_tmp),value_tmp>threshold],axis=0))
        mask_values     = np.all([~np.isnan(value_tmp),value_tmp>threshold],axis=0)

    value     = value_tmp[mask_values]
    # ======= CAUTION =======
    lat_in_1D = lat_in_1D[mask_values]  # here I make nan in values as the standard
    lon_in_1D = lon_in_1D[mask_values]
    #print("shape value = ", np.shape(value))
    #print("shape lat_in_1D = ", np.shape(lat_in_1D))
    #print("shape lon_in_1D = ", np.shape(lon_in_1D))
    # =======================
    #print("value =",value)
    #print("lat_in_1D =",lat_in_1D)
    #print("lon_in_1D =",lon_in_1D)
    Value = griddata((lon_in_1D, lat_in_1D), value, (lon_out_2D, lat_out_2D), method=method)

    return Value

def add_climatology(leap_year,window=31):

    # Check for leap year automatically if not provided
    print('leap_year',leap_year)
    if leap_year =='True':
        ntime = 366
        leap_common = "leap"
    elif leap_year =='False':
        ntime = 365
        leap_common = "common"

    # ========================= Setting up ===========================
    albedo_path    = "/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/"
    lai_path       = "/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/"
    landcover_path = "/g/data/w97/mm3972/data/MODIS/MODIS_landcover/"

    # ========================= Reading MODIS ===========================
    # ------- read MODIS albedo -------
    albedo_band_list = ['BSA_vis', 'BSA_nir', 'WSA_vis', 'WSA_nir']  # match with cable code
    # ssnow%albsoilsn(:,1), BSA_vis: direct beam visible albedo
    # ssnow%albsoilsn(:,2), BSA_nir: direct beam near-infrared albedo
    # ssnow%albsoilsn(:,3): WSA_vis: diffuse visible albedo
    # ssnow%albsoilsn(:,4): WSA_nir: diffuse near-infrared albedo

    albedo_in = np.zeros((4, ntime, 681, 841))

    for i, albedo_band in enumerate(albedo_band_list):
        var_name = 'Albedo_'+albedo_band
        file_albedo = f"{albedo_path}/regrid_2_AWAP_5km_climatology/MCD43A3.061_clim_{albedo_band}_{leap_common}_{window}day_smooth.nc"
        with nc.Dataset(file_albedo, 'r') as f_albedo:
            albedo_in[i, :, :, :] = f_albedo.variables[var_name][:, :, :]
        albedo_in = np.where(albedo_in>1., 0, albedo_in)
        albedo_in = np.where(albedo_in<0,  0, albedo_in)

    # ------- read MODIS LAI -------
    file_lai = f"{lai_path}/regrid_2_AWAP_5km_climatology/MCD15A3H.061_clim_5km_{leap_common}year_{window}day_smooth.nc"
    with nc.Dataset(file_lai, 'r') as f_lai:
        lai_in = f_lai.variables['LAI'][:, :, :]
        lai_in = np.where(lai_in>20, 0.001, lai_in) # set to the value of C%LAI_THRESH in CABLE
        lai_in = np.where(lai_in<0,  0.001, lai_in)

    # ------- read MODIS landcover -------
    file_lc = f"{landcover_path}/var_landcover_only/MCD12Q1.061_landcover_clim_PFT_5km.nc"
    with nc.Dataset(file_lc, 'r') as f_lc:
        lc_in = f_lc.variables['LC'][:, :]

    # ========================== Copy data =============================
    # Open the original NetCDF file
    f_landinfo       = "/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc"
    f_MODIS_landinfo = f'/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_clim_{leap_common}.nc'

    # Open the source NetCDF file in read mode
    f_in  = nc.Dataset(f_landinfo, 'r')
    f_out = nc.Dataset(f_MODIS_landinfo, 'w')

    # Copy global attributes
    for attr in f_in.ncattrs():
        f_out.setncattr(attr, f_in.getncattr(attr))

    # Copy dimensions
    for name, dimension in f_in.dimensions.items():
        f_out.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

    # Copy variables except LAI and Albedo
    for name, variable in f_in.variables.items():
        if name not in ['LAI', 'Albedo']:
            f_out_var = f_out.createVariable(name, variable.datatype, variable.dimensions)
            # Copy variable attributes
            f_out_var.setncatts({attr: variable.getncattr(attr) for attr in variable.ncattrs()})
            # Copy data
            f_out_var[:] = variable[:]

    # Create new dimensions for time and rad
    f_out.createDimension('time1', ntime)
    f_out.createDimension('rad1', 4)

    # Create new LAI variable with dimensions [365 or 366, latitude, longitude]
    lai_var            = f_out.createVariable('LAI', 'f4', ('time1', 'latitude', 'longitude'),fill_value=-9999)
    lai_var.units      = 'm2/m2'
    lai_var.long_name  = 'Leaf area index'
    lai_var[:, :, :]   = lai_in

    # Create new Albedo variable with dimensions [4, 365 or 366, latitude, longitude]
    albedo_var            = f_out.createVariable('Albedo', 'f4', ('rad1', 'time1', 'latitude', 'longitude'),fill_value=-9999)
    albedo_var.units      = 'dimensionless'
    albedo_var.long_name  = 'Surface albedo, 1: BSA_vis, 2: BSA_nir, 3: WSA_vis, 4: WSA_nir"'
    albedo_var[:, :, :, :] = albedo_in

    # copy new landcover
    f_out['iveg'][:]      = lc_in[:,:]
    f_in.close()
    f_out.close()

    return

def add_time_varying(year,window=31):

    # Check for leap year automatically if not provided
    if is_leap(year):
        ntime = 366
    else:
        ntime = 365

    # ========================= Setting up ===========================
    albedo_path    = "/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/"
    lai_path       = "/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/"
    landcover_path = "/g/data/w97/mm3972/data/MODIS/MODIS_landcover/"

    # ========================= Reading MODIS ===========================
    # ------- read MODIS albedo -------
    albedo_band_list = ['BSA_vis', 'BSA_nir', 'WSA_vis', 'WSA_nir']  # match with cable code
    # ssnow%albsoilsn(:,1), BSA_vis: direct beam visible albedo
    # ssnow%albsoilsn(:,2), BSA_nir: direct beam near-infrared albedo
    # ssnow%albsoilsn(:,3): WSA_vis: diffuse visible albedo
    # ssnow%albsoilsn(:,4): WSA_nir: diffuse near-infrared albedo

    albedo_in = np.zeros((4, ntime, 681, 841))

    for i, albedo_band in enumerate(albedo_band_list):
        var_name = 'Albedo_'+albedo_band
        file_albedo = f"{albedo_path}/regrid_2_AWAP_5km_daily/MCD43A3.061_albedo_{albedo_band}_5km_{year}_{window}day_smooth.nc"

        with nc.Dataset(file_albedo, 'r') as f_albedo:
            albedo_in[i, :, :, :] = f_albedo.variables[var_name][:, :, :]
        albedo_in = np.where(albedo_in>1., 0, albedo_in)
        albedo_in = np.where(albedo_in<0,  0, albedo_in)

    # ------- read MODIS LAI -------
    file_lai = f"{lai_path}/regrid_2_AWAP_5km_daily/MCD15A3H.061_500m_aid0001_LAI_regridded_daily_{year}_{window}day_smooth.nc"

    with nc.Dataset(file_lai, 'r') as f_lai:
        lai_in = f_lai.variables['LAI'][:, :, :]
        lai_in = np.where(lai_in>20, 0.001, lai_in) # set to the value of C%LAI_THRESH in CABLE
        lai_in = np.where(lai_in<0,  0.001, lai_in)

    # ------- read MODIS landcover -------
    if year == 2000:
        # for 2000 use the land cover of 2001
        file_lc = f"{landcover_path}var_landcover_only/MCD12Q1.061_landcover_2001_PFT_5km.nc"
    elif year == 2023:
        # for 2023 use the land cover of 2022
        file_lc = f"{landcover_path}var_landcover_only/MCD12Q1.061_landcover_2022_PFT_5km.nc"
    else:
        file_lc = f"{landcover_path}var_landcover_only/MCD12Q1.061_landcover_{year}_PFT_5km.nc"

    with nc.Dataset(file_lc, 'r') as f_lc:
        lc_in = f_lc.variables['LC'][:, :]

    # ========================== Copy data =============================
    # Open the original NetCDF file
    f_landinfo       = "/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc"
    f_MODIS_landinfo = f'/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_time_varying_{year}.nc'

    # Open the source NetCDF file in read mode
    f_in  = nc.Dataset(f_landinfo, 'r')
    f_out = nc.Dataset(f_MODIS_landinfo, 'w')

    # Copy global attributes
    for attr in f_in.ncattrs():
        f_out.setncattr(attr, f_in.getncattr(attr))

    # Copy dimensions
    for name, dimension in f_in.dimensions.items():
        f_out.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

    # Copy variables except LAI and Albedo
    for name, variable in f_in.variables.items():
        if name not in ['LAI', 'Albedo']:
            f_out_var = f_out.createVariable(name, variable.datatype, variable.dimensions)
            # Copy variable attributes
            f_out_var.setncatts({attr: variable.getncattr(attr) for attr in variable.ncattrs()})
            # Copy data
            f_out_var[:] = variable[:]

    # Create new dimensions for time and rad
    f_out.createDimension('time1', ntime)
    f_out.createDimension('rad1', 4)

    # Create new LAI variable with dimensions [365 or 366, latitude, longitude]
    lai_var            = f_out.createVariable('LAI', 'f4', ('time1', 'latitude', 'longitude'),fill_value=-9999)
    lai_var.units      = 'm2/m2'
    lai_var.long_name  = 'Leaf area index'
    lai_var[:, :, :]   = lai_in

    # Create new Albedo variable with dimensions [4, 365 or 366, latitude, longitude]
    albedo_var            = f_out.createVariable('Albedo', 'f4', ('rad1', 'time1', 'latitude', 'longitude'),fill_value=-9999)
    albedo_var.units      = 'dimensionless'
    albedo_var.long_name  = 'Surface albedo, 1: BSA_vis, 2: BSA_nir, 3: WSA_vis, 4: WSA_nir"'
    albedo_var[:, :, :, :] = albedo_in

    # copy new landcover
    f_out['iveg'][:]      = lc_in[:,:]
    f_in.close()
    f_out.close()

    return

def use_SM_ST_after_90year_spinup(leap_year=None, year=None):

    file_in = "/g/data/w97/mm3972/model/cable/runs/Land_drought_rainfall_runs/prepare_SM_ST_for_offline_run/outputs/cable_out_1999.nc"

    with nc.Dataset(file_in, 'r') as f_in:
        SM_in  = f_in.variables['SoilMoist'][:,:, :, :]
        ST_in  = f_in.variables['SoilTemp'][:,:, :, :]
        lat_in = f_in.variables['latitude'][:, :]
        lon_in = f_in.variables['longitude'][:, :]
        SM_in  = np.where(SM_in<0, np.nan, SM_in)
        ST_in  = np.where(ST_in<0, np.nan, ST_in)

    # Open the original NetCDF file
    if leap_year==None:
        file_out = f'/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_time_varying_{year}.nc'
    elif leap_year=='True':
        file_out = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_clim_leap.nc'
    elif leap_year=='False':
        file_out = '/g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_MODIS_LAI_albedo_lc_clim_common.nc'

    # Open the source NetCDF file in read mode
    f_out    = nc.Dataset(file_out, 'r+')
    lat_out  = f_out.variables['latitude'][:]
    lon_out  = f_out.variables['longitude'][:]
    watr_out = f_out.variables['watr'][:]

    for m in np.arange(12):
        for s in np.arange(6):
            SoilMoist_tmp = regrid_data(lat_in, lon_in, lat_out, lon_out, SM_in[m,s,:,:], method='nearest',threshold=None)
            f_out.variables['SoilMoist'][m,s,:,:] = np.where(SoilMoist_tmp > watr_out[s,:,:], SoilMoist_tmp, watr_out[s,:,:]+0.001) # to make sure SM is always larger than watr
            f_out.variables['SoilTemp'][m,s,:,:]  = regrid_data(lat_in, lon_in, lat_out, lon_out, ST_in[m,s,:,:], method='nearest',threshold=None)

    f_out.close()

    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate new landinfo file with MODIS LAI and albedo.')
    parser.add_argument('--leap_year', type=str, choices=['True', 'False'], default=None, help='Is it a leap year? If not provided, it will be calculated.')
    parser.add_argument('--window', type=int, default=None, help='Window?')

    args = parser.parse_args()

    print(args)

    if args.leap_year != None:
    #    add_climatology(args.leap_year, args.window)
       use_SM_ST_after_90year_spinup(args.leap_year)
    else:
       for year in np.arange(2000,2024,1):
           # add_time_varying(year, args.window) # args.window
           use_SM_ST_after_90year_spinup(year=year)
