import argparse
import netCDF4 as nc
import numpy as np
import pandas as pd


def day_in_year(year):
    """Check if a year is a leap year."""
    if (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0):
        return 366
    else:
        return 365

def extend_LAI_to_2000_by_climology(window_vary=15,window_clim=61):

    """
    Extend LAI data to the year 2000 by using climatology data.
    using window_vary day smoothing at the connection point (+- window_vary)
    """

    file_in          = f"/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily/MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2002-2024_{window_vary}day_smooth.nc"
    file_clim_common = f"/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_climatology/MCD15A3H.061_clim_5km_commonyear_{window_clim}day_smooth.nc"
    file_clim_leap   = f"/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_climatology/MCD15A3H.061_clim_5km_leapyear_{window_clim}day_smooth.nc"
    file_out         = f"/g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily/MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2000-2023_{window_vary}day_smooth.nc"

    # Read climatology
    with nc.Dataset(file_clim_common, 'r') as f_clim_common:
        lai_clim_common = f_clim_common.variables['LAI'][:, :, :].data

        #print('np.unique(f_clim_common.variables["LAI"][:, :, :].data)', np.unique(f_clim_common.variables['LAI'][:, :, :].data),
        #      'np.unique(f_clim_common.variables["LAI"][:, :, :])',np.unique(f_clim_common.variables['LAI'][:, :, :]))
        # summary: with .data or not will return the same values but the data types are different

    with nc.Dataset(file_clim_leap, 'r') as f_clim_leap:
        lai_clim_leap = f_clim_leap.variables['LAI'][:, :, :].data

    # Read smoothed daily LAI
    with nc.Dataset(file_in, 'r') as f_in:
        lai_in = f_in.variables['Lai_500m'][:, :, :].data
        lai_in = np.where(lai_in>20., 0.001, lai_in) # set to the value of C%LAI_THRESH in CABLE
        lat_in = f_in.variables['latitude'][:].data
        lon_in = f_in.variables['longitude'][:].data

    # test difference between with and without data

    nlat             = len(lat_in)
    nlon             = len(lon_in)
    ntime            = sum(day_in_year(year) for year in np.arange(2000, 2024))
    print('ntime', ntime)

    # Reading lai
    lai_out                  = np.zeros((ntime, nlat, nlon))
    # 2000
    lai_out[:366,:,:]        = lai_clim_leap
    # 2001
    lai_out[366:366+365,:,:] = lai_clim_common
    # 1st Jan 2002 ~ 3th July 2002
    lai_out[366+365:915,:,:] = lai_clim_common[:184,:,:]
    # 4th July 2002 ~ 31st Dec 2023
    lai_out[915:,:,:]        = lai_in[:ntime-915,:,:]

    # Running 11 day smoothing during June 30th ~ 9th July

    for i in np.arange(nlat):
        for j in np.arange(nlon):
            lai_out[915-window_vary+1:915+window_vary,i,j] = pd.Series(lai_out[915-window_vary+1:915+window_vary, i, j]).rolling(window=window_vary, min_periods=1, center=True).mean().values
            print("i,j,np.unique(lai_out[915-window_vary+1:915+window_vary,i,j])",i,j,np.unique(lai_out[915-window_vary+1:915+window_vary,i,j]))

    # Create output file
    f_out                    = nc.Dataset(file_out, 'w')

    # Copy global attributes
    for attr in f_in.ncattrs():
        f_out.setncattr(attr, f_in.getncattr(attr))

    # Create dimensions
    f_out.createDimension('time', ntime)
    f_out.createDimension('latitude', nlat)
    f_out.createDimension('longitude', nlon)

    Time               = f_out.createVariable('time', 'i4', ('time'))
    Time.standard_name = "time"
    Time.units         = "days since 2000-01-01 00:00:00"
    Time.calendar      = "gregorian"
    Time.axis          = "T"
    Time[:]            = np.arange(0,ntime,1)

    Longitude               = f_out.createVariable('longitude', 'f4', ('longitude'))
    Longitude.standard_name = "longitude"
    Longitude.long_name     = "longitude"
    Longitude.units         = "degrees_East"
    Longitude.axis          = "X"
    Longitude[:]            = lon_in[:]

    Latitude                = f_out.createVariable('latitude', 'f4', ('latitude'))
    Latitude.standard_name = "latitude"
    Latitude.long_name     = "latitude"
    Latitude.units         = "degrees_North"
    Latitude.axis          = "Y"
    Latitude[:]            = lat_in[:]

    var_out                = f_out.createVariable('LAI', 'f4', ('time', 'latitude', 'longitude'),fill_value=-9999.)
    var_out.long_name      = "MCD15A3H MODIS/Terra Gridded 500M Leaf Area Index LAI (4-day composite)"
    var_out.standard_name  = "LAI"
    var_out.units          = "m^2/m^2"
    var_out[:]             = lai_out[:]

    # Lai_500m               = f_out.createVariable('Lai_500m', 'i2', ('time', 'latitude', 'longitude'), fill_value=255)
    # Lai_500m.long_name     = "MCD15A3H MODIS/Terra Gridded 500M Leaf Area Index LAI (4-day composite)"
    # Lai_500m.units         = "m^2/m^2"
    # Lai_500m.add_offset    = 0.
    # Lai_500m.scale_factor  = 0.1
    # # Lai_500m._FillValue = 255
    # Lai_500m.missing_value = 255
    # Lai_500m.scale_factor_err = 0.
    # Lai_500m.add_offset_err= 0.
    # Lai_500m.calibrated_nt = 21
    # Lai_500m.MOD15A2_FILLVALUE_DOC = "MOD15A2 FILL VALUE LEGEND\n255 = _Fillvalue, assigned when:\n    * the MOD09GA suf. reflectance for channel VIS, NIR was assigned its _Fillvalue, or\n    * land cover pixel itself was assigned _Fillvalus 255 or 254.\n254 = land cover assigned as perennial salt or inland fresh water.\n253 = land cover assigned as barren, sparse vegetation (rock, tundra, desert.)\n252 = land cover assigned as perennial snow, ice.\n251 = land cover assigned as \"permanent\" wetlands/inundated marshlands.\n250 = land cover assigned as urban/built-up.\n249 = land cover assigned as \"unclassified\" or not able to determine.\n"
    # Lai_500m[:]            = lai_out[:]

    f_out.close()

    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Extending time-varying LAI to 2000')
    parser.add_argument('--window_vary', type=int, default=0, help='Window size for smoothing the time-varying file')
    parser.add_argument('--window_clim', type=int, default=0, help='Window size for smoothing the climatology')
    args = parser.parse_args()

    extend_LAI_to_2000_by_climology(args.window_vary, args.window_clim)
