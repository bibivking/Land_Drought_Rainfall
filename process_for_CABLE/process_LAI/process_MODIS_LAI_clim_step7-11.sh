#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=4:00:00
#PBS -l mem=300GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97

module use /g/data/hh5/public/modules
module load conda/analysis3

# # Step 7: cut the nc file to 2003-1-1 ~ 2023-12-31
# cd /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily
# echo "Step 7: cut the nc file to 2003-1-1 ~ 2023-12-31"
# cdo seldate,2003-01-01T00:00:00,2023-12-31T23:59:59 MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2002-2024.nc MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2003-2023.nc

# # Step 8: remove 29 Feb
# cd /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily
# echo "Step 8: remove 29 Feb"
# cdo del29feb MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2003-2023.nc MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2003-2023_no29Feb.nc

# # Step 9: to climatology
# cd /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily
# echo "Step 9: to climatology"
# cdo ydaymean MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2003-2023_no29Feb.nc ../regrid_2_AWAP_5km_climatology/MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2003-2023_climatology_common.nc

# # Step 10: save climatology to nc file used in CABLE
# cd /g/data/w97/mm3972/scripts/Drought/Post_drought_rainfall/process_for_CABLE/process_LAI
# echo "Step 10: save climatology to nc file used in CABLE"
# python process_MODIS_LAI_clim_step10_save_nc.py

# Step 11: calculate xxx-day running mean to smooth the climatology LAI
window=31
cd /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_climatology
echo "Step 11: calculate $window-day running mean to smooth the climatology LAI"

echo "Calculate common year"
cdo settaxis,2002-01-01,00:00:00,1day MCD15A3H.061_2003-2023_climatology_5km_common_year.nc MCD15A3H.061_2003-2023_climatology_5km_common_year_additional.nc
cdo cat MCD15A3H.061_2003-2023_climatology_5km_leap_year.nc MCD15A3H.061_2003-2023_climatology_5km_common_year.nc MCD15A3H.061_2003-2023_climatology_5km_common_year_additional.nc MCD15A3H.061_2003-2023_climatology_5km_common_year_extent.nc
cdo runmean,$window MCD15A3H.061_2003-2023_climatology_5km_common_year_extent.nc "MCD15A3H.061_2003-2023_climatology_5km_common_year_${window}smooth.nc"
cdo seldate,2001-01-01T00:00:00,2001-12-31T23:59:59 "MCD15A3H.061_2003-2023_climatology_5km_common_year_${window}smooth.nc" "MCD15A3H.061_clim_5km_commonyear_${window}day_smooth.nc"
rm MCD15A3H.061_2003-2023_climatology_5km_common_year_additional.nc MCD15A3H.061_2003-2023_climatology_5km_common_year_extent.nc "MCD15A3H.061_2003-2023_climatology_5km_common_year_${window}smooth.nc"

echo "Calculate leap year"
cdo settaxis,1999-01-01,00:00:00,1day MCD15A3H.061_2003-2023_climatology_5km_common_year.nc MCD15A3H.061_2003-2023_climatology_5km_common_year_additional.nc
cdo cat MCD15A3H.061_2003-2023_climatology_5km_common_year_additional.nc MCD15A3H.061_2003-2023_climatology_5km_leap_year.nc MCD15A3H.061_2003-2023_climatology_5km_common_year.nc MCD15A3H.061_2003-2023_climatology_5km_leap_year_extent.nc
cdo runmean,$window MCD15A3H.061_2003-2023_climatology_5km_leap_year_extent.nc "MCD15A3H.061_2003-2023_climatology_5km_leap_year_${window}smooth.nc"
cdo seldate,2000-01-01T00:00:00,2000-12-31T23:59:59 "MCD15A3H.061_2003-2023_climatology_5km_leap_year_${window}smooth.nc" "MCD15A3H.061_clim_5km_leapyear_${window}day_smooth.nc"
rm MCD15A3H.061_2003-2023_climatology_5km_common_year_additional.nc MCD15A3H.061_2003-2023_climatology_5km_leap_year_extent.nc "MCD15A3H.061_2003-2023_climatology_5km_leap_year_${window}smooth.nc"