#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=1:00:00
#PBS -l mem=100GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97

module use /g/data/hh5/public/modules
module load conda/analysis3-23.10

# Call the Python script with the file path and window size
albedo_band="WSA_nir"

# # Step 7: cut to 2003-2023
# cd /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/
# echo "Step 7: cut to 2003-2023"
# cdo seldate,2003-01-01T00:00:00,2023-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2003-2023_albedo_regridded_daily.nc

# # Step 8: remove 29 Feb
# cd /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/
# echo "Step 8: remove 29 Feb"
# cdo del29feb MCD43A3.061_500m_aid0001_${albedo_band}_2003-2023_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2003-2023_albedo_regridded_daily_no29Feb.nc

# # Step 9: to climatology
# cd /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/
# echo "Step 9: to climatology"
# cdo ydaymean MCD43A3.061_500m_aid0001_${albedo_band}_2003-2023_albedo_regridded_daily_no29Feb.nc ../regrid_2_AWAP_5km_climatology/MCD43A3.061_500m_aid0001_${albedo_band}_2003-2023_climatology_common.nc

# Step 10: save climatology to nc file used in CABLE
window=30
cd /g/data/w97/mm3972/scripts/Drought/Land_drought_rainfall/process_for_CABLE/process_albedo
echo "Step 10: save climatology to nc file used in CABLE"
cp process_MODIS_albedo_clim_step10_save_nc.py process_MODIS_albedo_clim_step10_save_nc_${albedo_band}.py
python process_MODIS_albedo_clim_step10_save_nc_${albedo_band}.py --albedo_band ${albedo_band} --window "${window}"
rm process_MODIS_albedo_clim_step10_save_nc_${albedo_band}.py

# Step 11: xxx day smooth of albedo
window=31
echo "Calculate common year"
cd /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_climatology
cdo settaxis,2002-01-01,00:00:00,1day MCD43A3.061_2003-2023_${albedo_band}_climatology_common.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_common_additional.nc
cdo cat MCD43A3.061_2003-2023_${albedo_band}_climatology_leap.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_common.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_common_additional.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_common_extent.nc
cdo runmean,${window} MCD43A3.061_2003-2023_${albedo_band}_climatology_common_extent.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_common_${window}day_smooth.nc
cdo seldate,2001-01-01T00:00:00,2001-12-31T23:59:59 MCD43A3.061_2003-2023_${albedo_band}_climatology_common_${window}day_smooth.nc MCD43A3.061_clim_${albedo_band}_common_${window}day_smooth.nc
rm MCD43A3.061_2003-2023_${albedo_band}_climatology_common_additional.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_common_extent.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_common_${window}day_smooth.nc

echo "Calculate leap year"
cd /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_climatology
cdo settaxis,1999-01-01,00:00:00,1day MCD43A3.061_2003-2023_${albedo_band}_climatology_common.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_common_additional.nc
cdo cat MCD43A3.061_2003-2023_${albedo_band}_climatology_common_additional.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_leap.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_common.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_leap_extent.nc
cdo runmean,${window} MCD43A3.061_2003-2023_${albedo_band}_climatology_leap_extent.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_leap_${window}day_smooth.nc
cdo seldate,2000-01-01T00:00:00,2000-12-31T23:59:59 MCD43A3.061_2003-2023_${albedo_band}_climatology_leap_${window}day_smooth.nc MCD43A3.061_clim_${albedo_band}_leap_${window}day_smooth.nc
rm MCD43A3.061_2003-2023_${albedo_band}_climatology_common_additional.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_leap_extent.nc MCD43A3.061_2003-2023_${albedo_band}_climatology_leap_${window}day_smooth.nc
