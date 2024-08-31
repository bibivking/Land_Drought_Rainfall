#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=1:00:00
#PBS -l mem=500GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97

module use /g/data/hh5/public/modules
module load conda/analysis3-23.10

# Call the Python script with the file path and window size
file_dir="/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/"
albedo_band="WSA_nir"
window=31

# # Step 9: cut to the year
# cd /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/
# echo "Step 9: cut to the year "
#
# echo "2000, starts from Feb 2000"
# cdo seldate,2000-02-24T00:00:00,2000-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Feb2000-Jan2001_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Feb-Dec2000_albedo_regridded_daily_gapfill.nc
#
# echo "2001"
# cdo seldate,2001-01-01T00:00:00,2001-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2000-Jan2002_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2001_albedo_regridded_daily_gapfill.nc
#
# echo "2002"
# cdo seldate,2002-01-01T00:00:00,2002-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2001-Jan2003_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2002_albedo_regridded_daily_gapfill.nc
#
# echo "2003"
# cdo seldate,2003-01-01T00:00:00,2003-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2002-Jan2004_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2003_albedo_regridded_daily_gapfill.nc
#
# echo "2004"
# cdo seldate,2004-01-01T00:00:00,2004-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2003-Jan2005_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2004_albedo_regridded_daily_gapfill.nc
#
# echo "2005"
# cdo seldate,2005-01-01T00:00:00,2005-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2004-Jan2006_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2005_albedo_regridded_daily_gapfill.nc
#
# echo "2006"
# cdo seldate,2006-01-01T00:00:00,2006-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2005-Jan2007_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2006_albedo_regridded_daily_gapfill.nc
#
# echo "2007"
# cdo seldate,2007-01-01T00:00:00,2007-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2006-Jan2008_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2007_albedo_regridded_daily_gapfill.nc
#
# echo "2008"
# cdo seldate,2008-01-01T00:00:00,2008-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2007-Jan2009_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2008_albedo_regridded_daily_gapfill.nc
#
# echo "2009"
# cdo seldate,2009-01-01T00:00:00,2009-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2008-Jan2010_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2009_albedo_regridded_daily_gapfill.nc
#
# echo "2010"
# cdo seldate,2010-01-01T00:00:00,2010-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2009-Jan2011_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2010_albedo_regridded_daily_gapfill.nc
#
# echo "2011"
# cdo seldate,2011-01-01T00:00:00,2011-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2010-Jan2012_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2011_albedo_regridded_daily_gapfill.nc
#
# echo "2012"
# cdo seldate,2012-01-01T00:00:00,2012-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2011-Jan2013_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2012_albedo_regridded_daily_gapfill.nc
#
# echo "2013"
# cdo seldate,2013-01-01T00:00:00,2013-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2012-Jan2014_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2013_albedo_regridded_daily_gapfill.nc
#
# echo "2014"
# cdo seldate,2014-01-01T00:00:00,2014-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2013-Jan2015_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2014_albedo_regridded_daily_gapfill.nc
#
# echo "2015"
# cdo seldate,2015-01-01T00:00:00,2015-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2014-Jan2016_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2015_albedo_regridded_daily_gapfill.nc
#
# echo "2016"
# cdo seldate,2016-01-01T00:00:00,2016-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2015-Jan2017_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2016_albedo_regridded_daily_gapfill.nc
#
# echo "2017"
# cdo seldate,2017-01-01T00:00:00,2017-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2016-Jan2018_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2017_albedo_regridded_daily_gapfill.nc
#
# echo "2018"
# cdo seldate,2018-01-01T00:00:00,2018-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2017-Jan2019_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2018_albedo_regridded_daily_gapfill.nc
#
# echo "2019"
# cdo seldate,2019-01-01T00:00:00,2019-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2018-Jan2020_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2019_albedo_regridded_daily_gapfill.nc
#
# echo "2020"
# cdo seldate,2020-01-01T00:00:00,2020-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2019-Jan2021_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2020_albedo_regridded_daily_gapfill.nc
#
# echo "2021"
# cdo seldate,2021-01-01T00:00:00,2021-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2020-Jan2022_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2021_albedo_regridded_daily_gapfill.nc
#
# echo "2022"
# cdo seldate,2022-01-01T00:00:00,2022-12-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2021-Jan2023_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2022_albedo_regridded_daily_gapfill.nc
#
# echo "2023, with Jan 2024"
# cdo seldate,2023-01-01T00:00:00,2024-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_Dec2022-Jan2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_2023-Jan2024_albedo_regridded_daily_gapfill.nc
#
# # Step 10: fill the date without nearby values by using climatology
# file_dir="/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/"
# cd /g/data/w97/mm3972/scripts/Drought/Land_drought_rainfall/process_for_CABLE/process_albedo
#
# # echo "Step 10: fill the date without nearby values by using climatology"
# echo "Feb-Dec2000"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_Feb-Dec2000_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2001"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2001_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2002"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2002_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2003"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2003_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2004"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2004_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2005"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2005_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2006"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2006_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2007"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2007_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2008"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2008_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2009"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2009_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2010"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2010_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2011"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2011_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2012"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2012_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2013"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2013_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2014"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2014_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2015"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2015_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2016"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2016_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2017"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2017_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2018"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2018_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2019"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2019_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2020"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2020_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2021"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2021_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2022"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2022_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"
# echo "2023-Jan2024"
# python process_MODIS_albedo_time_varying_step10_gap_fill_by_climatology.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_2023-Jan2024_albedo_regridded_daily_gapfill.nc" --albedo_band "$albedo_band"

# # Step 11: merge the gapfill data
# cd /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/
# echo "Step 11: merge the gapfill data"
# ls MCD43A3.061_500m_aid0001_${albedo_band}_*_albedo_regridded_daily_gapfill.nc
# cdo mergetime MCD43A3.061_500m_aid0001_${albedo_band}_*_albedo_regridded_daily_gapfill.nc MCD43A3.061_500m_aid0001_${albedo_band}_Feb2000-Jan2024_albedo_regridded_daily_gapfill.nc

# # Step 12: 31-day running mean
# cd /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/
# echo "Step 12: 31-day running mean"
# cdo runmean,31 MCD43A3.061_500m_aid0001_${albedo_band}_Feb2000-Jan2024_albedo_regridded_daily_gapfill.nc MCD43A3.061_500m_aid0001_${albedo_band}_Feb2000-Jan2024_albedo_regridded_daily_gapfill_31day_smooth.nc

# Step 13: fill Jan-Feb 2000 with climatology, save to yearly nc file
cd /g/data/w97/mm3972/scripts/Drought/Land_drought_rainfall/process_for_CABLE/process_albedo
echo "Step 13: save to yearly nc file"
python process_MODIS_albedo_time_varying_step13_save_annual_nc.py --albedo_band "$albedo_band" --window $window
