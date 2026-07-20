#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalbw
#PBS -l walltime=0:20:00
#PBS -l mem=80GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97

module use /g/data/hh5/public/modules
module load conda/analysis3-23.10

# Set the path to search
albedo_band="WSA_nir"

# Step 7: cut to Dec to next year Jan, for gap filling
cd /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/
echo "Step 7: cut from Dec year before to Jan year after"

echo "2000"
cdo seldate,2000-02-24T00:00:00,2001-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Feb2000-Jan2001_albedo_regridded_daily.nc

# echo "2001"
# cdo seldate,2000-12-01T00:00:00,2002-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2000-Jan2002_albedo_regridded_daily.nc

# echo "2002"
# cdo seldate,2001-12-01T00:00:00,2003-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2001-Jan2003_albedo_regridded_daily.nc

# echo "2003"
# cdo seldate,2002-12-01T00:00:00,2004-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2002-Jan2004_albedo_regridded_daily.nc

# echo "2004"
# cdo seldate,2003-12-01T00:00:00,2005-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2003-Jan2005_albedo_regridded_daily.nc

# echo "2005"
# cdo seldate,2004-12-01T00:00:00,2006-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2004-Jan2006_albedo_regridded_daily.nc

# echo "2006"
# cdo seldate,2005-12-01T00:00:00,2007-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2005-Jan2007_albedo_regridded_daily.nc

# echo "2007"
# cdo seldate,2006-12-01T00:00:00,2008-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2006-Jan2008_albedo_regridded_daily.nc

# echo "2008"
# cdo seldate,2007-12-01T00:00:00,2009-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2007-Jan2009_albedo_regridded_daily.nc

# echo "2009"
# cdo seldate,2008-12-01T00:00:00,2010-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2008-Jan2010_albedo_regridded_daily.nc

# echo "2010"
# cdo seldate,2009-12-01T00:00:00,2011-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2009-Jan2011_albedo_regridded_daily.nc

# echo "2011"
# cdo seldate,2010-12-01T00:00:00,2012-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2010-Jan2012_albedo_regridded_daily.nc

# echo "2012"
# cdo seldate,2011-12-01T00:00:00,2013-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2011-Jan2013_albedo_regridded_daily.nc

# echo "2013"
# cdo seldate,2012-12-01T00:00:00,2014-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2012-Jan2014_albedo_regridded_daily.nc

# echo "2014"
# cdo seldate,2013-12-01T00:00:00,2015-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2013-Jan2015_albedo_regridded_daily.nc

# echo "2015"
# cdo seldate,2014-12-01T00:00:00,2016-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2014-Jan2016_albedo_regridded_daily.nc

# echo "2016"
# cdo seldate,2015-12-01T00:00:00,2017-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2015-Jan2017_albedo_regridded_daily.nc

# echo "2017"
# cdo seldate,2016-12-01T00:00:00,2018-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2016-Jan2018_albedo_regridded_daily.nc

# echo "2018"
# cdo seldate,2017-12-01T00:00:00,2019-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2017-Jan2019_albedo_regridded_daily.nc

# echo "2019"
# cdo seldate,2018-12-01T00:00:00,2020-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2018-Jan2020_albedo_regridded_daily.nc

# echo "2020"
# cdo seldate,2019-12-01T00:00:00,2021-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2019-Jan2021_albedo_regridded_daily.nc

# echo "2021"
# cdo seldate,2020-12-01T00:00:00,2022-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2020-Jan2022_albedo_regridded_daily.nc

# echo "2022"
# cdo seldate,2021-12-01T00:00:00,2023-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2021-Jan2023_albedo_regridded_daily.nc

# echo "2023"
# cdo seldate,2022-12-01T00:00:00,2024-01-31T23:59:59 MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc MCD43A3.061_500m_aid0001_${albedo_band}_Dec2022-Jan2024_albedo_regridded_daily.nc

# Step 7: compress ncfile
# nccompress -r -o -np 8 /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc
