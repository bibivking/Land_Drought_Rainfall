#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=4:00:00
#PBS -l mem=400GB
#PBS -l ncpus=8
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97

module use /g/data/hh5/public/modules
module load conda/analysis3

# Step 5: to one file
cd /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km
echo 'Step 5: to one file'
cdo mergetime MCD15A3H.061_500m_aid0001_*-*_LAI_*regridded.nc MCD15A3H.061_500m_aid0001_LAI_regridded_2002-2024.nc

# Step 6: to daily
cd /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km
echo "Step 6: to daily"
cdo inttime,2002-07-04,00:00:00,1day /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km/MCD15A3H.061_500m_aid0001_LAI_regridded_2002-2024.nc /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily/MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2002-2024.nc
