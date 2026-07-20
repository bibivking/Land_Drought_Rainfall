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
albedo_band="BSA_nir"

# Step 5: to one file
cd /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS
echo "Step 5: to one file"
cdo mergetime ./regrid_2_AWAP_5km/MCD43A3.061_500m_aid0001_${albedo_band}_????_albedo_regridded.nc ./regrid_2_AWAP_5km/MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded.nc

# Step 6: to daily
cd /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS
echo "Step 6: to daily"
cdo inttime,2000-02-24,00:00:00,1day ./regrid_2_AWAP_5km/MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded.nc ./regrid_2_AWAP_5km_daily/MCD43A3.061_500m_aid0001_${albedo_band}_2000-2024_albedo_regridded_daily.nc
