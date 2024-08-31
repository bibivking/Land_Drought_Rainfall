#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalbw
#PBS -l walltime=4:00:00
#PBS -l mem=100GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97


module use /g/data/hh5/public/modules
module load conda/analysis3

# Step 1: to Gregorian calendar since cdo inttime cannot process Julian calendar
# ncatted -O -a calendar,time,m,c,gregorian MCD15A3H.061_500m_aid0001_2002-2005.nc

# Step 2: out put AWAP resolution file's lat lon information
# cdo griddes /g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc > target_grid.txt

# Step 3: extract LAI from files
# cdo selname,Lai_500m /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/MCD15A3H.061_500m_aid0001_2006-2010.nc /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/var_LAI_only/MCD15A3H.061_500m_aid0001_2006-2010_LAI.nc

# Step 4: resampling
# cdo remapbil,target_grid.txt /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/var_LAI_only/MCD15A3H.061_500m_aid0001_2006-2010_LAI.nc /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km/MCD15A3H.061_500m_aid0001_2006-2010_LAI_regridded.nc

