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

# to Gregorian calendar since cdo inttime cannot process Julian calendar
# ncatted -O -a calendar,time,m,c,gregorian MCD15A3H.061_500m_aid0001_2002-2005.nc

# out put AWAP resolution file's lat lon information
# cdo griddes /g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc > target_grid.txt

# extract LAI from files
# cdo selname,Lai_500m /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/MCD15A3H.061_500m_aid0001_2002-2005.nc /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/var_LAI_only/MCD15A3H.061_500m_aid0001_2002-2005_LAI.nc

# resampling
# cdo remapbil,target_grid.txt /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/var_LAI_only/MCD15A3H.061_500m_aid0001_2002-2005_LAI.nc /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km/MCD15A3H.061_500m_aid0001_2002-2005_LAI_regridded.nc

# to daily
cdo inttime,2002-07-04,00:00:00,1day /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km/MCD15A3H.061_500m_aid0001_2002-2005_LAI_regridded.nc /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily/MCD15A3H.061_500m_aid0001_2002-2005_LAI_regridded_daily.nc

# nccompress -r -o -np 8 /g/data/w35/mm3972/model/cable/runs/*/*/outputs

