#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=4:00:00
#PBS -l mem=500GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+gdata/zv2

module use /g/data/hh5/public/modules
module load conda/analysis3

cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/
# Step 1
#cdo mergetime /g/data/zv2/agcd/v1-0-2/precip/total/r005/01day/agcd_v1_precip_total_r005_daily_19{50..99}.nc /g/data/zv2/agcd/v1-0-2/precip/total/r005/01day/agcd_v1_precip_total_r005_daily_20{00..23}.nc ./nc_files/agcd_v1_precip_total_r005_daily_1950_2023.nc

# Step 2
#ncpdq -a lat,lon,time ./nc_files/agcd_v1_precip_total_r005_daily_1950_2023.nc ./nc_files/agcd_v1_precip_total_r005_daily_1950_2023_reordered.nc

## Step 3
#source /home/561/mm3972/myenv/bin/activate
#process_climate_indices --index spi  --periodicity daily --netcdf_precip ./nc_files/agcd_v1_precip_total_r005_daily_1950_2023_reordered.nc --var_name_precip precip --output_file_base ./nc_files --scales 30 90 --calibration_start_year 1950 --calibration_end_year 1999 --multiprocessing all

# Step 4
ncpdq -a time,lat,lon ./nc_files/spi_gamma_30.nc ./nc_files/spi_gamma_30_reorder.nc
ncpdq -a time,lat,lon ./nc_files/spi_gamma_90.nc ./nc_files/spi_gamma_90_reorder.nc
ncpdq -a time,lat,lon ./nc_files/spi_pearson_30.nc ./nc_files/spi_pearson_30_reorder.nc
ncpdq -a time,lat,lon ./nc_files/spi_pearson_90.nc ./nc_files/spi_pearson_90_reorder.nc
