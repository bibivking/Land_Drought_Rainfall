#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=4:00:00
#PBS -l mem=200GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+gdata/zv2

module use /g/data/hh5/public/modules
module load conda/analysis3

cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/nc_files
cdo seldate,2000-01-01,2023-12-31 spi_pearson_90_reorder_drought_periods.nc spi_pearson_90_reorder_drought_periods_2000-2023.nc