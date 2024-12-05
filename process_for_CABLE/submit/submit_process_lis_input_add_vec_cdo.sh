#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=6:00:00
#PBS -l mem=500GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97+gdata/ua8

module use /g/data/hh5/public/modules
module load conda/analysis3-23.10

cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/regrid_lis_input
# cdo griddes /scratch/w97/mm3972/model/NUWRF/Drought_breaking_rainfall/test_runs/LIS_offline/lis_input.d01.nc > lambert_grid.txt
cdo genbil,lambert_grid.txt /scratch/w97/mm3972/model/NUWRF/Drought_breaking_rainfall/test_runs/LIS_offline/lis_input.d01.nc weights.nc
cdo remap,lambert_grid.txt,weights.nc /g/data/w97/mm3972/scripts/wrf_scripts/make_LIS_landinfo/nc_file/Openlandmap_soilcomposition_CORDEX_180E_depth_varying.nc /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/nc_files/regrid_lis_input/Openlandmap_soilcomposition_CORDEX_180E_depth_varying_lis_domain1.nc
# cdo remapbil,/scratch/w97/mm3972/model/NUWRF/Drought_breaking_rainfall/test_runs/LIS_offline/lis_input.d01.nc /g/data/w97/mm3972/scripts/wrf_scripts/make_LIS_landinfo/nc_file/Openlandmap_soilcomposition_CORDEX_180E_depth_varying.nc Openlandmap_soilcomposition_CORDEX_180E_depth_varying_lis_domain1.nc