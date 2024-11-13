#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=4:00:00
#PBS -l mem=300GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97

module use /g/data/hh5/public/modules
module load conda/analysis3-23.10

window=31

cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE

# leap_year='True'
# python add_MODIS_LAI_SG_albedo_to_CABLE_land_info.py --leap_year $leap_year --window $window

# leap_year='False'
# python add_MODIS_LAI_SG_albedo_to_CABLE_land_info.py --leap_year $leap_year --window $window

python add_MODIS_LAI_SG_albedo_to_CABLE_land_info.py --window $window
