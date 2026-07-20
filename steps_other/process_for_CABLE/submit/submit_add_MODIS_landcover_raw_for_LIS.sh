#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=6:00:00
#PBS -l mem=500GB
#PBS -l ncpus=104
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97+gdata/ua8

module use /g/data/hh5/public/modules
module load conda/analysis3-23.10

cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE
python add_MODIS_landcover_raw_for_LIS.py
