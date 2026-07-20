#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=1:00:00
#PBS -l mem=100GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97

module use /g/data/hh5/public/modules
module load conda/analysis3

n=3 
is_clim='True'
case_name='test_runs'
start_date='20051130'
end_date='20070201'

echo $n $is_clim $case_name $start_date $end_date
cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE
python -X faulthandler process_lis_input_add_MODIS_albedo_lai.py --n $n --is_clim $is_clim --case_name $case_name --start_date $start_date --end_date $end_date
