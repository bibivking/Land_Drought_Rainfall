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

n=2
case_name='test_runs'

cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE
python process_lis_input_add_vec.py --case_name ${case_name} --n ${n}


n=3
case_name='test_runs'

cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE
python process_lis_input_add_vec.py --case_name ${case_name} --n ${n}
