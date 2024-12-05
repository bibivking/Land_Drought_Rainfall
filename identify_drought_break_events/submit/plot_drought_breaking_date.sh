#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=4:00:00
#PBS -l mem=500GB
#PBS -l ncpus=104
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+gdata/zv2+gdata/rt52

module use /g/data/hh5/public/modules
module load conda/analysis3

cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events/
python plot_drought_breaking_date.py
