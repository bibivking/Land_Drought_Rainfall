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
module load conda/analysis3

# Step 7: using climatology to derive data back to 2000, smooth the connection between clim and time-varying LAI, and save 2000-2023
cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/process_LAI/LAI_time_varying_step7-9_method_2_smooth_anomaly
echo "Step 7: using climatology to derive data back to 2000"
window_vary=31
window_clim=31
python process_MODIS_LAI_time_varying_step7_extend_to_2000_SG.py --window_vary "$window_vary" --window_clim "$window_clim"

# Step 8: remove high frequent varibility
echo "Step 8: remove high frequent varibility"
python process_MODIS_LAI_time_varying_step8_remove_high_frequent_varibility_SG.py

# Step 9: to one file
echo "Step 9: to one file"
cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/process_LAI/LAI_time_varying_step7-9_method_2_smooth_anomaly
python process_MODIS_LAI_time_varying_step9_save_yearly_nc_SG.py
