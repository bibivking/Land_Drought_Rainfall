#!/bin/bash
#PBS -m ae
#PBS -P w97
#PBS -q hugemem
#PBS -l walltime=3:00:00
#PBS -l mem=1000GB
#PBS -l ncpus=48
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97

module use /g/data/hh5/public/modules
module load conda/analysis3

cd /g/data/w97/mm3972/scripts/Drought/Post_drought_rainfall/process_for_CABLE/process_albedo

cp process_MODIS_albedo_time_varying_step8_gap_fill.py process_MODIS_albedo_time_varying_step8_gap_fill_WSA_vis_2000.py
if [ 2000 -eq 2000 ]; then
    python process_MODIS_albedo_time_varying_step8_gap_fill_WSA_vis_2000.py "/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/MCD43A3.061_500m_aid0001_WSA_vis_Feb2000-Jan2001_albedo_regridded_daily.nc" --albedo_band "WSA_vis" --window "30"
else
    python process_MODIS_albedo_time_varying_step8_gap_fill_WSA_vis_2000.py "/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/MCD43A3.061_500m_aid0001_WSA_vis_Dec1999-Jan2001_albedo_regridded_daily.nc" --albedo_band "WSA_vis" --window "30"
fi
rm process_MODIS_albedo_time_varying_step8_gap_fill_WSA_vis_2000.py

