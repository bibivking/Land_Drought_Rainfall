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

cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/process_LAI/remove_high_frequent_varibility
python process_MODIS_LAI_time_varying_step9_remove_high_frequent_varibility.py


# # Step 7: smoothing LAI:  Option 1 running smooth
  # # Note that: I checked the daily data, too much high-frequency varibility, need a >10-day running smooth,
  # # so I tried 11-day. However, 11-day smooth has too much high-frequency varibility. The varibility
  # # of 11-day-mean time-varying and 91-day-mean climatology LAI are too different in the output of step 8.
  # # Thus, I increase time-varying's smoothing window to 15 day and reduce the climatology's window to 61 day.
  # # I then tried 15 day mean in time-varying vs 61 day mean in climatology. I realized if they use different
  # # smoothing windows, when I compare simulations using them, the differences between simulations may come from
  # # the high frequent variablity rather than the longer term trend (drought vs wet year). So using 31 day mean
  # # for both climatology and time-varying data.

# window=31
# cd /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily
# echo "Step 7: $window day running mean"
# cdo runmean,$window MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2002-2024.nc "MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2002-2024_${window}day_smooth.nc"

# # Step 7: smoothing LAI:  Option 2 cubic spline smooth
# window=31
# cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/process_for_CABLE/process_LAI/remove_high_frequent_varibility
# echo "Step 7: $window day running smooth"
# python run_cubic_spline_smooth.py

# # Step 8: using climatology to derive data back to 2000, smooth the connection between clim and time-varying LAI, and save 2000-2023
# cd /g/data/w97/mm3972/scripts/Drought/Post_drought_rainfall/process_for_CABLE/process_LAI
# echo "Step 8: using climatology to derive data back to 2000"
# window_vary=31
# window_clim=31
# python process_MODIS_LAI_time_varying_step8_extend_to_2000.py --window_vary "$window_vary" --window_clim "$window_clim"


# # Step 9: remove high frequent varibility
# cd /g/data/w97/mm3972/scripts/Drought/Post_drought_rainfall/process_for_CABLE/process_LAI/remove_high_frequent_varibility
# echo "Step 9: remove high frequent varibility "
# window=15
# python process_MODIS_LAI_time_varying_step9_remove_high_frequent_varibility.py --window "$window"



# # Step 10: to one file
# window=31
# echo "Step 9: to one file"
# cd /g/data/w97/mm3972/scripts/Drought/Post_drought_rainfall/process_for_CABLE/process_LAI
# python process_MODIS_LAI_time_varying_step9_save_yearly_nc.py --window "$window"

# # Step 7: compress ncfile
# # echo "Step 7: compress ncfile"
# #nccompress -r -o -np 8 /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km_daily/MCD15A3H.061_500m_aid0001_LAI_regridded_daily_2002-2024.nc
# #nccompress -r -o -np 8 /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/var_LAI_only/MCD15A3H.061_500m_aid0001_2002-2005_LAI.nc
# #nccompress -r -o -np 8 /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/var_LAI_only/MCD15A3H.061_500m_aid0001_2006-2010_LAI.nc
# #nccompress -r -o -np 8 /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/var_LAI_only/MCD15A3H.061_500m_aid0001_2011-2015_LAI.nc
# # nccompress -r -o -np 8 /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/var_LAI_only/MCD15A3H.061_500m_aid0001_2016-2020_LAI_fix.nc
# # nccompress -r -o -np 8 /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/var_LAI_only/MCD15A3H.061_500m_aid0001_2021-2024_LAI.nc
