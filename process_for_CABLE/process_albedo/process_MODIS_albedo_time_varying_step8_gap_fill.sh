#!/bin/bash
# Call the Python script with the file path and window size
file_dir="/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km_daily/"
albedo_band="WSA_nir"
window=30

# # Loop through all files in the path
# # for year in {2003..2023}; do
# for year in {2001..2002}; do
#     pre_year=$((year-1))
#     aft_year=$((year+1))
#     echo "$year"
#     echo "$pre_year"
#     echo "$aft_year"
year=2000
aft_year=2001
pre_year=1999

cat > process_MODIS_albedo_time_varying_step8_gap_fill_${albedo_band}_${year}.sh << EOF_process_albedo
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

cp process_MODIS_albedo_time_varying_step8_gap_fill.py process_MODIS_albedo_time_varying_step8_gap_fill_${albedo_band}_${year}.py
if [ $year -eq 2000 ]; then
    python process_MODIS_albedo_time_varying_step8_gap_fill_${albedo_band}_${year}.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_Feb${year}-Jan${aft_year}_albedo_regridded_daily.nc" --albedo_band "$albedo_band" --window "$window"
else
    python process_MODIS_albedo_time_varying_step8_gap_fill_${albedo_band}_${year}.py "${file_dir}MCD43A3.061_500m_aid0001_${albedo_band}_Dec${pre_year}-Jan${aft_year}_albedo_regridded_daily.nc" --albedo_band "$albedo_band" --window "$window"
fi
rm process_MODIS_albedo_time_varying_step8_gap_fill_${albedo_band}_${year}.py

EOF_process_albedo

cd /g/data/w97/mm3972/scripts/Drought/Post_drought_rainfall/process_for_CABLE/process_albedo

qsub process_MODIS_albedo_time_varying_step8_gap_fill_${albedo_band}_${year}.sh

# done
