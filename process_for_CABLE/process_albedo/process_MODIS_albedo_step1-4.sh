#!/bin/bash

# Set the path to search
albedo_band="WSA_vis"
MODIS_albedo_path="/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/${albedo_band}/"

# Loop through all files in the path
for file in $(find $MODIS_albedo_path -type f -name "*.nc"); do

  # Extract the file name
  file_name=$(basename "$file")

  # Extract the year from the file name
  year="${file_name##*_}"
  year="${year%%.*}"

  # Print the site name to the console
  echo "$year"

cd /g/data/w97/mm3972/scripts/Drought/Post_drought_rainfall/process_for_CABLE/process_albedo/

cat > process_albedo_step1-4_${albedo_band}_${year}.sh << EOF_process_albedo
#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalbw
#PBS -l walltime=4:00:00
#PBS -l mem=100GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97

module use /g/data/hh5/public/modules
module load conda/analysis3

cd /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/${albedo_band}

# Step 1: to Gregorian calendar since cdo inttime cannot process Julian calendar
# ncatted -O -a calendar,time,m,c,gregorian MCD43A3.061_500m_aid0001_${albedo_band}_${year}.nc

# Step 2: out put AWAP resolution file's lat lon information
# cdo griddes /g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc > target_grid.txt

# Step 3: extract albedo from files
# cdo selname,Albedo_${albedo_band} /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/${albedo_band}/MCD43A3.061_500m_aid0001_${albedo_band}_${year}.nc /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/var_albedo_only/MCD43A3.061_500m_aid0001_${albedo_band}_${year}_albedo.nc

# Step 4: resampling to 5km
cd /g/data/w97/mm3972/scripts/Drought/Post_drought_rainfall/process_for_CABLE/process_albedo
cdo remapbil,target_grid.txt /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/var_albedo_only/MCD43A3.061_500m_aid0001_${albedo_band}_${year}_albedo.nc /g/data/w97/mm3972/data/MODIS/MODIS_Albedo/AUS/regrid_2_AWAP_5km/MCD43A3.061_500m_aid0001_${albedo_band}_${year}_albedo_regridded.nc

EOF_process_albedo

cd /g/data/w97/mm3972/scripts/Drought/Post_drought_rainfall/process_for_CABLE/process_albedo

qsub process_albedo_step1-4_${albedo_band}_${year}.sh

done
