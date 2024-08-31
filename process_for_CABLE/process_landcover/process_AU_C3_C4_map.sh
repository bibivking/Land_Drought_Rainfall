#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=5:00:00
#PBS -l mem=300GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97

module use /g/data/hh5/public/modules
module load conda/analysis3

var_name="C4_proportion"
tiff_path="/g/data/w97/mm3972/data/AU_C3_C4_plant_map/Proportional_C4_vegation.tif"
netcdf_path="/g/data/w97/mm3972/data/AU_C3_C4_plant_map/nc_file/Proportional_C4_vegation.nc"
output_path="/g/data/w97/mm3972/data/AU_C3_C4_plant_map/nc_file/Proportional_C4_vegation_NVIS_res.nc"
resolution="NVIS_res"

# # Step 0: check tiff file info
# echo "Step 0: check tiff file info"
# gdalinfo $tiff_path > gdalinfo_tif_file_info

# # Step 1: translate tiff file to netcdf file
# echo "Step 1: translate tiff file to netcdf file"
# python process_AU_C3_C4_step1_from_tif_to_nc.py --var_name $var_name --tiff_path $tiff_path --netcdf_path $netcdf_path

# # Step 2: out put AWAP resolution file's lat lon information
# Option 1: Get AWAP lat-lon info
# cdo griddes /g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc > target_grid.txt
# Option 2: Get NVIS lat-lon info
# cdo griddes /g/data/w97/mm3972/data/AU_NVIS_landcover/nc_file/AU_NVIS_landcover.nc > AU_NVIS_landcover_grid.txt

# # # Step 3: resampling
# echo "Step 3: resampling, Method 1: python script to MODIS res"
# python process_AU_C3_C4_step3_regrid_to_coarse_res.py --var_name $var_name --netcdf_path $netcdf_path --output_path $output_path --resolution $resolution
# # echo "Step 3: resampling by bilinear method, Method 2: use cdo to 5km res"
# # cdo remapbil,target_grid.txt $netcdf_path $output_path
#echo "Step 3: resampling by nearest method, Method 3: use cdo to NVIS res"
#cdo remapnn,AU_NVIS_landcover_grid.txt $netcdf_path $output_path

# Step 4: make dominant veg type (C3 or C4) map
echo "Step 4: make dominant veg type (C3 or C4) map"
python process_AU_C3_C4_step4_dominant_type.py --resolution $resolution
