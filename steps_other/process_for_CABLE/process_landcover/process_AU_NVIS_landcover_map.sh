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

var_name="LC"
tiff_path="/g/data/w97/mm3972/data/AU_NVIS_landcover/Australia_veg.tif"
netcdf_path="/g/data/w97/mm3972/data/AU_NVIS_landcover/nc_file/AU_NVIS_landcover.nc"
resolution="5km"

# # Step 0: check tiff file info
# echo "Step 0: check tiff file info"
# gdalinfo $tiff_path > gdalinfo_tif_file_info

# Step 1: translate tiff file to netcdf file
echo "Step 1: translate tiff file to netcdf file"
python process_AU_NVIS_landcover_step1_from_tif_to_nc.py --var_name $var_name --tiff_path $tiff_path --netcdf_path $netcdf_path

# # Step 2: regrid NVIS landcover to CABLE PFT with C3/C4 classification to 5 km
# echo "Step 4: regrid to 5 km"
# cd /g/data/w97/mm3972/scripts/Drought/Post_drought_rainfall/process_for_CABLE/process_landcover
# python process_MODIS_landcover_step4_regrid_to_5km.py
