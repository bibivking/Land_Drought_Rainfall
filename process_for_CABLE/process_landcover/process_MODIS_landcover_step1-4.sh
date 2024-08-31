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

# =========== Setting ===========
# LC_Type1: Annual International Geosphere-Biosphere Programme (IGBP) classification
# LC_Type5:	Annual Plant Functional Types classification
lc_varname="LC_Type1"

# =========== Processing ===========
# # Step 1: to Gregorian calendar since cdo inttime cannot process Julian calendar
# cd /g/data/w97/mm3972/data/MODIS/MODIS_landcover
# ncatted -O -a calendar,time,m,c,gregorian MCD12Q1.061_500m_aid0001_2001-2005.nc
# ncatted -O -a calendar,time,m,c,gregorian MCD12Q1.061_500m_aid0001_2006-2010.nc
# ncatted -O -a calendar,time,m,c,gregorian MCD12Q1.061_500m_aid0001_2011-2015.nc
# ncatted -O -a calendar,time,m,c,gregorian MCD12Q1.061_500m_aid0001_2016-2020.nc
# ncatted -O -a calendar,time,m,c,gregorian MCD12Q1.061_500m_aid0001_2021-2022.nc

# # Step 2: to one file
# cd /g/data/w97/mm3972/data/MODIS/MODIS_landcover
# echo 'Step 2: to one file'
# cdo mergetime MCD12Q1.061_500m_aid0001_????-????.nc MCD12Q1.061_500m_aid0001_2001-2022.nc

# # Step 3: extract landcover [LC_Type1 or LC_Type5] from files
# if [ "$lc_varname" == "LC_Type1" ]; then
#     lc_name="IGBP"
# elif [ "$lc_varname" == "LC_Type5" ]; then
#     lc_name="PFT"
# fi
# cd /g/data/w97/mm3972/data/MODIS/MODIS_landcover
# echo "Step 3: extract landcover ${lc_name} from files"
# cdo selname,$lc_varname /g/data/w97/mm3972/data/MODIS/MODIS_landcover/MCD12Q1.061_500m_aid0001_2001-2022.nc "/g/data/w97/mm3972/data/MODIS/MODIS_landcover/var_landcover_only/MCD12Q1.061_500m_aid0001_2001-2022_${lc_name}.nc"

# Step 4: regrid IGBP/PFT to CABLE PFT and to 5 km
echo "Step 4: regrid to 5 km"
cd /g/data/w97/mm3972/scripts/Drought/Post_drought_rainfall/process_for_CABLE/process_landcover
python process_MODIS_landcover_step4_regrid_to_5km.py

# Step 5: resampling
# cdo remapbil,target_grid.txt /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/var_LAI_only/MCD15A3H.061_500m_aid0001_2002-2005_LAI.nc /g/data/w97/mm3972/data/MODIS/MODIS_LAI/AUS/regrid_2_AWAP_5km/MCD15A3H.061_500m_aid0001_2002-2005_LAI_regridded.nc
