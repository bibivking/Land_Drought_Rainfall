bash

# Case study name, used in wrf src file, etc
# edit in input_scripts.txt, 


# Set experiment, DYN, DYN_SM, CLIM
export MMY_EXP_NAME="DYN"

export MMY_ENSEMBLE_NAME="A1"

export MMY_CASE_NAME="test_runs"

# Set the name of exemplary file
export MMY_EXEMPLARY_FILE_NAME="test_runs"

# Starting year of WRF run 
export MMY_SYEAR="2006" 

# Starting month of WRF run
export MMY_SMONTH="12"

# Starting day of WRF run
export MMY_SDAY="28"

# Ending year of WRF run
export MMY_NYEAR="2007" 

# Ending month of WRF run
export MMY_NMONTH="01"

# Ending day of WRF run
export MMY_NDAY="02"

### Domain setting
# Edit for each case study in namelist.wps
export MMY_I_PARENT_START_AR="1, 84, 41"
export MMY_J_PARENT_START_AR="1, 14, 63"
export MMY_E_WE_AR="187, 271, 607"
export MMY_E_SN_AR="172, 253, 505"


# The same for each study [Don't change]
export MMY_STAND_LON=""
export MMY_POLE_LON=""
export MMY_POLE_LAT=""
export MMY_REF_LAT=""
export MMY_REF_LON=""
export MMY_TRUELAT1=""
export MMY_TRUELAT2=""

### WRF physics
export MMY_mp_physics="4, 4, 4,"
export MMY_ra_lw_physics="5, 5, 5,"
export MMY_ra_sw_physics="5, 5, 5,"
export MMY_bl_pbl_physics="1, 1, 1,"
export MMY_cu_physics="1, 1, 0,"

### SET UP LIS_config ### 

# set up case file in Land_Drought_Rainfall_wrf_runs
mkdir -p "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Land_Drought_Rainfall_wrf_runs/${MMY_CASE_NAME}"

# set up ensemble folder in Land_Drought_Rainfall_wrf_runs
mkdir -p "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Land_Drought_Rainfall_wrf_runs/${MMY_CASE_NAME}/${MMY_EXP_NAME}_${MMY_ENSEMBLE_NAME}"

# copy WRF src code
cp -r /g/data/w97/mm3972/model/wrf/NUWRF/nuwrf_drought_breaking_rainfall/installs/time_varying_SM_LAI_albedo NUWRF_exe

# copy nu-wrf_scripts from exemplary file
cp "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Land_Drought_Rainfall_wrf_runs/${MMY_EXEMPLARY_FILE_NAME}/set_up_nu-wrf_scripts.sh" .
./set_up_nu-wrf_scripts.sh

# edit deck



# Run prepare simulation
cd /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Land_Drought_Rainfall_wrf_runs/${MMY_CASE_NAME}/nu-wrf_scripts/prepare_decks/codes
./run_prepare_simulation.sh


# Run geogrid
cd /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Land_Drought_Rainfall_wrf_runs/${MMY_CASE_NAME}/nu-wrf_scripts/prepare_decks/templates
qsub run_geogrid.sh

