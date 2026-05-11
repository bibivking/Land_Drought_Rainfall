#!/bin/bash
module purge
module use /g/data/hh5/public/modules
module load conda/analysis3-24.04
module load fftw3/3.3.10
module load netcdf/4.7.3
module load hdf5/1.12.2
module load gcc/14.2.0
module load eccodes/2.34.1

#export PATH=/home/561/mm3972/.local/lib/python3.6/site-packages/genshi:/home/561/mm3972/.local/lib/python3.6/site-packages/eccodes:$PATH 

#export PYTHONPATH=/home/561/mm3972/.local/lib/python3.10/site-packages:$PYTHONPATH

cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/back_trajectory/flex_extract
#python mmy_test_web_api.py
#python mmy_test_cds_api_2.py
./setup_local.sh
