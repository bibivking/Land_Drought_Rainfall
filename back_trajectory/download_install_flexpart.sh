#!/bin/bash

# -----------------------------
# Set installation directories
# -----------------------------
INSTALL_DIR="/srv/ccrc/LandAP/z5218916/script/Land_Drought_Rainfall/back_trajectory/flexpart"
FLEXPART_VERSION="flexpart_v10.4"

# -----------------------------
# Create Conda Environment
# -----------------------------
ENV_NAME="flexpart_env"

# Create conda env with required packages (you can adjust if needed)
conda create -n $ENV_NAME -y gfortran_linux-64 make wget curl libnetcdf

# Activate the environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# Check if conda environment is activated
if [ "$CONDA_DEFAULT_ENV" != "$ENV_NAME" ]; then
    echo "Failed to activate conda environment '$ENV_NAME'. Exiting."
    exit 1
fi

# -----------------------------
# Download and extract FLEXPART
# -----------------------------
mkdir -p $INSTALL_DIR
cd $INSTALL_DIR
wget https://www.flexpart.eu/downloads/71 -O ${FLEXPART_VERSION}.tar.gz

tar -xvzf ${FLEXPART_VERSION}.tar.gz
cd $FLEXPART_VERSION

# -----------------------------
# Set environment variables for build
# -----------------------------
export FLEXPART_HOME=$INSTALL_DIR/$FLEXPART_VERSION
export FC=gfortran
export NETCDF=$(conda run -n $ENV_NAME nc-config --prefix)  # Get NetCDF path from Conda

# Check if NetCDF path exists
if [ ! -d "$NETCDF" ]; then
    echo "NetCDF not found in conda environment. Exiting."
    exit 1
fi

# -----------------------------
# Build FLEXPART
# -----------------------------
make -j4  # Use 4 cores for faster build

# -----------------------------
# Add FLEXPART to PATH for convenience
# -----------------------------
echo "conda activate $ENV_NAME" >> ~/.bashrc
echo "export PATH=$FLEXPART_HOME:\$PATH" >> ~/.bashrc

echo "FLEXPART installed successfully in $FLEXPART_HOME"
echo "Every time you want to use FLEXPART, run: conda activate $ENV_NAME"

# -----------------------------
# Reminder
# -----------------------------
echo ""
echo "✅ All done! To start using FLEXPART, run:"
echo "    conda activate $ENV_NAME"
echo "    cd $FLEXPART_HOME"
echo "    ./FLEXPART"

