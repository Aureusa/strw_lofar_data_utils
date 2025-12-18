#!/bin/bash
# -----------------------------------------------------------------------------------
# Script to get mosaic coverage information for a list of fields in the STRW Cluster
# -----------------------------------------------------------------------------------
# Set paths and parameters below:

# Path to strw_lofar_data_utils on your system (probably looks like /home/username/strw_lofar_data_utils)
export PROJECTPATH=/home/penchev/strw_lofar_data_utils # /path/to/your/strw_lofar_data_utils

# Path to save the output coverage file
export OUTPUTPATH=/home/penchev/strw_lofar_data_utils/mosaic_coverage.csv # /path/to/save/mosaic_coverage.csv

# Whether or not to save the output to a CSV file
export SAVE=true

# Whether or not to print verbose output
export VERBOSE=true

# -----------------------------------------------------------------------------------
VERBOSE_FLAG=""
if [ "$VERBOSE" = false ]; then
    VERBOSE_FLAG="--no-verbose"
fi

SAVE_FLAG=""
if [ "$SAVE" = false ]; then
    SAVE_FLAG="--no-save-csv"
fi

python3 $PROJECTPATH/src/core/get_mosaic_coverage.py --save-path $OUTPUTPATH $VERBOSE_FLAG $SAVE_FLAG
# -----------------------------------------------------------------------------------
