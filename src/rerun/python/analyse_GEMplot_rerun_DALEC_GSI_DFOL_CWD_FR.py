#===============================================================================
# analyse_GEMplot_rerun_DALEC_GSI_DFOL_CWD_FR.py
#-------------------------------------------------------------------------------
# This script analyses output from the GEM plot CARDAMOM runs using the model
# DALEC_GSI_DFOL_CWD_FR. Analysis includes:
# - diagnostic plots against GEM plot survey data
# - (to be added as required)
#-------------------------------------------------------------------------------
# Author: D.T.Milodowski
# Date: October 2017
#===============================================================================
import numpy as np
from matplotlib import pyplot as plt
from netCDF4 import Dataset

path2project = '/home/dmilodow/DataStore_DTM/BALI/CARDAMOM_BALI/projects/'
project = 'BALI_GEMplots_daily'
run = '001'

# find NetCDF_file and load
NetCDF_file = '%s%s/rerun/%s' % (path2project,project,run)
dataset = Dataset(NetCDF_file)

# Get the spatial information from the layer
Lat = dataset.variables[lat_var][:]
Long = dataset.variables[lon_var][:]

