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

import sys
sys.path.append('/home/dmilodow/DataStore_DTM/BALI/CARDAMOM_BALI/CARDAMOM_DTM/src/rerun/python')

import plot_CARDAMOM_output as pCAR

path2project = '/home/dmilodow/DataStore_DTM/BALI/CARDAMOM_BALI/projects/'
project = 'BALI_GEMplots_daily'
run = '001'
filename = 'BALI_GEMplots_daily_2011_2016_DALEC_GSI_DFOL_CWD_FR.nc'

# find NetCDF_file and load
NetCDF_file = '%s%s/rerun/%s/%s' % (path2project,project,run,filename)
mod = Dataset(NetCDF_file)

site_index = 0
mod_pools = {}
mod_pools['Cwoo']=mod.variables['Cwoo'][:,:,0,0]
mod_pools['Croo']=mod.variables['Croo'][:,:,0,0]
mod_pools['Clit']=mod.variables['Clit'][:,:,0,0]
mod_pools['Csom']=mod.variables['Csom'][:,:,0,0]
mod_pools['Ccwd']=mod.variables['Ccwd'][:,:,0,0]
mod_pools['Clab']=mod.variables['Clab'][:,:,0,0]
mod_pools['Cfol']=mod.variables['Cfol'][:,:,0,0]
mod_pools['time']=mod.variables['time']

obs_pools={}

# plot Carbon stocks
pCAR.plot_carbon_pools_ts(mod_pools,obs_pools)
