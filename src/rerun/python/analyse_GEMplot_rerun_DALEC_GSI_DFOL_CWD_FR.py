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
# Get model output
mod_pools = {}
# carbon pools
model['Cwoo']=mod.variables['Cwoo'][:,:,0,0]
model['Croo']=mod.variables['Croo'][:,:,0,0]
model['Clit']=mod.variables['Clit'][:,:,0,0]
model['Csom']=mod.variables['Csom'][:,:,0,0]
model['Ccwd']=mod.variables['Ccwd'][:,:,0,0]
model['Clab']=mod.variables['Clab'][:,:,0,0]
model['Cfol']=mod.variables['Cfol'][:,:,0,0]
model['Cbio']=mod.variables['Cbio'][:,:,0,0]
model['time']=mod.variables['time']

# Litter/foliage related components
model['lai']=mod.variables['lai'][:,:,0,0]
model['flux_fol_lit']=mod.variables['flux_fol_lit'][:,:,0,0]
model['flux_root_lit']=mod.variables['flux_root_lit'][:,:,0,0]
model['flux_cwd_lit']=mod.variables['flux_cwd_lit'][:,:,0,0]
model['flux_wood_cwd']=mod.variables['flux_wood_cwd'][:,:,0,0]
model['']=mod.variables[''][:,:,0,0]

# fluxes
model['Reco']=mod.variables['Reco'][:,:,0,0]
model['Rauto']=mod.variables['Rauto'][:,:,0,0]
model['Rhet']=mod.variables['Rhet'][:,:,0,0]
model['Rh_lit']=mod.variables['Rh_lit'][:,:,0,0]
model['gpp']=mod.variables['gpp'][:,:,0,0]
model['nee']=mod.variables['nee'][:,:,0,0]
model['decomp_lit']=mod.variables['decomp_lit'][:,:,0,0]

# GSI
model['gsi']=mod.variables['gsi'][:,:,0,0]
model['gsi_iphoto']=mod.variables['gsi_iphoto'][:,:,0,0]
model['gsi_ivpd']=mod.variables['gsi_ivpd'][:,:,0,0]
model['gsi_itemp']=mod.variables['gsi_itemp'][:,:,0,0]

# Get corresponding observations

obs_pools={}

# plot Carbon stocks
pCAR.plot_carbon_pools_ts(mod_pools,obs_pools)
