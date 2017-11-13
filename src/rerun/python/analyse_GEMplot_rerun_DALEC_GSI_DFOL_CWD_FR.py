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

# Project data file
data_dir = "/home/dmilodow/DataStore_DTM/BALI/CARDAMOM_BALI/npydata/"
project_met = "BALI_GEMplots_daily_drivers.npy"
project_par = "BALI_GEMplots_daily_params.npy"
project_obs = "BALI_GEMplots_daily_obs.npy"

project = 'BALI_GEMplots_daily'
run = '001'
filename = 'BALI_GEMplots_daily_2011_2016_DALEC_GSI_DFOL_CWD_FR.nc'

# find NetCDF_file for rerun and load
NetCDF_file = '%s%s/rerun/%s/%s' % (path2project,project,run,filename)
mod = Dataset(NetCDF_file)

site = 0
# Get model output
model = {}
# carbon pools
model['Cwoo']=mod.variables['Cwoo'][:,:,site]
model['Croo']=mod.variables['Croo'][:,:,site]
model['Clit']=mod.variables['Clit'][:,:,site]
model['Csom']=mod.variables['Csom'][:,:,site]
model['Ccwd']=mod.variables['Ccwd'][:,:,site]
model['Clab']=mod.variables['Clab'][:,:,site]
model['Cfol']=mod.variables['Cfol'][:,:,site]
model['Cbio']=mod.variables['Cbio'][:,:,site]
model['time']=mod.variables['time']

# Litter/foliage related components
model['lai']=mod.variables['lai'][:,:,site]
model['flux_fol_lit']=mod.variables['flux_fol_lit'][:,:,site]
model['flux_root_lit']=mod.variables['flux_root_lit'][:,:,site]
model['flux_cwd_lit']=mod.variables['flux_cwd_lit'][:,:,site]
model['flux_wood_cwd']=mod.variables['flux_wood_cwd'][:,:,site]

# fluxes
model['Reco']=mod.variables['Reco'][:,:,site]
model['Rauto']=mod.variables['Rauto'][:,:,site]
model['Rhet']=mod.variables['Rhet'][:,:,site]
model['Rh_lit']=mod.variables['Rh_lit'][:,:,site]
model['gpp']=mod.variables['gpp'][:,:,site]
model['nee']=mod.variables['nee'][:,:,site]
model['decomp_lit']=mod.variables['decomp_lit'][:,:,site]

# GSI
model['gsi']=mod.variables['gsi'][:,:,site]
model['gsi_iphoto']=mod.variables['gsi_iphoto'][:,:,site]
model['gsi_ivpd']=mod.variables['gsi_ivpd'][:,:,site]
model['gsi_itemp']=mod.variables['gsi_itemp'][:,:,site]

# Get corresponding observations
obs_in = np.load('%s%s' % (data_dir,project_obs))
obs_in[obs_in==-9999]=np.nan
obs={}
obs['time'] = mod.variables['time']
obs['gpp'] = obs_in[site,:,0]
obs['gpp_u'] = obs_in[site,:,11]
obs['nee'] = obs_in[site,:,2]
obs['nee_u'] = obs_in[site,:,13]
obs['lai'] = obs_in[site,:,1]
obs['lai_u'] = obs_in[site,:,12]
obs['Cwoo'] = obs_in[site,:,6]
obs['Cwoo_u'] = obs_in[site,:,17]
obs['Cfol'] = obs_in[site,:,5]
obs['Cfol_u'] = obs_in[site,:,16]
obs['Croo'] = obs_in[site,:,7]
obs['Croo_u'] = obs_in[site,:,18]
obs['Clit'] = obs_in[site,:,8]
obs['Clit_u'] = obs_in[site,:,19]
obs['Csom'] = obs_in[site,:,9]
obs['Csom_u'] = obs_in[site,:,20]
obs['Reco'] = obs_in[site,:,4]
obs['Reco_u'] = obs_in[site,:,15]
obs['flux_fol_lit'] = obs_in[site,:,32]
obs['flux_fol_lit_u'] = obs_in[site,:,33]

# plot Carbon stocks
pCAR.plot_carbon_pools_ts(model,obs)
pCAR.plot_litter_components_ts(model,obs)
