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
run = '025'
filename = 'BALI_GEMplots_daily_2011_2017_DALEC_GSI_DFOL_CWD_FR.nc'

# find NetCDF_file for rerun and load
NetCDF_file = '%s%s/rerun/%s/%s' % (path2project,project,run,filename)
mod = Dataset(NetCDF_file)

# Get corresponding observations
obs_in = np.load('%s%s/%s' % (data_dir,run,project_obs))
obs_in[obs_in==-9999]=np.nan

n_sites = 6
sites = ['MLA01','MLA02','SAF04','SAF05','SAF02','SAF01']
params = {}
model = {}
obs = {}
start_date = np.datetime64('2011-01-01','D')
for i in range(0,n_sites):
    model[sites[i]] = {}
    obs[sites[i]] = {}
    
    param_file = '%s%s/rerun/%s/0000%i.csv' % (path2project,project,run,i+1)
    ll_file = '%s%s/rerun/%s/0000%i_ll.csv' % (path2project,project,run,i+1)
    params[sites[i]] = np.genfromtxt(param_file,skiprows = 1,delimiter=',')[:,1:]
    ll=np.genfromtxt(param_file,skiprows = 1,delimiter=',')[1,1:]
    print sites[i],params.shape[sites[i]],ll.shape
    
    # Get model output
    # carbon pools
    model[sites[i]]['Cwoo']=mod.variables['Cwoo'][:,:,i]
    model[sites[i]]['Croo']=mod.variables['Croo'][:,:,i]
    model[sites[i]]['Clit']=mod.variables['Clit'][:,:,i]
    model[sites[i]]['Csom']=mod.variables['Csom'][:,:,i]
    model[sites[i]]['Ccwd']=mod.variables['Ccwd'][:,:,i]
    model[sites[i]]['Clab']=mod.variables['Clab'][:,:,i]
    model[sites[i]]['Cfol']=mod.variables['Cfol'][:,:,i]
    model[sites[i]]['Cbio']=mod.variables['Cbio'][:,:,i]
    model[sites[i]]['time']=start_date+(mod.variables['time'][:]-1)
    
    # Litter/foliage related components
    model[sites[i]]['lai']=mod.variables['lai'][:,:,i]
    model[sites[i]]['flux_fol_lit']=mod.variables['flux_fol_lit'][:,:,i]
    model[sites[i]]['flux_root_lit']=mod.variables['flux_root_lit'][:,:,i]
    model[sites[i]]['flux_cwd_lit']=mod.variables['flux_cwd_lit'][:,:,i]
    model[sites[i]]['flux_wood_cwd']=mod.variables['flux_wood_cwd'][:,:,i]
    model[sites[i]]['ll']= ll.copy()
    
    # fluxes
    model[sites[i]]['Reco']=mod.variables['Reco'][:,:,i]
    model[sites[i]]['Rauto']=mod.variables['Rauto'][:,:,i]
    model[sites[i]]['Rhet']=mod.variables['Rhet'][:,:,i]
    model[sites[i]]['Rh_lit']=mod.variables['Rh_lit'][:,:,i]
    model[sites[i]]['gpp']=mod.variables['gpp'][:,:,i]
    model[sites[i]]['npp']=mod.variables['npp'][:,:,i]
    model[sites[i]]['nee']=mod.variables['nee'][:,:,i]
    model[sites[i]]['decomp_lit']=mod.variables['decomp_lit'][:,:,i]
    
    # GSI
    model[sites[i]]['gsi']=mod.variables['gsi'][:,:,i]
    model[sites[i]]['gsi_iphoto']=mod.variables['gsi_iphoto'][:,:,i]
    model[sites[i]]['gsi_ivpd']=mod.variables['gsi_ivpd'][:,:,i]
    model[sites[i]]['gsi_itemp']=mod.variables['gsi_itemp'][:,:,i]
    
    obs[sites[i]]['time'] =start_date+(mod.variables['time'][:]-1)
    obs[sites[i]]['gpp'] = obs_in[i,:,0]
    obs[sites[i]]['gpp_u'] = obs_in[i,:,11]
    obs[sites[i]]['nee'] = obs_in[i,:,2]
    obs[sites[i]]['nee_u'] = obs_in[i,:,13]
    obs[sites[i]]['lai'] = obs_in[i,:,1]
    obs[sites[i]]['lai_u'] = obs_in[i,:,12]
    obs[sites[i]]['Cwoo'] = obs_in[i,:,6]
    obs[sites[i]]['Cwoo_u'] = obs_in[i,:,17]
    obs[sites[i]]['Cfol'] = obs_in[i,:,5]
    obs[sites[i]]['Cfol_u'] = obs_in[i,:,16]
    obs[sites[i]]['Croo'] = obs_in[i,:,7]
    obs[sites[i]]['Croo_u'] = obs_in[i,:,18]
    #obs[sites[i]]['Clit'] = obs_in[i,:,8]
    #obs[sites[i]]['Clit_u'] = obs_in[i,:,19]
    obs[sites[i]]['Csom'] = obs_in[i,:,9]
    obs[sites[i]]['Csom_u'] = obs_in[i,:,20]
    obs[sites[i]]['Reco'] = obs_in[i,:,4]
    obs[sites[i]]['Reco_u'] = obs_in[i,:,15]
    obs[sites[i]]['flux_fol_lit'] = obs_in[i,:,8]
    obs[sites[i]]['flux_fol_lit_u'] = obs_in[i,:,19]
    obs[sites[i]]['lit_acc_days'] = obs_in[i,:,28]

    # plot Carbon stocks
    #pCAR.plot_carbon_pools_ts(model[sites[i]],obs[sites[i]],figname='carbon_pools_ts_%s_%s.png' % (run, sites[i]))
    #pCAR.plot_carbon_fluxes_ts(model[sites[i]],obs[sites[i]],figname='carbon_fluxes_ts_%s_%s.png' % (run,sites[i]))
    #pCAR.plot_litter_components_ts(model[sites[i]],obs[sites[i]],figname='litter_components_ts_%s_%s.png' % (run,sites[i]))
    #pCAR.plot_litter_trap_comparison(model[sites[i]],obs[sites[i]],figname='litter_trap_components_%s_%s.png' % (run,sites[i]))
    #pCAR.plot_parameters(params[sites[i]],figname='parameters_%s_%s.png' % (run,sites[i]))
    pCAR.plot_summary_ts(model[sites[i]],obs[sites[i]],figname='summary_ts_%s_%s.png' % (run,sites[i]))
    plt.show()
    plt.close('all')
    
# Summarise the plots with temporal average and lower and upper quartiles
Cwoo = np.zeros((6,3))
Croo = np.zeros((6,3))
Cfol = np.zeros((6,3))
Csom = np.zeros((6,3))
Clit = np.zeros((6,3))
Clab = np.zeros((6,3))

GPP = np.zeros((6,3))
NPP = np.zeros((6,3))
NEE = np.zeros((6,3))
Litterfall = np.zeros((6,3))
Raut = np.zeros((6,3))
Rhet = np.zeros((6,3))
LAI = np.zeros((6,3))

Sites = ["MLA01","MLA02","SAF04","SAF05","SAF02","SAF01"]
    
Cwoo[:,0] = np.mean(mod.variables['Cwoo'][:,1,:],axis=0)
Croo[:,0] = np.mean(mod.variables['Croo'][:,1,:],axis=0)
Cfol[:,0] = np.mean(mod.variables['Cfol'][:,1,:],axis=0)
Csom[:,0] = np.mean(mod.variables['Csom'][:,1,:],axis=0)
Clit[:,0] = np.mean(mod.variables['Clit'][:,1,:],axis=0)
Clab[:,0] = np.mean(mod.variables['Clab'][:,1,:],axis=0)

GPP[:,0] = np.mean(mod.variables['gpp'][:,0,:],axis=0)
NPP[:,0] = np.mean(mod.variables['npp'][:,0,:],axis=0)
NEE[:,0] = np.mean(mod.variables['nee'][:,0,:],axis=0)
Litterfall[:,0] = np.mean(mod.variables['flux_fol_lit'][:,0,:],axis=0)
Raut[:,0] = np.mean(mod.variables['Rauto'][:,0,:],axis=0)
Rhet[:,0] = np.mean(mod.variables['Rhet'][:,0,:],axis=0)
LAI[:,0] = np.mean(mod.variables['lai'][:,0,:],axis=0)

Cwoo[:,1] = np.mean(mod.variables['Cwoo'][:,-2,:],axis=0)
Croo[:,1] = np.mean(mod.variables['Croo'][:,-2,:],axis=0)
Cfol[:,1] = np.mean(mod.variables['Cfol'][:,-2,:],axis=0)
Csom[:,1] = np.mean(mod.variables['Csom'][:,-2,:],axis=0)
Clit[:,1] = np.mean(mod.variables['Clit'][:,-2,:],axis=0)
Clab[:,1] = np.mean(mod.variables['Clab'][:,-2,:],axis=0)

GPP[:,1] = np.mean(mod.variables['gpp'][:,-2,:],axis=0)
NPP[:,1] = np.mean(mod.variables['npp'][:,-2,:],axis=0)
NEE[:,1] = np.mean(mod.variables['nee'][:,-2,:],axis=0)
Litterfall[:,1] = np.mean(mod.variables['flux_fol_lit'][:,-2,:],axis=0)
Raut[:,1] = np.mean(mod.variables['Rauto'][:,-2,:],axis=0)
Rhet[:,1] = np.mean(mod.variables['Rhet'][:,-2,:],axis=0)
LAI[:,1] = np.mean(mod.variables['lai'][:,-2,:],axis=0)

Cwoo[:,2] = np.mean(mod.variables['Cwoo'][:,-1,:],axis=0)
Croo[:,2] = np.mean(mod.variables['Croo'][:,-1,:],axis=0)
Cfol[:,2] = np.mean(mod.variables['Cfol'][:,-1,:],axis=0)
Csom[:,2] = np.mean(mod.variables['Csom'][:,-1,:],axis=0)
Clit[:,2] = np.mean(mod.variables['Clit'][:,-1,:],axis=0)
Clab[:,2] = np.mean(mod.variables['Clab'][:,-1,:],axis=0)

GPP[:,2] = np.mean(mod.variables['gpp'][:,-1,:],axis=0)
NPP[:,2] = np.mean(mod.variables['npp'][:,-1,:],axis=0)
NEE[:,2] = np.mean(mod.variables['nee'][:,-1,:],axis=0)
Litterfall[:,2] = np.mean(mod.variables['flux_fol_lit'][:,-1,:],axis=0)
Raut[:,2] = np.mean(mod.variables['Rauto'][:,-1,:],axis=0)
Rhet[:,2] = np.mean(mod.variables['Rhet'][:,-1,:],axis=0)
LAI[:,2] = np.mean(mod.variables['lai'][:,-1,:],axis=0)

# Write report to file
out = open('BALI_stocks_summary_%s.csv' % run,'w')
# First set up header
out.write('Site, Cwoo_50, Cwoo_25, Cwoo_75, Croo_50, Croo_25, Croo_75,  Cfol_50, Cfol_25, Cfol_75, Clab_50, Clab_25, Clab_75, Clit_50, Clit_25, Clit_75, Csom_50, Csom_25, Csom_75\n')
for i in range(0,6):
    out.write('%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n' % (Sites[i],Cwoo[i,0],Cwoo[i,1],Cwoo[i,2],Croo[i,0],Croo[i,1],Croo[i,2],Cfol[i,0],Cfol[i,1],Cfol[i,2],Clab[i,0],Clab[i,1],Clab[i,2],Clit[i,0],Clit[i,1],Clit[i,2],Csom[i,0],Csom[i,1],Csom[i,2]))
out.close()  


    
out = open('BALI_fluxes_summary_%s.csv' % run,'w')
# First set up header
out.write('Site, GPP_50, GPP_25, GPP_75, NPP_50, NPP_25, NPP_75, NEE_50, NEE_25, NEE_75,  Ra_50, Ra_25, Ra_75, Rh_50, Rh_25, Rh_75, Litterfall_50, Litterfall_25, Litterfall_75, LAI_50, LAI_25, LAI_75\n')
for i in range(0,6):
    out.write('%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n' % (Sites[i],GPP[i,0],GPP[i,1],GPP[i,2],NPP[i,0],NPP[i,1],NPP[i,2],NEE[i,0],NEE[i,1],NEE[i,2],Raut[i,0],Raut[i,1],Raut[i,2],Rhet[i,0],Rhet[i,1],Rhet[i,2],Litterfall[i,0],Litterfall[i,1],Litterfall[i,2],LAI[i,0],LAI[i,1],LAI[i,2]))
out.close()  


## Violin plots
# pull out allocation fractions
import pandas

# allocation fractions
n_params = 0

for pp in range(0,n_sites):
    n_params += params[sites[pp]].shape[1]
    
all_wood = np.zeros(n_params)*np.nan    
all_canopy = np.zeros(n_params)*np.nan
all_root = np.zeros(n_params)*np.nan
rt_wood = np.zeros(n_params)*np.nan    
rt_fol = np.zeros(n_params)*np.nan
rt_root = np.zeros(n_params)*np.nan
LMA = np.zeros(n_params)*np.nan
Narea = np.zeros(n_params)*np.nan
CN = np.zeros(n_params)*np.nan
plot_array = np.empty(n_params,dtype='S5')
ff=0
for pp in range(0,n_sites):
    
    ii=ff
    ff=ii+params[sites[pp]].shape[1]
    all_root[ii:ff]=params[sites[pp]][3,:].copy()
    all_canopy[ii:ff]=params[sites[pp]][12,:]*(1-params[sites[pp]][3,:])
    all_wood[ii:ff]=(1-params[sites[pp]][12,:])*(1-params[sites[pp]][3,:])
    LMA[ii:ff]=params[sites[pp]][16,:].copy()
    Narea[ii:ff]=10**params[sites[pp]][10,:]
    CN[ii:ff]=LMA[ii:ff]/Narea[ii:ff]
    plot_array[ii:ff]=sites[pp]
    rt_fol[ii:ff]=model[sites[pp]]['ll']
    rt_wood[ii:ff]=1/params[sites[pp]][5,:]/365.
    rt_root[ii:ff]=1/params[sites[pp]][6,:]/365.
    

#put into pandas data frame
header_names=['root','canopy','wood','plot','LMA','Narea','CNratio','wood_rt','fol_rt','root_rt']
combined_array = np.asarray((all_root,all_canopy,all_wood,plot_array,LMA,Narea,CN,rt_wood,rt_fol,rt_root)).transpose()
df = pandas.DataFrame(data=combined_array,columns=header_names)
for col in df.keys():
    if col != 'plot':
        df[col] = df[col].astype('float64')
    else:
        df[col] = df[col].astype('category')
        
pCAR.plot_allocation_fractions(df)
pCAR.plot_leaf_traits(df)
pCAR.plot_residence_times(df)
