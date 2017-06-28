# This script sets up the CARDAMOM simulations for the BALI GEM plot data.
# It is set up to run DALEC
# DTM 27/06/2017

import numpy as np
import datetime as dt
import os
from CARDAMOM_DTM import *
import load_data as data

# Project data file (to load in if already generated
data_dir = "/home/dmilodow/DataStore_DTM/BALI/CARDAMOM_DTM/npydata/"
project_met_npydata = "BALI_GEMplots_daily_drivers.npy"
project_par_npydata = "BALI_GEMplots_daily_params.npy"
project_obs_npydata = "BALI_GEMplots_daily_obs.npy"

# File list containing drivers, observations and plot coordinates
coordinate_file = "BALI_plot_coordinates.csv"
met_file = "/home/dmilodow/DataStore_DTM/BALI/CSPA_BALI_data_and_analysis/scripts/construct_drivers/BALI_gapfilled_met_station_daily_v1.csv"
par_file = "<SOMETHING>.csv"
obs_file = "<SOMETHING>.csv"

# first load in coordinates
plot, latitude, longitude = data.load_plot_coordinates(coordinate_file)
met_data = data.load_met_data(met_file)
par_data = data.load_par_data(par_file)
dates_obs,obs_data = data.load_obs_data(obs_file)

# start date ### read from data file
d0 = met_data['date'][0]
DoY = met_data['date'] - met_data['date'].astype('timedelta64[Y]')
tstep = met_data['tstep_days']#dates_met - dates_met[0] + np.timedelta[1,'D']
sim_length = tstep[-1]
nosites = len(plot)

#-----------------------------------------------------------------------------
# First deal with the met data 
met = np.zeros((nosites,sim_length,14))

if project_met_npydata not in os.listdir(data_dir):
    
    # First load the met data into the met array - use same driving data for now
    met[pp,:,0] = tstep.copy()                             #  0 = run day
    met[pp,:,1] = met_data['mn2t']                         #  1 = minimum temperature oC
    met[pp,:,2] = met_data['mx2t']                         #  2 = maximum temperature oC
    met[pp,:,3] = met_data['ssrd']                         #  3 = surface shortwave radiation in MJ.m-2.day-1
    met[pp,:,4] = 400.                                     #  4 = atmospheric CO2 concentration ppm
    met[pp,:,5] = DoY.copy()                               #  5 = day of year
    met[pp,:,6] = met_data['pptn']                         #  6 = lagged precipitation
    met[pp,:,7] = -9999                                    #  7 = fire burned deforestation fraction - not applicable 
    met[pp,:,8] = -9999                                    #  8 = fire burned deforestation fraction - not applicable 
    met[pp,:,9] = met_data['mn2t_21d']                     #  9 = 21 day average min temperature K
    met[pp,:,10] = met_data['mx2t_21d']                    # 10 = 21 day average max temperature K
    met[pp,:,11] = met_data['vpd_21d']                     # 11 = 21 day average vpd Pa
    met[pp,:,12] = -9999                                   # 12 = forest management practice to accompany any clearing - not applicable
    met[pp,:,13] = (met_data['mn2t']+met_data['mx2t'])/2.  # 13 = mean temperature oC ???

    np.save(data_dir + project_met_npydata,met)

else:
    print "Loading met drivers"
    met = np.load(data_dir+project_met_npydata)

#-----------------------------------------------------------------------------
# Now deal with the parameters 
parprior = np.zeros((nosites,100))-9999.
parprior_unc = np.zeros((nosites,100))-9999.
otherprior = np.zeros((nosites,100))-9999. # What are these for????
otherprior_unc = np.zeros((nosites,100))-9999.

if project_par_npydata not in os.listdir(data_dir):
    
    for pp in range(0,nosites):
        
        parprior[pp,0] = par_data[plot[pp]]['m_r']      # Litter to SOM conversion rate 
        parprior_unc[pp,0] = par_data[plot[pp]]['m_r_u']
        parprior[pp,1] = par_data[plot[pp]]['f_a']      # Fraction of NPP respired
        parprior_unc[pp,1] = par_data[plot[pp]]['f_a_u']    
        parprior[pp,2] = par_data[plot[pp]]['f_f']      # Fraction of NPP allocated to foliage 
        parprior_unc[pp,2] = par_data[plot[pp]]['f_f_u'] 
        parprior[pp,3] = par_data[plot[pp]]['f_r']      # Fraction of NPP allocated to roots
        parprior_unc[pp,3] = par_data[plot[pp]]['f_r_u']  
        parprior[pp,4] = par_data[plot[pp]]['L_f']      # max leaf turnover (GSI) ! Leaf lifespan (CDEA)
        parprior_unc[pp,4] = par_data[plot[pp]]['L_f_u']  
        parprior[pp,5] = par_data[plot[pp]]['t_w']      # Turnover rate of wood
        parprior_unc[pp,5] = par_data[plot[pp]]['t_w_u'] 
        parprior[pp,6] = par_data[plot[pp]]['t_r']      # Turnover rate of roots
        parprior_unc[pp,6] = par_data[plot[pp]]['t_r_u']  
        parprior[pp,7] = par_data[plot[pp]]['t_l']      # Litter turnover rate
        parprior_unc[pp,7] = par_data[plot[pp]]['t_l_u'] 
        parprior[pp,8] = par_data[plot[pp]]['t_S']      # SOM turnover rate
        parprior_unc[pp,8] = par_data[plot[pp]]['t_S_u']   
        parprior[pp,9] = par_data[plot[pp]]['theta']    # Parameter in exponential term of temperature
        parprior_unc[pp,9] = par_data[plot[pp]]['theta_u']   
        parprior[pp,10] = par_data[plot[pp]]['C_eff']   # Canopy efficiency parameter (part of ACM) 
        parprior_unc[pp,10] = par_data[plot[pp]]['C_eff_u'] 
        parprior[pp,11] = par_data[plot[pp]]['B_day']   # max labile turnover(GSI) ! date of Clab release (CDEA)
        parprior_unc[pp,11] = par_data[plot[pp]]['B_day_u']
        parprior[pp,12] = par_data[plot[pp]]['f_l']     # Fraction allocated to Clab
        parprior_unc[pp,12] = par_data[plot[pp]]['f_l_u']   
        parprior[pp,13] = par_data[plot[pp]]['R_l']     # min temp threshold (GSI) ! lab release duration period (CDEA)
        parprior_unc[pp,13] = par_data[plot[pp]]['R_l_u'] 
        parprior[pp,14] = par_data[plot[pp]]['F_day']   # max temp threshold (GSI)! date of leaf fall
        parprior_unc[pp,14] = par_data[plot[pp]]['F_day_u']
        parprior[pp,15] = par_data[plot[pp]]['mn_photo']# min photoperiod threshold (GSI)
        parprior_unc[pp,15] = par_data[plot[pp]]['mn_photo_u']
        parprior[pp,16] = par_data[plot[pp]]['LMA']     # LMA
        parprior_unc[pp,16] = par_data[plot[pp]]['LMA_u']    
        parprior[pp,17] = par_data[plot[pp]]['Clab_i']  # initial C stock for labile
        parprior_unc[pp,17] = par_data[plot[pp]]['Clab_i_u']
        parprior[pp,18] = par_data[plot[pp]]['Cfol_i']  # initial C stock for foliar
        parprior_unc[pp,18] = par_data[plot[pp]]['Cfol_i_u']
        parprior[pp,19] = par_data[plot[pp]]['Croo_i']  # initial C stock for roots
        parprior_unc[pp,19] = par_data[plot[pp]]['Croo_i_u']
        parprior[pp,20] = par_data[plot[pp]]['Cwoo_i']  # initial C stock for litter
        parprior_unc[pp,20] = par_data[plot[pp]]['Cwoo_i_u'] 
        parprior[pp,21] = par_data[plot[pp]]['Clit_i']  # initial C stock for litter 
        parprior_unc[pp,21] = par_data[plot[pp]]['Clit_i_u']
        parprior[pp,22] = par_data[plot[pp]]['Csom_i']  # initial C stock for som
        parprior_unc[pp,22] = par_data[plot[pp]]['Csom_i_u']
        parprior[pp,23] = par_data[plot[pp]]['mx_photo']# max photoperiod threshold (GSI)
        parprior_unc[pp,23] = par_data[plot[pp]]['mx_photo_u']
        parprior[pp,24] = par_data[plot[pp]]['mn_vpd']  # min VPD threshold (GSI)
        parprior_unc[pp,24] = par_data[plot[pp]]['mn_vpd_u'] 
        parprior[pp,25] = par_data[plot[pp]]['mx_vpd']  # max VPD threshold (GSI)
        parprior_unc[pp,25] = par_data[plot[pp]]['mx_vpd_u'] 
        parprior[pp,26] = par_data[plot[pp]]['mn_gain'] # minimum GPP benefit of increased LAI for labile allocation to be allowed
        parprior_unc[pp,26] = par_data[plot[pp]]['mn_gain_u']
        parprior[pp,27] = par_data[plot[pp]]['F_br']    # fraction of Cwood which is Cbranch
        parprior_unc[pp,27] = par_data[plot[pp]]['F_br_u']  
        parprior[pp,28] = par_data[plot[pp]]['F_croo']  # fraction of Cwood which is Ccoarseroot
        parprior_unc[pp,28] = par_data[plot[pp]]['F_croo_u']
        parprior[pp,29] = par_data[plot[pp]]['lab_rpl'] # labile replanting
        parprior_unc[pp,29] = par_data[plot[pp]]['lab_rpl_u']
        parprior[pp,30] = par_data[plot[pp]]['fol_rpl'] # foliar replanting
        parprior_unc[pp,30] = par_data[plot[pp]]['fol_rpl_u']
        parprior[pp,31] = par_data[plot[pp]]['roo_rpl'] # fine root replanting
        parprior_unc[pp,31] = par_data[plot[pp]]['roo_rpl_u']
        parprior[pp,32] = par_data[plot[pp]]['woo_rpl'] # wood replanting
        parprior_unc[pp,32] = par_data[plot[pp]]['woo_rpl_u']

        parprior[pp,37] = par_data[plot[pp]]['Ccwd_i']   # initial C stock for CWD
        parprior_unc[pp,37] = par_data[plot[pp]]['Ccwd_i_u']

    np.save(data_dir + project_par_npydata,[parprior,parprior_unc])
   
else:
    print "Loading parameters"
    parprior, parprior_unc = np.load(data_dir+project_par_npydata)

#-----------------------------------------------------------------------------
# Now deal with the observations 
obs = np.zeros((nosites,sim_length,34))-9999.  # check number of observations and their uncertainties

if project_obs_npydata not in os.listdir(data_dir):
    
    for pp in range(0,nosites):
        obs[pp,:,0] = obs_data[plot[pp]]['GPP']       # GPP
        obs[pp,:,1] = obs_data[plot[pp]]['LAI']       # LAI
        obs[pp,:,2] = obs_data[plot[pp]]['NEE']       # NEE
        obs[pp,:,3] = obs_data[plot[pp]]['woo']       # woody increment
        obs[pp,:,4] = obs_data[plot[pp]]['Reco']      # Reco
        obs[pp,:,5] = obs_data[plot[pp]]['Cfol']      # Cfol
        obs[pp,:,6] = obs_data[plot[pp]]['Cwoo']      # Cwood
        obs[pp,:,7] = obs_data[plot[pp]]['Croo']      # Croot
        obs[pp,:,8] = obs_data[plot[pp]]['Clit']      # Clit
        obs[pp,:,9] = obs_data[plot[pp]]['Csom']      # Csom
        obs[pp,:,10] = obs_data[plot[pp]]['Cagb']     # Cagb
        obs[pp,:,22] = obs_data[plot[pp]]['Cstem']    # Cstem
        obs[pp,:,24] = obs_data[plot[pp]]['Cbranch']  # Cbranch
        obs[pp,:,26] = obs_data[plot[pp]]['Ccroo']    # Ccoarseroot
        obs[pp,:,28] = obs_data[plot[pp]]['Cfol_max'] # maximum Cfol
        obs[pp,:,30] = obs_data[plot[pp]]['Evap']     # Evapotranspiration
        obs[pp,:,32] = obs_dataplot[pp]]['flit']      # Litter flux

        obs[pp,:,11] = obs_data[plot[pp]]['GPP_u']       # GPP
        obs[pp,:,12] = obs_data[plot[pp]]['LAI_u']       # LAI
        obs[pp,:,13] = obs_data[plot[pp]]['NEE_u']       # NEE
        obs[pp,:,14] = obs_data[plot[pp]]['woo_u']       # woody increment
        obs[pp,:,15] = obs_data[plot[pp]]['Reco_u']      # Reco
        obs[pp,:,16] = obs_data[plot[pp]]['Cfol_u']      # Cfol
        obs[pp,:,17] = obs_data[plot[pp]]['Cwoo_u']      # Cwood
        obs[pp,:,18] = obs_data[plot[pp]]['Croo_u']      # Croot
        obs[pp,:,19] = obs_data[plot[pp]]['Clit_u']      # Clit
        obs[pp,:,20] = obs_data[plot[pp]]['Csom_u']      # Csom
        obs[pp,:,21] = obs_data[plot[pp]]['Cagb_u']     # Cagb
        obs[pp,:,23] = obs_data[plot[pp]]['Cstem_u']    # Cstem
        obs[pp,:,25] = obs_data[plot[pp]]['Cbranch_u']  # Cbranch
        obs[pp,:,27] = obs_data[plot[pp]]['Ccroo_u']    # Ccoarseroot
        obs[pp,:,29] = obs_data[plot[pp]]['Cfol_max_u'] # maximum Cfol
        obs[pp,:,31] = obs_data[plot[pp]]['Evap_u']     # Evapotranspiration
        obs[pp,:,33] = obs_dataplot[pp]]['flit_u']      # Litter flux
                
    np.save(data_dir + project_obs_npydata,[obs,obs_unc])

else:
    print "Loading observations"
    obs = np.load(data_dir+project_obs_npydata)

prj=CARDAMOM(project_name="BALI_GEMplots_daily")
prj.setup(latitude,longitude,met,obs,parprior,parpriorunc,otherprior,otherpriorunc)
