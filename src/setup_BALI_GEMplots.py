# This script sets up the CARDAMOM simulations for the BALI GEM plot data.
# It is set up to run DALEC
# DTM 27/06/2017

import numpy as np
import datetime as dt
import os
from CARDAMOM_DTM import *
import load_data as data

# Project data file (to load in if already generated
data_dir = ""
project_met_npydata = "BALI_GEMplots_daily_drivers.npy"
project_par_npydata = "BALI_GEMplots_daily_params.npy"
project_obs_npydata = "BALI_GEMplots_daily_obs.npy"

# File list containing drivers, observations and plot coordinates
coordinate_file = "BALI_plot_coordinates.csv"
met_file = "<SOMETHING>.csv"
par_file = "<SOMETHING>.csv"
obs_file = "<SOMETHING>.csv"

# first load in coordinates
plot, latitude, longitude = data.load_plot_coordinates(coordinate_file)
dates_met,met_data = data.load_met_data(met_file)
par_data = data.load_obs_data(par_file)
dates_obs,obs_data = data.load_obs_data(obs_file)

# start date ### read from data file
d0 = dates_met[0]
sim_length = dates_met[-1]-dates_met[0]
DoY = dates_met - dates_met.astype('timedelta64[D]')
tstep = dates_met - dates_met[0] + np.timedelta[1,'D']
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
parprior = np.zeros((nosites,50))-9999.
parprior_unc = np.zeros((nosites,50))-9999.

if project_par_npydata not in os.listdir(data_dir):
    
    for pp in range(0,nosites):
        parprior[pp,0] = par_data[plot[pp]]['m_r']      # Litter to SOM conversion rate 
        parprior[pp,1] = par_data[plot[pp]]['f_a']      # Fraction of NPP respired
        parprior[pp,2] = par_data[plot[pp]]['f_f']      # Fraction of NPP allocated to foliage
        parprior[pp,3] = par_data[plot[pp]]['f_r']      # Fraction of NPP allocated to roots
        parprior[pp,4] = par_data[plot[pp]]['L_f']      # max leaf turnover (GSI) ! Leaf lifespan (CDEA)
        parprior[pp,5] = par_data[plot[pp]]['t_w']      # Turnover rate of wood
        parprior[pp,6] = par_data[plot[pp]]['t_r']      # Turnover rate of roots
        parprior[pp,7] = par_data[plot[pp]]['t_l']      # Litter turnover rate
        parprior[pp,8] = par_data[plot[pp]]['t_S']      # SOM turnover rate
        parprior[pp,9] = par_data[plot[pp]]['theta']    # Parameter in exponential term of temperature
        parprior[pp,10] = par_data[plot[pp]]['C_eff']   # Canopy efficiency parameter (part of ACM)
        parprior[pp,11] = par_data[plot[pp]]['B_day']   # max labile turnover(GSI) ! date of Clab release (CDEA)
        parprior[pp,12] = par_data[plot[pp]]['f_l']     # Fraction allocated to Clab
        parprior[pp,13] = par_data[plot[pp]]['R_l']     # min temp threshold (GSI) ! lab release duration period (CDEA)
        parprior[pp,14] = par_data[plot[pp]]['F_day']   # max temp threshold (GSI)! date of leaf fall
        parprior[pp,15] = par_data[plot[pp]]['mn_photo']# min photoperiod threshold (GSI)
        parprior[pp,16] = par_data[plot[pp]]['LMA']     # LMA
        parprior[pp,17] = par_data[plot[pp]]['Clab_i']  # initial C stock for labile
        parprior[pp,18] = par_data[plot[pp]]['Cfol_i']  # initial C stock for foliar
        parprior[pp,19] = par_data[plot[pp]]['Croo_i']  # initial C stock for roots
        parprior[pp,20] = par_data[plot[pp]]['Cwoo_i']  # initial C stock for litter
        parprior[pp,21] = par_data[plot[pp]]['Clit_i']  # initial C stock for litter
        parprior[pp,22] = par_data[plot[pp]]['Csom_i']  # initial C stock for som
        parprior[pp,23] = par_data[plot[pp]]['mx_photo']# max photoperiod threshold (GSI)
        parprior[pp,24] = par_data[plot[pp]]['mn_vpd']  # min VPD threshold (GSI)
        parprior[pp,25] = par_data[plot[pp]]['mx_vpd']  # max VPD threshold (GSI)
        parprior[pp,26] = par_data[plot[pp]]['mn_gain'] # minimum GPP benefit of increased LAI for labile allocation to be allowed
        parprior[pp,27] = par_data[plot[pp]]['F_br']    # fraction of Cwood which is Cbranch
        parprior[pp,28] = par_data[plot[pp]]['F_croo']  # fraction of Cwood which is Ccoarseroot
        parprior[pp,29] = par_data[plot[pp]]['lab_rpl'] # labile replanting
        parprior[pp,30] = par_data[plot[pp]]['fol_rpl'] # foliar replanting
        parprior[pp,31] = par_data[plot[pp]]['roo_rpl'] # fine root replanting
        parprior[pp,32] = par_data[plot[pp]]['woo_rpl'] # wood replanting
        #parprior[pp,33] = par_data[plot[pp]]['']
        #parprior[pp,34] = par_data[plot[pp]]['']
        #parprior[pp,35] = par_data[plot[pp]]['']
        #parprior[pp,36] = par_data[plot[pp]]['']
        parprior[pp,37] = par_data[plot[pp]]['Ccwd_i']   # initial C stock for CWD
        #parprior[pp,38] = par_data[plot[pp]]['']
        #parprior[pp,39] = par_data[plot[pp]]['']
        #parprior[pp,40] = par_data[plot[pp]]['']
        #parprior[pp,41] = par_data[plot[pp]]['']
        #parprior[pp,42] = par_data[plot[pp]]['']
        #parprior[pp,43] = par_data[plot[pp]]['']
        #parprior[pp,45] = par_data[plot[pp]]['']
        #parprior[pp,46] = par_data[plot[pp]]['']
        #parprior[pp,47] = par_data[plot[pp]]['']
        #parprior[pp,48] = par_data[plot[pp]]['']
        #parprior[pp,49] = par_data[plot[pp]]['']

    np.save(data_dir + project_par_npydata,pars)


#-----------------------------------------------------------------------------
# Now deal with the observations 

   
else:
    print "Loading parameters"
    met = np.load(data_dir+project_par_npydata)
