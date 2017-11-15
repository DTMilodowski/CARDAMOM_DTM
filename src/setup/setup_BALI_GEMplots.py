# This script sets up the CARDAMOM simulations for the BALI GEM plot data.
# It is set up to run DALEC
# DTM 27/06/2017

import numpy as np
import datetime as dt
import os
import sys
sys.path.append('/exports/csce/datastore/geos/users/dmilodow/BALI/CARDAMOM_BALI/CARDAMOM_DTM/src')
import CARDAMOM_DTM as CAR
sys.path.append('/exports/csce/datastore/geos/users/dmilodow/BALI/CARDAMOM_BALI/CARDAMOM_DTM/src/setup')
import load_data as data

# Project data file (to load in if already generated
run = '002'
data_dir = "/home/dmilodow/DataStore_DTM/BALI/CARDAMOM_BALI/npydata/"
if data_dir + run not in os.listdir(data_dir):
    os.mkdir(data_dir + run)
project_met_npydata = "BALI_GEMplots_daily_drivers.npy"
project_par_npydata = "BALI_GEMplots_daily_params.npy"
project_obs_npydata = "BALI_GEMplots_daily_obs.npy"

# File list containing drivers, observations and plot coordinates
coordinate_file = "/exports/csce/datastore/geos/users/dmilodow/BALI/CARDAMOM_BALI/parameter_files/BALI_plot_coordinates.csv"
met_file = "/exports/csce/datastore/geos/users/dmilodow/BALI/CARDAMOM_BALI/met_data/BALI_metstation_daily_v1.csv"
# note that plot parameters may change each run, so use specified parameter file
par_file = "/exports/csce/datastore/geos/users/dmilodow/BALI/CARDAMOM_BALI/parameter_files/BALI_GEM_plot_params_"+run+".csv"

# first load in coordinates and other data
obs_dir = '/exports/csce/datastore/geos/users/dmilodow/BALI/CARDAMOM_BALI/CARDAMOM_DTM/src/setup/obs/'
obs_data = {}

plot, latitude, longitude = data.load_plot_coordinates(coordinate_file)
met_data = data.load_met_data(met_file)
par_data = data.load_par_data(par_file)
for i in range(0,len(plot)):
    obs_file = obs_dir+"CARDAMOM_obs_"+plot[i]+".csv"
    obs_data[plot[i]] = data.load_obs_data(obs_file)

# start date ### read from data file
dates = met_data['date'].astype('datetime64[D]')
d0 = dates[0]
DoY = (dates-dates.astype('datetime64[Y]')+1).astype('float')
tstep = met_data['tstep_days']
sim_length = tstep[-1]+1
nosites = len(plot)

#-----------------------------------------------------------------------------
# First deal with the met data 
met = np.zeros((nosites,sim_length,14))

if project_met_npydata not in os.listdir(data_dir+run):
    
    # First load the met data into the met array - use same driving data for now
    met[:,:,0] = tstep+1                                  #  0 = run day
    met[:,:,1] = met_data['mn2t']                         #  1 = minimum temperature oC
    met[:,:,2] = met_data['mx2t']                         #  2 = maximum temperature oC
    met[:,:,3] = met_data['ssrd']                         #  3 = surface shortwave radiation in MJ.m-2.day-1
    met[:,:,4] = 400.                                     #  4 = atmospheric CO2 concentration ppm
    met[:,:,5] = DoY.copy()                               #  5 = day of year
    met[:,:,6] = met_data['pptn']                         #  6 = lagged precipitation
    met[:,:,7] = -9999                                    #  7 = fire burned fraction - not applicable 
    met[:,:,8] = -9999                                    #  8 = deforestation fraction - not applicable 
    met[:,:,9] = met_data['mn2t_21d']                     #  9 = 21 day average min temperature oC
    met[:,:,11] = met_data['vpd_21d']                     # 11 = 21 day average vpd Pa
    met[:,:,12] = -9999                                   # 12 = forest management practice to accompany any clearing - not applicable
    met[:,:,13] = (met_data['mn2t']+met_data['mx2t'])/2.  # 13 = mean temperature oC ???
    
    # Next calculate photoperiod in s for each day in time series
    dayl = np.zeros(met.shape[:-1])
    dec = - np.arcsin( np.sin( 23.45 * np.pi/180. ) * np.cos( 2.0 * np.pi * ( DoY + 10.0 ) / 365.0 ) )
    ones = np.ones(sim_length)
    for ss in range(0,nosites):
        mult = latitude[ss]*np.pi/180.
        sinld = np.sin( mult ) * np.sin( dec )
        cosld = np.cos( mult ) * np.cos( dec )
        aob = np.max((-1.0*ones,np.min((ones,sinld / cosld), axis=0)),axis=0)
        met[ss,:,10] = 12.*60.*60.*(1.+2.*np.arcsin(aob)/np.pi) # 10 = photoperiod in s
    
    print '\tnodata test: ', np.sum(met==-9999,axis=0).sum(axis=0)
    np.save(data_dir+run+'/'+project_met_npydata,met)

else:
    print "Loading met drivers"
    met = np.load(data_dir+run+'/'+project_met_npydata)

#-----------------------------------------------------------------------------
# Now deal with the parameters 
parprior = np.zeros((nosites,100))-9999.
parprior_unc = np.zeros((nosites,100))-9999.
otherprior = np.zeros((nosites,100))-9999. # What are these for????
otherprior_unc = np.zeros((nosites,100))-9999.

if project_par_npydata not in os.listdir(data_dir+run):
    
    for pp in range(0,nosites):
        ii = par_data['plot']==plot[pp]
        parprior[pp,0] = par_data['m_r'][ii]      # Litter to SOM conversion rate 
        parprior_unc[pp,0] = par_data['m_r_u'][ii]
        parprior[pp,1] = par_data['f_a'][ii]       # Fraction of NPP respired
        parprior_unc[pp,1] = par_data['f_a_u'][ii]
        parprior[pp,2] = par_data['f_f'][ii]       # Fraction of NPP allocated to foliage 
        parprior_unc[pp,2] = par_data['f_f_u'][ii]
        parprior[pp,3] = par_data['f_r'][ii]       # Fraction of NPP allocated to roots
        parprior_unc[pp,3] = par_data['f_r_u'][ii]
        parprior[pp,4] = par_data['L_f'][ii]       # max leaf turnover (GSI) ! Leaf lifespan (CDEA)
        parprior_unc[pp,4] = par_data['L_f_u'][ii]
        parprior[pp,5] = par_data['t_w'][ii]       # Turnover rate of wood
        parprior_unc[pp,5] = par_data['t_w_u'][ii]
        parprior[pp,6] = par_data['t_r'][ii]       # Turnover rate of roots
        parprior_unc[pp,6] = par_data['t_r_u'][ii]
        parprior[pp,7] = par_data['t_l'][ii]       # Litter turnover rate
        parprior_unc[pp,7] = par_data['t_l_u'][ii]
        parprior[pp,8] = par_data['t_S'][ii]       # SOM turnover rate
        parprior_unc[pp,8] = par_data['t_S_u'][ii]
        parprior[pp,9] = par_data['theta'][ii]     # Parameter in exponential term of temperature
        parprior_unc[pp,9] = par_data['theta_u'][ii]
        parprior[pp,10] = par_data['C_eff'][ii]    # Canopy efficiency parameter (part of ACM) 
        parprior_unc[pp,10] = par_data['C_eff_u'][ii]
        parprior[pp,11] = par_data['B_day'][ii]    # max labile turnover(GSI) ! date of Clab release (CDEA)
        parprior_unc[pp,11] = par_data['B_day_u'][ii]
        parprior[pp,12] = par_data['f_l'][ii]      # Fraction allocated to Clab
        parprior_unc[pp,12] = par_data['f_l_u'][ii]
        parprior[pp,13] = par_data['R_l'][ii]      # min temp threshold (GSI) ! lab release duration period (CDEA)
        parprior_unc[pp,13] = par_data['R_l_u'][ii]
        parprior[pp,14] = par_data['F_day'][ii]    # max temp threshold (GSI)! date of leaf fall
        parprior_unc[pp,14] = par_data['F_day_u'][ii]
        parprior[pp,15] = par_data['mn_photo'][ii] # min photoperiod threshold (GSI)
        parprior_unc[pp,15] = par_data['mn_photo_u'][ii]
        parprior[pp,16] = par_data['LMA'][ii]      # LMA
        parprior_unc[pp,16] = par_data['LMA_u'][ii]
        parprior[pp,17] = par_data['Clab_i'][ii]   # initial C stock for labile
        parprior_unc[pp,17] = par_data['Clab_i_u'][ii]
        parprior[pp,18] = par_data['Cfol_i'][ii]   # initial C stock for foliar
        parprior_unc[pp,18] = par_data['Cfol_i_u'][ii]
        parprior[pp,19] = par_data['Croo_i'][ii]   # initial C stock for roots
        parprior_unc[pp,19] = par_data['Croo_i_u'][ii]
        parprior[pp,20] = par_data['Cwoo_i'][ii]   # initial C stock for litter
        parprior_unc[pp,20] = par_data['Cwoo_i_u'][ii]
        parprior[pp,21] = par_data['Clit_i'][ii]   # initial C stock for litter 
        parprior_unc[pp,21] = par_data['Clit_i_u'][ii]
        parprior[pp,22] = par_data['Csom_i'][ii]   # initial C stock for som
        parprior_unc[pp,22] = par_data['Csom_i_u'][ii]
        parprior[pp,23] = par_data['mx_photo'][ii] # max photoperiod threshold (GSI)
        parprior_unc[pp,23] = par_data['mx_photo_u'][ii]
        parprior[pp,24] = par_data['mn_vpd'][ii]   # min VPD threshold (GSI)
        parprior_unc[pp,24] = par_data['mn_vpd_u'][ii]
        parprior[pp,25] = par_data['mx_vpd'][ii]   # max VPD threshold (GSI)
        parprior_unc[pp,25] = par_data['mx_vpd_u'][ii]
        parprior[pp,26] = par_data['mn_gain'][ii]  # minimum GPP benefit of increased LAI for labile allocation to be allowed
        parprior_unc[pp,26] = par_data['mn_gain_u'][ii]
        parprior[pp,27] = par_data['F_br'][ii]     # fraction of Cwood which is Cbranch
        parprior_unc[pp,27] = par_data['F_br_u'][ii]
        parprior[pp,28] = par_data['F_croo'][ii]   # fraction of Cwood which is Ccoarseroot
        parprior_unc[pp,28] = par_data['F_croo_u'][ii]
        parprior[pp,29] = par_data['lab_rpl'][ii]  # labile replanting
        parprior_unc[pp,29] = par_data['lab_rpl_u'][ii]
        parprior[pp,30] = par_data['fol_rpl'][ii]  # foliar replanting
        parprior_unc[pp,30] = par_data['fol_rpl_u'][ii]
        parprior[pp,31] = par_data['roo_rpl'][ii]  # fine root replanting
        parprior_unc[pp,31] = par_data['roo_rpl_u'][ii]
        parprior[pp,32] = par_data['woo_rpl'][ii]  # wood replanting
        parprior_unc[pp,32] = par_data['woo_rpl_u'][ii]

        parprior[pp,37] = par_data['Ccwd_i'][ii]    # initial C stock for CWD
        parprior_unc[pp,37] = par_data['Ccwd_i_u'][ii]

    np.save(data_dir+run+'/' + project_par_npydata,[parprior,parprior_unc])
   
else:
    print "Loading parameters"
    parprior, parprior_unc = np.load(data_dir+run+'/'+project_par_npydata)

#-----------------------------------------------------------------------------
# Now deal with the observations 
obs = np.zeros((nosites,sim_length,34))-9999.  # check number of observations and their uncertainties

if project_obs_npydata not in os.listdir(data_dir+run):
    
    for pp in range(0,nosites):
        obs[pp,:,0] = obs_data[plot[pp]]['GPP']         # GPP
        obs[pp,:,1] = obs_data[plot[pp]]['LAI']         # LAI
        obs[pp,:,2] = obs_data[plot[pp]]['NEE']         # NEE
        obs[pp,:,3] = obs_data[plot[pp]]['woo']         # woody increment
        obs[pp,:,4] = obs_data[plot[pp]]['Reco']        # Reco
        obs[pp,:,5] = obs_data[plot[pp]]['Cfol']        # Cfol
        obs[pp,:,6] = obs_data[plot[pp]]['Cwoo']        # Cwood
        obs[pp,:,7] = obs_data[plot[pp]]['Croo']        # Croot
        obs[pp,:,8] = obs_data[plot[pp]]['Clit']        # Clit
        obs[pp,:,9] = obs_data[plot[pp]]['Csom']        # Csom
        obs[pp,:,10] = obs_data[plot[pp]]['Cagb']       # Cagb
        obs[pp,:,22] = obs_data[plot[pp]]['Cstem']      # Cstem
        obs[pp,:,24] = obs_data[plot[pp]]['Cbranch']    # Cbranch
        obs[pp,:,26] = obs_data[plot[pp]]['Ccroo']      # Ccoarseroot
        obs[pp,:,28] = obs_data[plot[pp]]['Cfol_max']   # maximum Cfol
        obs[pp,:,30] = obs_data[plot[pp]]['Evap']       # Evapotranspiration
        obs[pp,:,32] = obs_data[plot[pp]]['flit']       # Litter flux

        obs[pp,:,11] = obs_data[plot[pp]]['GPP_u']      # GPP
        obs[pp,:,12] = obs_data[plot[pp]]['LAI_u']      # LAI
        obs[pp,:,13] = obs_data[plot[pp]]['NEE_u']      # NEE
        obs[pp,:,14] = obs_data[plot[pp]]['woo_u']      # woody increment
        obs[pp,:,15] = obs_data[plot[pp]]['Reco_u']     # Reco
        obs[pp,:,16] = obs_data[plot[pp]]['Cfol_u']     # Cfol
        obs[pp,:,17] = obs_data[plot[pp]]['Cwoo_u']     # Cwood
        obs[pp,:,18] = obs_data[plot[pp]]['Croo_u']     # Croot
        obs[pp,:,19] = obs_data[plot[pp]]['Clit_u']     # Clit
        obs[pp,:,20] = obs_data[plot[pp]]['Csom_u']     # Csom
        obs[pp,:,21] = obs_data[plot[pp]]['Cagb_u']     # Cagb
        obs[pp,:,23] = obs_data[plot[pp]]['Cstem_u']    # Cstem
        obs[pp,:,25] = obs_data[plot[pp]]['Cbranch_u']  # Cbranch
        obs[pp,:,27] = obs_data[plot[pp]]['Ccroo_u']    # Ccoarseroot
        obs[pp,:,29] = obs_data[plot[pp]]['Cfol_max_u'] # maximum Cfol
        obs[pp,:,31] = obs_data[plot[pp]]['Evap_u']     # Evapotranspiration
        obs[pp,:,33] = obs_data[plot[pp]]['flit_u']     # Litter flux
                
    np.save(data_dir+run+'/' + project_obs_npydata,obs)

else:
    print "Loading observations"
    obs = np.load(data_dir+run+'/'+project_obs_npydata)

prj=CAR.CARDAMOM(project_name="BALI_GEMplots_daily")
prj.setup(latitude,longitude,met,obs,parprior,parprior_unc,otherprior,otherprior_unc)
