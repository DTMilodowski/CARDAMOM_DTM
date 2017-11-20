# Some scripts to load in the various data sets

import numpy as np

def load_plot_coordinates(coordinate_file):
    dtypes = {'names': ('plot','lat','long'),'formats':('S8', 'float32', 'float32')}
    info = np.genfromtxt(coordinate_file,skiprows=1,delimiter=',',dtype=dtypes)
    return info['plot'],info['lat'],info['long']

def load_met_data(met_file):
    met_dtype = {'names': ('tstep_days','date','mn2t','mx2t','vpd','ssrd','pptn','mn2t_21d','mx2t_21d','vpd_21d','ssrd_21d','pptn_21d'),'formats':('int32', 'S10', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32')}
    met_data = np.genfromtxt(met_file,skiprows=1,delimiter=',',dtype=met_dtype)
    return met_data
    
def load_obs_data(obs_file):
    obs_dtype = {'names': ('tstep_days','date','GPP','GPP_u','LAI','LAI_u','NEE','NEE_u','woo','woo_u','Reco','Reco_u','Cfol','Cfol_u','Cwoo','Cwoo_u','Croo','Croo_u','Clit','Clit_u','Csom','Csom_u','Cagb','Cagb_u','Cstem','Cstem_u','Cbranch','Cbranch_u','Ccroo','Ccroo_u','Cfol_max','Cfol_max_u','Evap','Evap_u','flit','flit_u','flit_acc_days'),'formats':('int','S10', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float','int')}
    obs_data = np.genfromtxt(obs_file,skiprows=1,delimiter=',',dtype=obs_dtype)
    return obs_data

def load_par_data(par_file):
    par_dtype = {'names': ('plot','m_r','m_r_u','f_a','f_a_u','f_f','f_f_u','f_r','f_r_u','L_f','L_f_u','t_w','t_w_u','t_r','t_r_u','t_l','t_l_u','t_S','t_S_u','theta','theta_u','C_eff','C_eff_u','B_day','B_day_u','f_l','f_l_u','R_l','R_l_u','F_day','F_day_u','mn_photo','mn_photo_u','LMA','LMA_u','Clab_i','Clab_i_u','Cfol_i','Cfol_i_u','Croo_i','Croo_i_u','Cwoo_i','Cwoo_i_u','Clit_i','Clit_i_u','Csom_i','Csom_i_u','mx_photo','mx_photo_u','mn_vpd','mn_vpd_u','mx_vpd','mx_vpd_u','mn_gain','mn_gain_u','F_br','F_br_u','F_croo','F_croo_u','lab_rpl','lab_rpl_u','fol_rpl','fol_rpl_u','roo_rpl','roo_rpl_u','woo_rpl','woo_rpl_u','Ccwd_i','Ccwd_i_u'), 'formats':('S8', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32')}
    par_data = np.genfromtxt(par_file,skiprows=1,delimiter=',',dtype=par_dtype)   
    return par_data
