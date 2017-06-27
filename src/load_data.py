# Some scripts to load in the various data sets

import numpy as np

def load_met_data(met_file):
    met_dtype = {'names': ('tstep_days','date','mn2t','mx2t','vpd','ssrd','pptn','mn2t_21d','mx2t_21d','vpd_21d','ssrd_21d','pptn_21d'),'formats':('int', 'datetime64[D]', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float')}
    met_data = np.genfromtxt(met_file,skiprows=1,delimiter=',',dtype=met_dtype)
    return met_data
    
