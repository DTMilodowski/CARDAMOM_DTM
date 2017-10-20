#===============================================================================
# plot_CARDAMOM_output.py
#-------------------------------------------------------------------------------
# This script contains functions to produce standard plots of CARDAMOM output,
# such as carbon pools through time in addition to data where present
#-------------------------------------------------------------------------------
# Author: D.T.Milodowski
# Date: October 2017
#===============================================================================
import numpy as np
from matplotlib import pyplot as plt

#---------------------
# plot_carbon_pools_ts
# --------------------
# Plot a time series of C stocks, with observations if available.
# The function reads in two dictionaries:
# - the model output (median, ulim, llim)
# - the obserations where available
# Note that these dictionaries should have uniform naming structures for the
# carbon pools (woo,fol,lab,roo,lit,cwd,som)
# Also takes optional arguments for start and end timestep. initially these will
# be index references (i.e. model timestep) but this will ultimately be altered
# to give options to specify date ranges.
def plot_carbon_pools_ts(model,obs,start_tstep=0,end_tstep=-1):

    return 0
