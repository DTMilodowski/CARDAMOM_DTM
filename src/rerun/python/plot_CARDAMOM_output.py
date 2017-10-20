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
from matplotlib import rcParams
# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2

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
def plot_carbon_pools_ts(model,obs,start_tstep=False,end_tstep=False):
    fig = plt.figure(1, facecolor='White',figsize=[8,14])

    # Plot a -> Cwood
    ax1a = plt.subplot2grid((7,1),(0,0))
    ax1a.annotate('a - C$_{wood}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1a.set_ylabel('C$_{wood}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    # Plot b -> Cfol
    ax1b = plt.subplot2grid((7,1),(1,0),sharex=ax1a)
    ax1b.annotate('b - C$_{foliar}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1b.set_ylabel('C$_{fol}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    # Plot c -> Croot
    ax1c = plt.subplot2grid((7,1),(2,0),sharex=ax1a)
    ax1c.annotate('c - C$_{root}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1c.set_ylabel('C$_{root}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    # Plot d -> Clab
    ax1d = plt.subplot2grid((7,1),(3,0),sharex=ax1a)
    ax1d.annotate('d - C$_{labile}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1d.set_ylabel('C$_{lab}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    # Plot e -> Clit
    ax1e = plt.subplot2grid((7,1),(4,0),sharex=ax1a)
    ax1e.annotate('e - C$_{litter}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1e.set_ylabel('C$_{lit}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    # Plot f -> Ccwd
    ax1f = plt.subplot2grid((7,1),(5,0),sharex=ax1a)
    ax1f.annotate('f - C$_{cwd}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1f.set_ylabel('C$_{cwd}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    # Plot g -> Csom
    ax1g = plt.subplot2grid((7,1),(6,0),sharex=ax1a)
    ax1g.annotate('g - C$_{SOM}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1g.set_ylabel('C$_{som}$ / g(C) m$^{-2}$',fontsize = axis_size)
    ax1g.set_xlabel('timestep',fontsize = axis_size)

    # set xlimits if desired
    if start_tstep!=False:
        ax1a.set_xlim(xmin=start_tstep)
    if end_tstep!=False:
        ax1a.set_xlim(xmax=end_tstep)

    plt.tight_layout()
    plt.show()
    return 0
