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
colour = ['#46E900','#1A2BCE','#E0007F']

#---------------------
# plot_carbon_pools_ts
# --------------------
# Plot a time series of C stocks, with observations if available.
# The function reads in two dictionaries:
# - the model output (median, ulim, llim)
# - the obserations where available
# Note that these dictionaries should have uniform naming structures for the
# carbon pools (Cwoo,Cfol,Clab,Croo,Clit,Ccwd,Csom)
# Also takes optional arguments for start and end timestep. initially these will
# be index references (i.e. model timestep) but this will ultimately be altered
# to give options to specify date ranges.
def plot_carbon_pools_ts(model,obs,start_tstep=False,end_tstep=False):
    fig = plt.figure(1, facecolor='White',figsize=[8,14])

    # Plot a -> Cwood
    ax1a = plt.subplot2grid((7,1),(0,0))
    ax1a.annotate('a - C$_{wood}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1a.set_ylabel('C$_{wood}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    ax1a.fill_between(model['time'],model['Cwoo'][:,3],model['Cwoo'][:,4],color=colour[0],alpha=0.2)
    ax1a.plot(model['time'],model['Cwoo'][:,1],'-',color=colour[0])

    if 'Cwoo' in obs.keys(): # check for observations
        if 'Cwoo_u' in obs.keys(): # check for uncertainty bounds
            ax1a.error_bar(obs['time'],obs['Cwoo'],yerr=obs['Cwoo_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1a.plot(obs['time'],obs['Cwoo'],marker='o',c='black',mec='black',mfc='black')

    
    # Plot b -> Cfol
    ax1b = plt.subplot2grid((7,1),(1,0),sharex=ax1a)
    ax1b.annotate('b - C$_{foliar}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1b.set_ylabel('C$_{fol}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    ax1b.fill_between(model['time'],model['Cfol'][:,3],model['Cfol'][:,4],color=colour[0],alpha=0.2)
    ax1b.plot(model['time'],model['Cfol'][:,1],'-',color=colour[0])
    
    if 'Cfol' in obs.keys(): # check for observations
        if 'Cfol_u' in obs.keys(): # check for uncertainty bounds
            ax1b.error_bar(obs['time'],obs['Cfol'],yerr=obs['Cfol_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1b.plot(obs['time'],obs['Cfol'],marker='o',c='black',mec='black',mfc='black')
    
    # Plot c -> Croot
    ax1c = plt.subplot2grid((7,1),(2,0),sharex=ax1a)
    ax1c.annotate('c - C$_{root}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1c.set_ylabel('C$_{root}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    ax1c.fill_between(model['time'],model['Croo'][:,3],model['Croo'][:,4],color=colour[0],alpha=0.2)
    ax1c.plot(model['time'],model['Croo'][:,1],'-',color=colour[0])
    
    if 'Croo' in obs.keys(): # check for observations
        if 'Cwro_u' in obs.keys(): # check for uncertainty bounds
            ax1c.error_bar(obs['time'],obs['Croo'],yerr=obs['Croo_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1c.plot(obs['time'],obs['Croo'],marker='o',c='black',mec='black',mfc='black')
    
    # Plot d -> Clab
    ax1d = plt.subplot2grid((7,1),(3,0),sharex=ax1a)
    ax1d.annotate('d - C$_{labile}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1d.set_ylabel('C$_{lab}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    ax1d.fill_between(model['time'],model['Clab'][:,3],model['Clab'][:,4],color=colour[0],alpha=0.2)
    ax1d.plot(model['time'],model['Clab'][:,1],'-',color=colour[0])

    if 'Clab' in obs.keys(): # check for observations
        if 'Clab_u' in obs.keys(): # check for uncertainty bounds
            ax1d.error_bar(obs['time'],obs['Clab'],yerr=obs['Clab_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1d.plot(obs['time'],obs['Clab'],marker='o',c='black',mec='black',mfc='black')
            
    # Plot e -> Clit
    ax1e = plt.subplot2grid((7,1),(4,0),sharex=ax1a)
    ax1e.annotate('e - C$_{litter}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1e.set_ylabel('C$_{lit}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    ax1e.fill_between(model['time'],model['Clit'][:,3],model['Clit'][:,4],color=colour[1],alpha=0.2)
    ax1e.plot(model['time'],model['Clit'][:,1],'-',color=colour[1])
    
    if 'Clit' in obs.keys(): # check for observations
        if 'Clit_u' in obs.keys(): # check for uncertainty bounds
            ax1e.error_bar(obs['time'],obs['Clit'],yerr=obs['Clit_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1e.plot(obs['time'],obs['Clit'],marker='o',c='black',mec='black',mfc='black')
    
    # Plot f -> Ccwd
    ax1f = plt.subplot2grid((7,1),(5,0),sharex=ax1a)
    ax1f.annotate('f - C$_{cwd}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1f.set_ylabel('C$_{cwd}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    ax1f.fill_between(model['time'],model['Ccwd'][:,3],model['Ccwd'][:,4],color=colour[1],alpha=0.2)
    ax1f.plot(model['time'],model['Ccwd'][:,1],'-',color=colour[1])

    if 'Ccwd' in obs.keys(): # check for observations
        if 'Ccwd_u' in obs.keys(): # check for uncertainty bounds
            ax1f.error_bar(obs['time'],obs['Ccwd'],yerr=obs['Ccwd_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1f.plot(obs['time'],obs['Ccwd'],marker='o',c='black',mec='black',mfc='black')    
    
    # Plot g -> Csom
    ax1g = plt.subplot2grid((7,1),(6,0),sharex=ax1a)
    ax1g.annotate('g - C$_{SOM}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1g.set_ylabel('C$_{som}$ / g(C) m$^{-2}$',fontsize = axis_size)
    ax1g.set_xlabel('timestep',fontsize = axis_size)
    
    ax1g.fill_between(model['time'],model['Csom'][:,3],model['Csom'][:,4],color=colour[1],alpha=0.2)
    ax1g.plot(model['time'],model['Csom'][:,1],'-',color=colour[1])
    
    if 'Csom' in obs.keys(): # check for observations
        if 'Csom_u' in obs.keys(): # check for uncertainty bounds
            ax1g.error_bar(obs['time'],obs['Csom'],yerr=obs['Csom_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1g.plot(obs['time'],obs['Csom'],marker='o',c='black',mec='black',mfc='black')
    

    # set xlimits if desired
    if start_tstep!=False:
        ax1a.set_xlim(xmin=start_tstep)
    if end_tstep!=False:
        ax1a.set_xlim(xmax=end_tstep)
    else:
        ax1a.set_xlim(xmax=model['time'].size)

    plt.tight_layout()
    plt.show()
    return 0


#----------------------
# plot_carbon_fluxes_ts
# ---------------------
# Plot a time series of C fluxes into and out of ecosystem, with observations if
# available. LAI time series also plotted to provide context
# The function reads in two dictionaries:
# - the model output (median, ulim, llim)
# - the obserations where available
# Note that these dictionaries should have uniform naming structures for the
# carbon fluxes
# - lai : leaf area index
# - gpp : Gross Primary Production
# - nee : Net ecosystem exchange
# - Rauto : autotrophic respiration
# - Rhet : heterotrophic respiration
# Also takes optional arguments for start and end timestep. initially these will
# be index references (i.e. model timestep) but this will ultimately be altered
# to give options to specify date ranges.

#----------------------
# plot_litter_components_ts
# ---------------------
# Plot a time series of the model output relating specifically to the litter
# component of DALEC
# The function reads in two dictionaries:
# - the model output (median, ulim, llim)
# - the obserations where available
# Note that these dictionaries should have uniform naming structures
# - lai : leaf area index
# - gsi : growth season index
# - Clit : litter stocks
# - Croo : root carbon stocks
# - Ccwd : coarse woody debris stocks
# - flux_fol_lit : carbon flux from foliar pool to litter (litterfall)
# - flux_cwd_lit : carbon flux from CWD pool to litter
# - flux_root_lit : carbon flux from root pool to litter
# Also takes optional arguments for start and end timestep. initially these will
# be index references (i.e. model timestep) but this will ultimately be altered
# to give options to specify date ranges.
def plot_litter_components_ts(model,obs,start_tstep=False,end_tstep=False):
    fig = plt.figure(3, facecolor='White',figsize=[8,14])

    # Plot a -> LAI
    ax1a = plt.subplot2grid((8,1),(0,0))
    ax1a.annotate('a - LAI$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1a.set_ylabel('LAI / m$^2$m$^{-2}$',fontsize = axis_size)
    
    ax1a.fill_between(model['time'],model['lai'][:,3],model['lai'][:,4],color=colour[0],alpha=0.2)
    ax1a.plot(model['time'],model['lai'][:,1],'-',color=colour[0])

    if 'lai' in obs.keys(): # check for observations
        if 'lai' in obs.keys(): # check for uncertainty bounds
            ax1a.error_bar(obs['time'],obs['lai'],yerr=obs['lai_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1a.plot(obs['time'],obs['lai'],marker='o',c='black',mec='black',mfc='black')

    
    # Plot b -> gsi
    ax1b = plt.subplot2grid((8,1),(1,0),sharex=ax1a)
    ax1b.annotate('b - Growth Season Index (GSI)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1b.set_ylabel('GSI',fontsize = axis_size)
    
    ax1b.fill_between(model['time'],model['gsi'][:,3],model['gsi'][:,4],color=colour[0],alpha=0.2)
    ax1b.plot(model['time'],model['gsi'][:,1],'-',color=colour[0])
    
    
    # Plot c -> Clit
    ax1c = plt.subplot2grid((8,1),(2,0),sharex=ax1a)
    ax1c.annotate('e - C$_{litter}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1c.set_ylabel('C$_{lit}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    ax1c.fill_between(model['time'],model['Clit'][:,3],model['Clit'][:,4],color=colour[1],alpha=0.2)
    ax1c.plot(model['time'],model['Clit'][:,1],'-',color=colour[1])
    
    if 'Clit' in obs.keys(): # check for observations
        if 'Clit_u' in obs.keys(): # check for uncertainty bounds
            ax1c.error_bar(obs['time'],obs['Clit'],yerr=obs['Clit_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1c.plot(obs['time'],obs['Clit'],marker='o',c='black',mec='black',mfc='black')

            
    # Plot d -> Croot
    ax1d = plt.subplot2grid((8,1),(3,0),sharex=ax1a)
    ax1d.annotate('c - C$_{root}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1d.set_ylabel('C$_{root}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    ax1d.fill_between(model['time'],model['Croo'][:,3],model['Croo'][:,4],color=colour[0],alpha=0.2)
    ax1d.plot(model['time'],model['Croo'][:,1],'-',color=colour[0])
    
    if 'Croo' in obs.keys(): # check for observations
        if 'Cwro_u' in obs.keys(): # check for uncertainty bounds
            ax1d.error_bar(obs['time'],obs['Croo'],yerr=obs['Croo_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1d.plot(obs['time'],obs['Croo'],marker='o',c='black',mec='black',mfc='black')
    
    # Plot e -> Ccwd
    ax1e = plt.subplot2grid((8,1),(4,0),sharex=ax1a)
    ax1e.annotate('d - C$_{cwd}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1e.set_ylabel('C$_{cwd}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    ax1e.fill_between(model['time'],model['Ccwd'][:,3],model['Ccwd'][:,4],color=colour[0],alpha=0.2)
    ax1e.plot(model['time'],model['Ccwd'][:,1],'-',color=colour[0])

    if 'Ccwd' in obs.keys(): # check for observations
        if 'Ccwd_u' in obs.keys(): # check for uncertainty bounds
            ax1e.error_bar(obs['time'],obs['Ccwd'],yerr=obs['Ccwd_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1e.plot(obs['time'],obs['Ccwd'],marker='o',c='black',mec='black',mfc='black')
            
    # Plot f -> Cfol into Clit
    ax1f = plt.subplot2grid((8,1),(5,0),sharex=ax1a)
    ax1f.annotate('f - litterfall', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1f.set_ylabel('litter flux from foliage / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    
    ax1f.fill_between(model['time'],model['flux_fol_lit'][:,3],model['flux_fol_lit'][:,4],color=colour[1],alpha=0.2)
    ax1f.plot(model['time'],model['flux_fol_lit'][:,1],'-',color=colour[1])
    
    if 'flux_fol_lit' in obs.keys(): # check for observations
        if 'flux_fol_lit_u' in obs.keys(): # check for uncertainty bounds
            ax1f.error_bar(obs['time'],obs['flux_fol_lit'],yerr=obs['flux_fol_lit_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1f.plot(obs['time'],obs['flux_fol_lit'],marker='o',c='black',mec='black',mfc='black')

    # Plot g -> cwd fluxes into litter pool
    ax1g = plt.subplot2grid((8,1),(6,0),sharex=ax1a)
    ax1g.annotate('f - litter flux from CWD', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1g.set_ylabel('litter flux from CWD / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    
    ax1g.fill_between(model['time'],model['flux_cwd_lit'][:,3],model['flux_cwd_lit'][:,4],color=colour[1],alpha=0.2)
    ax1g.plot(model['time'],model['flux_cwd_lit'][:,1],'-',color=colour[1])


    # Plot h -> root fluxes into litter pool
    ax1h = plt.subplot2grid((8,1),(6,0),sharex=ax1a)
    ax1h.annotate('f - litter flux from CWD', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1h.set_ylabel('litter flux from CWD / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    
    ax1h.fill_between(model['time'],model['flux_root_lit'][:,3],model['flux_root_lit'][:,4],color=colour[1],alpha=0.2)
    ax1h.plot(model['time'],model['flux_root_lit'][:,1],'-',color=colour[1])
    ax1h = plt.subplot2grid((8,1),(7,0),sharex=ax1a)
    ax1h.annotate('f - litter loss fluxes', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1h.set_ylabel('litter flux from roots / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    
    ax1h.fill_between(model['time'],model['flux_root_lit'][:,3],model['flux_root_lit'][:,4],color=colour[1],alpha=0.2)
    ax1h.plot(model['time'],model['flux_root_lit'][:,1],'-',color=colour[1])

    # set xlimits if desired
    if start_tstep!=False:
        ax1a.set_xlim(xmin=start_tstep)
    if end_tstep!=False:
        ax1a.set_xlim(xmax=end_tstep)
    else:
        ax1a.set_xlim(xmax=model['time'].size)

    plt.tight_layout()
    plt.show()
