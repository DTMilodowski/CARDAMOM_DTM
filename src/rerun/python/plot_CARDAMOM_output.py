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
from scipy import stats
import datetime as dt
import pandas as pd
import seaborn as sns
sns.set()

# Set up some basic parameters for the plots
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
def plot_carbon_pools_ts(model,obs,start_tstep=False,end_tstep=False,figname=''):
    fig = plt.figure(1, facecolor='White',figsize=[8,14])

    # Plot a -> Cwood
    ax1a = plt.subplot2grid((7,1),(0,0))
    ax1a.annotate('a - C$_{wood}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1a.set_ylabel('C$_{wood}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    ax1a.fill_between(model['time'],model['Cwoo'][:,3],model['Cwoo'][:,4],color=colour[0],alpha=0.2)
    ax1a.plot(model['time'],model['Cwoo'][:,1],'-',color=colour[0])

    if 'Cwoo' in obs.keys(): # check for observations
        if 'Cwoo_u' in obs.keys(): # check for uncertainty bounds
            ax1a.errorbar(obs['time'],obs['Cwoo'],yerr=obs['Cwoo_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black',elinewidth=0.1)
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
            ax1b.errorbar(obs['time'],obs['Cfol'],yerr=obs['Cfol_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1b.plot(obs['time'],obs['Cfol'],marker='o',c='black',mec='black',mfc='black')
    
    # Plot c -> Croot
    ax1c = plt.subplot2grid((7,1),(2,0),sharex=ax1a)
    ax1c.annotate('c - C$_{root}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1c.set_ylabel('C$_{root}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    ax1c.fill_between(model['time'],model['Croo'][:,3],model['Croo'][:,4],color=colour[0],alpha=0.2)
    ax1c.plot(model['time'],model['Croo'][:,1],'-',color=colour[0])
    
    if 'Croo' in obs.keys(): # check for observations
        if 'Croo_u' in obs.keys(): # check for uncertainty bounds
            ax1c.errorbar(obs['time'],obs['Croo'],yerr=obs['Croo_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
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
            ax1d.errorbar(obs['time'],obs['Clab'],yerr=obs['Clab_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
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
            ax1e.errorbar(obs['time'],obs['Clit'],yerr=obs['Clit_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
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
            ax1f.errorbar(obs['time'],obs['Ccwd'],yerr=obs['Ccwd_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
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
            ax1g.errorbar(obs['time'],obs['Csom'],yerr=obs['Csom_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1g.plot(obs['time'],obs['Csom'],marker='o',c='black',mec='black',mfc='black')
    

    # set xlimits if desired
    if start_tstep!=False:
        ax1a.set_xlim(xmin=start_tstep)
    if end_tstep!=False:
        ax1a.set_xlim(xmax=end_tstep)
    else:
        ax1a.set_xlim(xmax=model['time'].size)

    ax1a.set_ylim(0,40000)
    ax1b.set_ylim(0,1000)
    ax1c.set_ylim(0,1000)
    ax1d.set_ylim(0,1000)
    ax1e.set_ylim(0,1000)
    ax1f.set_ylim(0,3000)
    ax1g.set_ylim(0,40000)
        
    plt.tight_layout()
    if len(figname)>0:
        plt.savefig(figname)
    #plt.show()
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
# - gpp : Gross Primary Production
# - npp : Net Primary Production
# - nee : Net ecosystem exchange
# Also takes optional arguments for start and end timestep. initially these will
# be index references (i.e. model timestep) but this will ultimately be altered
# to give options to specify date ranges.
def plot_carbon_fluxes_ts(model,obs,start_tstep=False,end_tstep=False,figname=''):
    fig = plt.figure(2, facecolor='White',figsize=[8,6])

    # Plot a -> GPP
    ax1a = plt.subplot2grid((3,1),(0,0))
    ax1a.annotate('a - Gross Primary Productivity', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1a.set_ylabel('GPP / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    
    ax1a.fill_between(model['time'],model['gpp'][:,3],model['gpp'][:,4],color=colour[0],alpha=0.2)
    ax1a.plot(model['time'],model['gpp'][:,1],'-',color=colour[0])

    if 'gpp' in obs.keys(): # check for observations
        if 'gpp_u' in obs.keys(): # check for uncertainty bounds
            ax1a.errorbar(obs['time'],obs['gpp'],yerr=obs['gpp_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1a.plot(obs['time'],obs['gpp'],marker='o',c='black',mec='black',mfc='black')

    
    # Plot b -> NPP
    ax1b = plt.subplot2grid((3,1),(1,0),sharex=ax1a)
    ax1b.annotate('b - Net Primary Productivity', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1b.set_ylabel('NPP / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    
    ax1b.fill_between(model['time'],model['npp'][:,3],model['npp'][:,4],color=colour[0],alpha=0.2)
    ax1b.plot(model['time'],model['npp'][:,1],'-',color=colour[0])
    
    if 'npp' in obs.keys(): # check for observations
        if 'npp_u' in obs.keys(): # check for uncertainty bounds
            ax1b.errorbar(obs['time'],obs['npp'],yerr=obs['npp_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1b.plot(obs['time'],obs['npp'],marker='o',c='black',mec='black',mfc='black')
    
    # Plot c -> NEE
    ax1c = plt.subplot2grid((3,1),(2,0),sharex=ax1a)
    ax1c.annotate('c - Net Ecosystem Exchange', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1c.set_ylabel('NEE / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    
    ax1c.fill_between(model['time'],model['nee'][:,3],model['nee'][:,4],color=colour[0],alpha=0.2)
    ax1c.plot(model['time'],model['nee'][:,1],'-',color=colour[0])
    
    if 'nee' in obs.keys(): # check for observations
        if 'nee_u' in obs.keys(): # check for uncertainty bounds
            ax1c.errorbar(obs['time'],obs['nee'],yerr=obs['nee_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1c.plot(obs['time'],obs['nee'],marker='o',c='black',mec='black',mfc='black')
    
    # set xlimits if desired
    if start_tstep!=False:
        ax1a.set_xlim(xmin=start_tstep)
    if end_tstep!=False:
        ax1a.set_xlim(xmax=end_tstep)
    else:
        ax1a.set_xlim(xmax=model['time'].size)

    ax1a.set_ylim(0,20)
    ax1b.set_ylim(-10,1)
    ax1c.set_ylim(-7,7)
        
    plt.tight_layout()
    if len(figname)>0:
        plt.savefig(figname)
    #plt.show()
    return 0


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
def plot_litter_components_ts(model,obs,start_tstep=False,end_tstep=False,figname=''):
    fig = plt.figure(3, facecolor='White',figsize=[8,14])

    # Plot a -> LAI
    ax1a = plt.subplot2grid((9,1),(0,0))
    ax1a.annotate('a - LAI', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1a.set_ylabel('LAI / m$^2$m$^{-2}$',fontsize = axis_size)
    
    ax1a.fill_between(model['time'],model['lai'][:,3],model['lai'][:,4],color=colour[0],alpha=0.2)
    ax1a.plot(model['time'],model['lai'][:,1],'-',color=colour[0])

    if 'lai' in obs.keys(): # check for observations
        if 'lai' in obs.keys(): # check for uncertainty bounds
            ax1a.errorbar(obs['time'],obs['lai'],yerr=obs['lai_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1a.plot(obs['time'],obs['lai'],marker='o',c='black',mec='black',mfc='black')

    
    # Plot b -> gsi
    ax1b = plt.subplot2grid((9,1),(1,0),sharex=ax1a)
    ax1b.annotate('b - Growth Season Index (GSI)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1b.set_ylabel('GSI',fontsize = axis_size)
    
    ax1b.fill_between(model['time'],model['gsi'][:,3],model['gsi'][:,4],color=colour[2],alpha=0.2)
    ax1b.plot(model['time'],model['gsi'][:,1],'-',color=colour[2])
    
    # Plot c -> gsi photo
    ax1c = plt.subplot2grid((9,1),(2,0),sharex=ax1a,sharey=ax1b)
    ax1c.annotate('c - Growth Season Index (GSI) - photoperiod component', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1c.set_ylabel('GSI',fontsize = axis_size)
    
    ax1c.fill_between(model['time'],model['gsi_iphoto'][:,3],model['gsi_iphoto'][:,4],color=colour[2],alpha=0.2)
    ax1c.plot(model['time'],model['gsi_iphoto'][:,0],'-',color=colour[2])

    # Plot d -> gsi temp
    ax1d = plt.subplot2grid((9,1),(3,0),sharex=ax1a,sharey=ax1b)
    ax1d.annotate('d - Growth Season Index (GSI) - temperature component', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1d.set_ylabel('GSI',fontsize = axis_size)
    
    ax1d.fill_between(model['time'],model['gsi_itemp'][:,3],model['gsi_itemp'][:,4],color=colour[2],alpha=0.2)
    ax1d.plot(model['time'],model['gsi_itemp'][:,1],'-',color=colour[2])

    # Plot e -> gsi vpd
    ax1e = plt.subplot2grid((9,1),(4,0),sharex=ax1a,sharey=ax1b)
    ax1e.annotate('e - Growth Season Index (GSI) - VPD component', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1e.set_ylabel('GSI',fontsize = axis_size)
    
    ax1e.fill_between(model['time'],model['gsi_ivpd'][:,3],model['gsi_ivpd'][:,4],color=colour[2],alpha=0.2)
    ax1e.plot(model['time'],model['gsi_ivpd'][:,1],'-',color=colour[2])

    ax1b.set_ylim(0,1)
    
    # Plot f -> Clit
    ax1f = plt.subplot2grid((9,1),(5,0),sharex=ax1a)
    ax1f.annotate('f - C$_{litter}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1f.set_ylabel('C$_{lit}$ / g(C) m$^{-2}$',fontsize = axis_size)
    
    ax1f.fill_between(model['time'],model['Clit'][:,3],model['Clit'][:,4],color=colour[0],alpha=0.2)
    ax1f.plot(model['time'],model['Clit'][:,1],'-',color=colour[0])
    
    if 'Clit' in obs.keys(): # check for observations
        if 'Clit_u' in obs.keys(): # check for uncertainty bounds
            ax1f.errorbar(obs['time'],obs['Clit'],yerr=obs['Clit_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
        else:
            ax1f.plot(obs['time'],obs['Clit'],marker='o',c='black',mec='black',mfc='black')

    # Plot g -> Cfol into Clit
    ax1g = plt.subplot2grid((9,1),(6,0),sharex=ax1a)
    ax1g.annotate('f - litterfall', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1g.set_ylabel('litter flux / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    
    ax1g.fill_between(model['time'],model['flux_fol_lit'][:,3],model['flux_fol_lit'][:,4],color=colour[1],alpha=0.2)
    ax1g.fill_between(model['time'],model['flux_fol_lit'][:,5],model['flux_fol_lit'][:,3],color=colour[1],alpha=0.1)
    ax1g.fill_between(model['time'],model['flux_fol_lit'][:,4],model['flux_fol_lit'][:,6],color=colour[1],alpha=0.1)
    ax1g.plot(model['time'],model['flux_fol_lit'][:,1],'-',color=colour[1])

    if 'flux_fol_lit' in obs.keys(): # check for observations
        obs['flux_fol_lit']/obs['lit_acc_days']
        mask = np.isfinite(obs['flux_fol_lit'])
        x1 = ((obs['time'][mask]-obs['lit_acc_days'][mask].astype('int').astype('timedelta64[D]'))).astype(dt.datetime)
        x2 = obs['time'][mask].astype(dt.datetime)
        xmid = (obs['time'][mask]-(obs['lit_acc_days'][mask]/2).astype('int').astype('timedelta64[D]')).astype(dt.datetime)
        y=obs['flux_fol_lit'][mask]/obs['lit_acc_days'][mask]
        if 'flux_fol_lit_u' in obs.keys(): # check for uncertainty bounds
            y_e=obs['flux_fol_lit_u'][mask]/obs['lit_acc_days'][mask]
            for ii in range(0,mask.sum()):
                ax1g.plot([x1[ii],x2[ii]],[y[ii],y[ii]],'-',linewidth=5,color='black')
                ax1g.errorbar(xmid[ii],y[ii],yerr=y_e[ii],marker=None,c='black',ecolor='black',elinewidth=0.5)
        else:
            for ii in range(0,mask.sum()):
                ax1g.plot([x1[ii],x2[ii]],[y,y],'-',linewidth=1,color='black',ms=5)
    """
    if 'flux_fol_lit' in obs.keys(): # check for observations
        obs['flux_fol_lit']/obs['lit_acc_days']
        mask = np.isfinite(obs['flux_fol_lit'])
        x1 = obs['time'][mask]-obs['lit_acc_days'][mask]
        x2 = obs['time'][mask]
        y=obs['flux_fol_lit'][mask]/obs['lit_acc_days'][mask]
        if 'flux_fol_lit_u' in obs.keys(): # check for uncertainty bounds
            y_e=obs['flux_fol_lit_u'][mask]/obs['lit_acc_days'][mask]
            for ii in range(0,mask.sum()):
                ax1g.plot([x1[ii],x2[ii]],[y[ii],y[ii]],'-',linewidth=1,color='black')
                ax1g.errorbar((x1[ii]+x2[ii])/2.,y[ii],yerr=y_e[ii],marker='.',c='black',ecolor='black')
        else:
            for ii in range(0,mask.sum()):
                ax1g.plot([x1[ii],x2[ii]],[y,y],'-',linewidth=1,color='black')
    """
    # Plot h -> cwd fluxes into litter pool
    ax1h = plt.subplot2grid((9,1),(7,0),sharex=ax1a)
    ax1h.annotate('f - litter flux from CWD', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1h.set_ylabel('litter flux / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    
    ax1h.fill_between(model['time'],model['flux_cwd_lit'][:,3],model['flux_cwd_lit'][:,4],color=colour[1],alpha=0.2)
    ax1h.plot(model['time'],model['flux_cwd_lit'][:,1],'-',color=colour[1])


    # Plot i -> root fluxes into litter pool
    ax1i = plt.subplot2grid((9,1),(8,0),sharex=ax1a)
    ax1i.annotate('f - litter flux from roots', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1i.set_ylabel('litter flux / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    
    ax1i.fill_between(model['time'],model['flux_root_lit'][:,3],model['flux_root_lit'][:,4],color=colour[1],alpha=0.2)
    ax1i.plot(model['time'],model['flux_root_lit'][:,1],'-',color=colour[1])
    
    # set xlimits if desired
    if start_tstep!=False:
        ax1a.set_xlim(xmin=start_tstep)
    if end_tstep!=False:
        ax1a.set_xlim(xmax=end_tstep)
    else:
        ax1a.set_xlim(xmax=model['time'][-1])

    ax1a.set_ylim(0,12)
    ax1f.set_ylim(0,1000)
    ax1g.set_ylim(0,1.5)
    ax1h.set_ylim(0,5)
    ax1i.set_ylim(0,1)    
        
    plt.tight_layout()
    if len(figname)>0:
        plt.savefig(figname)

#---------------------------
# Plot second figure that compares the litter trap obserations against observed litter accumulation.
def plot_litter_trap_comparison(model,obs,start_tstep=False,end_tstep=False,figname=''):
    fig = plt.figure(4, facecolor='White',figsize=[5,5])
    n_steps = model['time'].size
    acc_fol_flux_mod = np.zeros(model['flux_fol_lit'].shape)*np.nan

    for tt in range(0,n_steps):
        if obs['lit_acc_days'][tt]>0:
            tt_start = tt-obs['lit_acc_days'][tt]
            #print tt_start,tt, obs['lit_acc_days'][tt], np.sum(model['flux_fol_lit'][tt_start:tt+1,:],axis=0)
            acc_fol_flux_mod[tt,:] = np.sum(model['flux_fol_lit'][tt_start:tt+1,:],axis=0)

    ax4 = plt.subplot2grid((1,1),(0,0))
    ax4.errorbar(obs['flux_fol_lit']/obs['lit_acc_days'],acc_fol_flux_mod[:,1]/obs['lit_acc_days'],xerr=obs['flux_fol_lit_u']/obs['lit_acc_days'],marker='o',c='black',mec='black',mfc='black',ecolor='black')
    ax4.set_ylabel('modelled litter influx / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    ax4.set_xlabel('observed litter influx / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)

    mask = np.isfinite(obs['flux_fol_lit'])
    m, c, r, p, err = stats.linregress(obs['flux_fol_lit'][mask]/obs['lit_acc_days'][mask],acc_fol_flux_mod[:,1][mask]/obs['lit_acc_days'][mask])
    pstr=''
    if p<0.001:
        p1str='***'
    elif p<0.01:
        pstr='** '
    elif p<0.05:
        pstr='*  '
    elif p<0.1:
        pstr='$^.$  '
    else:
        pstr='   '

    stats_str = 'R$^2$=' + '%.3f%s' % (r**2, pstr)
    ax4.annotate(stats_str, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=rcParams['font.size'])

    fit_x = np.array([np.min(obs['flux_fol_lit'][mask]/obs['lit_acc_days'][mask]), np.max(obs['flux_fol_lit'][mask]/obs['lit_acc_days'][mask])])
    fit_y = m*fit_x+c
    ax4.plot(fit_x,fit_y,':k')
    ax4.axis('equal')
    ax4.set_xlim(0,6)
    ax4.set_ylim(0,6)
    plt.tight_layout()

    if len(figname)>0:
        plt.savefig(figname)



#================================================
# plot_parameters
#------------------------------------------------
# a function to plot histograms of the accepted parameters
# inputs:
# - a parameter array
# - figname (optional)
def plot_parameters(params,figname=''):
    fig = plt.figure(5, facecolor='White',figsize=[18,9])

    par_list = ['Decomposition rate',
                'Fraction GPP respired',
                'GSI sensitivity growth',
                'Root allocation',
                'GSI max leaf turnover',
                'TOR wood / d$^{-1}$',
                'TPR roots / d$^{-1}$',
                'TOR litter / d$^{-1}$',
                'TOR SOM / d$^{-1}$',
                'Q10',
                'foliar N / gN m${-2}$',
                'GSI max labile turnover',
                'Labile allocation',
                'GSI min temp / $^o$C',
                'GSI max temp / $^o$C',
                'GSI min photoperiod / s',
                'LMA / gC m$^{-2}$',
                'Clab$_i$',
                'Cfol$_i$',
                'Croot$_i$',
                'Cwood$_i$',
                'Clit$_i$',
                'Csom$_i$',
                'GSI max photoperiod / s',
                'GSI min VPD / Pa',
                'GSI max VPD / Pa',
                'crit. GPP for LAI',
                'frac Cwood branch',
                'frac Cwood c. root',
                'replant Clab',
                'replant Cfol',
                'replant Croot',
                'replant Cwood at age 1',
                'GSI sensitivity senescence',
                'GSI - have I just left growing state',
                'GSI$_i$.',
                'Ccwd$_i$',
                'TOR CWD / d$^{-1}$',
                'min TOR fol / d$^{-1}$',
                'likelihood']
    
    n_par = params.shape[0]
    rows = 5
    cols = 8
    par = 0
    for rr in range(0,rows):
        for cc in range(0,cols):
            if par >= n_par:
                break
            else:
                axpp = plt.subplot2grid((rows,cols),(rr,cc))
                if par == 11:
                    axpp.hist(10**params[par], bins = 20, normed=True, fc='#46E900', ec='#46E900')
                elif np.all((par >= 29, par<33)):
                    axpp.hist(10**params[par], bins = 20, normed=True, fc='0.5', ec='0.5')
                elif par==n_par-1:
                    axpp.hist(params[par], bins = 20, normed=True, fc=colour[2], ec=colour[2])
                    
                else:
                    axpp.hist(params[par], bins = 20, normed=True, fc='#46E900', ec='#46E900')
                axpp.annotate(par_list[par], xy=(0.5,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='center', verticalalignment='top', fontsize=rcParams['font.size'])
                axpp.xaxis.set_major_locator(plt.MaxNLocator(3))
                axpp.yaxis.set_major_locator(plt.NullLocator())
            par += 1
                
    plt.tight_layout()
    if len(figname)>0:
        plt.savefig(figname)
            
    #plt.show()


#----------------------
# plot_summary_ts
# ---------------------
# Plot a time series of the model output, summarising key components (for poster)
# The function reads in two dictionaries:
# - the model output (median, ulim, llim)
# - the obserations where available
# Note that these dictionaries should have uniform naming structures
# - lai : leaf area index
# - gsi : growth season index
# - gpp : Gross Primary Production
# - npp : Net Primary Production
# - Cwoo : woody carbon stocks
# - flux_fol_lit : carbon flux from foliar pool to litter (litterfall)
# Also takes optional arguments for start and end timestep. initially these will
# be index references (i.e. model timestep) but this will ultimately be altered
# to give options to specify date ranges.
def plot_summary_ts(model,obs,met,start_tstep=False,end_tstep=False,figname=''):
    sns.set()
    date = model['time'].astype('O')
    fig = plt.figure(3, facecolor='White',figsize=[11,12])

    # Plot a -> Precipitation
    ax1a = plt.subplot2grid((4,2),(0,0))
    ax1a.annotate('a - precipitation', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1a.set_ylabel('daily precipitation / mm',fontsize = axis_size)
    
    ax1a.fill_between(date,np.zeros(date.size),met['prcp'],color=colour[1], linewidth=0.0)

    # Plot b -> VPD
    ax1b = plt.subplot2grid((4,2),(1,0),sharex=ax1a)
    ax1b.annotate('c - Vapour pressure deficit', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1b.set_ylabel('VPD / Pa',fontsize = axis_size)
    
    ax1b.plot(date,met['VPD'],color=colour[1])
    
    # Plot a -> LAI
    ax1c = plt.subplot2grid((4,2),(0,1),sharex=ax1a)
    ax1c.annotate('b - LAI', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1c.set_ylabel('LAI / m$^2$m$^{-2}$',fontsize = axis_size)
    
    ax1c.fill_between(date,model['lai'][:,3],model['lai'][:,4],color=colour[0],alpha=0.5, linewidth=0.0)
    ax1c.fill_between(date,model['lai'][:,3],model['lai'][:,5],color=colour[0],alpha=0.25, linewidth=0.0)
    ax1c.fill_between(date,model['lai'][:,6],model['lai'][:,4],color=colour[0],alpha=0.25, linewidth=0.0)
    ax1c.plot(date,model['lai'][:,1],'-',color=colour[0])

    if 'lai' in obs.keys(): # check for observations
        if 'lai' in obs.keys(): # check for uncertainty bounds
            ax1c.errorbar(date,obs['lai'],yerr=obs['lai_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black',elinewidth=0.5,ms=5)
        else:
            ax1c.plot(date,obs['lai'],marker='o',c='black',mec='black',mfc='black',ms=5)
    
    # Plot b -> gsi
    ax1d = plt.subplot2grid((4,2),(2,0),sharex=ax1a)
    ax1d.annotate('e - Growth Season Index (GSI)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1d.set_ylabel('GSI',fontsize = axis_size)
    ax1d.fill_between(date,model['gsi'][:,3],model['gsi'][:,4],color=colour[2],alpha=0.5, linewidth=0.0)
    ax1d.fill_between(date,model['gsi'][:,3],model['gsi'][:,5],color=colour[2],alpha=0.25, linewidth=0.0)
    ax1d.fill_between(date,model['gsi'][:,6],model['gsi'][:,4],color=colour[2],alpha=0.25, linewidth=0.0)
    ax1d.plot(date,model['gsi'][:,1],'-',color=colour[2])
    # Plot c -> gpp
    ax1e = plt.subplot2grid((4,2),(1,1))#,sharex=ax1a)
    ax1e.annotate('d - Gross Primary Productivity', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1e.set_ylabel('GPP / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    ax1e.fill_between(date,model['gpp'][:,3],model['gpp'][:,4],color=colour[0],alpha=0.5, linewidth=0.0)
    ax1e.fill_between(date,model['gpp'][:,3],model['gpp'][:,5],color=colour[0],alpha=0.25, linewidth=0.0)
    ax1e.fill_between(date,model['gpp'][:,6],model['gpp'][:,4],color=colour[0],alpha=0.25, linewidth=0.0)
    ax1e.plot(date,model['gpp'][:,1],'-',color=colour[0])

    if 'gpp' in obs.keys(): # check for observations
        if 'gpp_u' in obs.keys(): # check for uncertainty bounds
            ax1e.errorbar(date,obs['gpp'],yerr=obs['gpp_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black',elinewidth=0.5,ms=5)
        else:
            ax1e.plot(date,obs['gpp'],marker='o',c='black',mec='black',mfc='black',ms=5)

    # Plot d -> NPP
    ax1f = plt.subplot2grid((4,2),(2,1),sharex=ax1a)
    ax1f.annotate('f - Net Primary Productivity', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1f.set_ylabel('NPP / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    
    ax1f.fill_between(date,model['npp'][:,3],model['npp'][:,4],color=colour[0],alpha=0.5, linewidth=0.0)
    ax1f.fill_between(date,model['npp'][:,3],model['npp'][:,5],color=colour[0],alpha=0.25, linewidth=0.0)
    ax1f.fill_between(date,model['npp'][:,6],model['npp'][:,4],color=colour[0],alpha=0.25, linewidth=0.0)
    ax1f.plot(date,model['npp'][:,1],'-',color=colour[0])
    
    if 'npp' in obs.keys(): # check for observations
        if 'npp_u' in obs.keys(): # check for uncertainty bounds
            ax1f.errorbar(date,obs['npp'],yerr=obs['npp_u'],marker='o',c='black',mec='black',mfc='black',ecolor='black',elinewidth=0.5,ms=5)
        else:
            ax1f.plot(date,obs['npp'],marker='o',c='black',mec='black',mfc='black',ms=5)
    
    # Plot e -> litterfall
    ax1g = plt.subplot2grid((4,2),(3,0),sharex=ax1a)
    ax1g.annotate('g - litterfall', xy=(0.05,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='bottom', fontsize=10)
    ax1g.set_ylabel('litter flux / g(C) m$^{-2}$ d$^{-1}$',fontsize = axis_size)
    
    ax1g.fill_between(date,model['flux_fol_lit'][:,3],model['flux_fol_lit'][:,4],color=colour[1],alpha=0.5, linewidth=0.0)
    ax1g.fill_between(date,model['flux_fol_lit'][:,3],model['flux_fol_lit'][:,5],color=colour[1],alpha=0.25, linewidth=0.0)
    ax1g.fill_between(date,model['flux_fol_lit'][:,6],model['flux_fol_lit'][:,4],color=colour[1],alpha=0.25, linewidth=0.0)
    ax1g.plot(date,model['flux_fol_lit'][:,1],'-',color=colour[1])
    
    if 'flux_fol_lit' in obs.keys(): # check for observations
        obs['flux_fol_lit']/obs['lit_acc_days']
        mask = np.isfinite(obs['flux_fol_lit'])
        x1 = ((obs['time'][mask]-obs['lit_acc_days'][mask].astype('int').astype('timedelta64[D]'))).astype(dt.datetime)
        x2 = obs['time'][mask].astype(dt.datetime)
        xmid = (obs['time'][mask]-(obs['lit_acc_days'][mask]/2).astype('int').astype('timedelta64[D]')).astype(dt.datetime)
        y=obs['flux_fol_lit'][mask]/obs['lit_acc_days'][mask]
        if 'flux_fol_lit_u' in obs.keys(): # check for uncertainty bounds
            y_e=obs['flux_fol_lit_u'][mask]/obs['lit_acc_days'][mask]
            for ii in range(0,mask.sum()):
                ax1g.plot([x1[ii],x2[ii]],[y[ii],y[ii]],'-',linewidth=5,color='black')
                ax1g.errorbar(xmid[ii],y[ii],yerr=y_e[ii],marker=None,c='black',ecolor='black',elinewidth=0.5)
        else:
            for ii in range(0,mask.sum()):
                ax1g.plot([x1[ii],x2[ii]],[y,y],'-',linewidth=1,color='black',ms=5)
    
    
    # Plot f -> Woody biomass
    ax1h = plt.subplot2grid((4,2),(3,1),sharex=ax1a)
    ax1h.annotate('h - C$_{wood}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1h.set_ylabel('C$_{wood}$ / Mg(C) ha$^{-2}$',fontsize = axis_size)
    
    ax1h.fill_between(date,model['Cwoo'][:,3]/100.,model['Cwoo'][:,4]/100.,color=colour[0],alpha=0.5, linewidth=0.0)
    ax1h.fill_between(date,model['Cwoo'][:,3]/100.,model['Cwoo'][:,5]/100.,color=colour[0],alpha=0.25, linewidth=0.0)
    ax1h.fill_between(date,model['Cwoo'][:,6]/100.,model['Cwoo'][:,4]/100.,color=colour[0],alpha=0.25, linewidth=0.0)
    ax1h.plot(date,model['Cwoo'][:,1]/100.,'-',color=colour[0])

    if 'Cwoo' in obs.keys(): # check for observations
        if 'Cwoo_u' in obs.keys(): # check for uncertainty bounds
            ax1h.errorbar(date,obs['Cwoo']/100.,yerr=obs['Cwoo_u']/100.,marker='o',c='black',mec='black',mfc='black',ecolor='black',elinewidth=0.5,ms=5)
        else:
            ax1h.plot(date,obs['Cwoo']/100.,marker='o',c='black',mec='black',mfc='black',ms=5)

    # set xlimits if desired
    if start_tstep!=False:
        ax1a.set_xlim(xmin=start_tstep)
    if end_tstep!=False:
        ax1a.set_xlim(xmax=end_tstep)
    else:
        ax1a.set_xlim(xmin=date[0],xmax=date[-1])
  
    ax1g.set_ylim(0,1.5) 
    ax1d.set_ylim(0,1) 
    """
    ax1a.tick_params(axis='x',labelbottom='off')
    ax1b.tick_params(axis='x',labelbottom='off')
    ax1c.tick_params(axis='x',labelbottom='off')
    ax1d.tick_params(axis='x',labelbottom='off')
    """
    plt.tight_layout()
    if len(figname)>0:
        plt.savefig(figname)


# PLOT ALLOCATION FRACTIONS
# Plot second figure that compares the litter trap obserations against observed litter accumulation.
def plot_allocation_fractions(df,figname=''):
    sns.set_style("white")
    fig = plt.figure(6, facecolor='White',figsize=[10,3])
    # axis one - allocation to wood
    ax6a = plt.subplot2grid((1,3),(0,0))
    ax6a.annotate('a - Fraction of NPP allocated to roots', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax6a.set_ylabel('allocation fraction',fontsize = axis_size)
    sns.violinplot(x='plot',y='root',data=df,inner=None,color='white',linewidth=0.5)
    sns.stripplot(x='plot',y='root',data=df,jitter=0.05,alpha=0.005,linewidth=None)

    # axis two - allocation to foliage
    ax6b = plt.subplot2grid((1,3),(0,1),sharey=ax6a)
    ax6b.annotate('b - Fraction NPP allocated to canopy', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

    sns.violinplot(x='plot',y='canopy',data=df,inner=None,color='white',linewidth=0.5)
    sns.stripplot(x='plot',y='canopy',data=df,jitter=0.05,alpha=0.005,linewidth=None)
    
    # axis three - allocation to foliage
    ax6c = plt.subplot2grid((1,3),(0,2),sharey=ax6a)
    ax6c.annotate('c - Fraction of NPP allocated to wood', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    
    ax6b.tick_params(axis='y',labelbottom='off')
    ax6c.tick_params(axis='y',labelbottom='off')

    sns.violinplot(x='plot',y='wood',data=df,inner=None,color='white',linewidth=0.5)
    sns.stripplot(x='plot',y='wood',data=df,jitter=0.05,alpha=0.005,linewidth=None)

    ax6a.set_ylim((0,1))
    
    ax6a.set_xlabel('')
    ax6b.set_xlabel('')
    ax6c.set_xlabel('')
    for lab  in ax6a.get_xticklabels():
        lab.set_rotation(45)
    for lab  in ax6b.get_xticklabels():
        lab.set_rotation(45)
    for lab  in ax6c.get_xticklabels():
        lab.set_rotation(45)
    plt.tight_layout()
    
    if len(figname)>0:
        plt.savefig(figname)

    plt.show()


# PLOT LEAF TRAITS
# Plot second figure that compares the litter trap obserations against observed litter accumulation.
def plot_leaf_traits(df,figname=''):
    sns.set_style("white")
    fig = plt.figure(7, facecolor='White',figsize=[10,3])
    # axis one - allocation to wood
    ax7a = plt.subplot2grid((1,3),(0,0))
    ax7a.annotate('a - LMA', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    sns.violinplot(x='plot',y='LMA',data=df,inner=None,color='white',linewidth=0.5)
    sns.stripplot(x='plot',y='LMA',data=df,jitter=0.05,alpha=0.005,linewidth=None)
    ax7a.set_ylabel('LMA / g(C) m$^{-2}$',fontsize = axis_size)
    ax7a.set_xlabel('')

    # axis two - Narea
    ax7b = plt.subplot2grid((1,3),(0,1))
    ax7b.annotate('b - Leaf Nitrogen, [N]$_{area}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

    sns.violinplot(x='plot',y='Narea',data=df,inner=None,color='white',linewidth=0.5)
    sns.stripplot(x='plot',y='Narea',data=df,jitter=0.05,alpha=0.005,linewidth=None)
    
    ax7b.set_ylabel('[N]$_{area}$ / g(N) m$^{-2}$',fontsize=10)
    ax7b.set_xlabel('')
    
    # axis three - C:N ratio
    ax7c = plt.subplot2grid((1,3),(0,2))
    ax7c.annotate('c - C:N ratio', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

    sns.violinplot(x='plot',y='CNratio',data=df,inner=None,color='white',linewidth=0.5)
    sns.stripplot(x='plot',y='CNratio',data=df,jitter=0.05,alpha=0.005,linewidth=None)
    ax7c.set_ylabel('C:N ratio')
    ax7c.set_xlabel('')

    ax7a.set_ylim((0,200))
    ax7b.set_ylim((0,10**0.7))
    ax7c.set_ylim(ymin=0)
    
    for lab  in ax7a.get_xticklabels():
        lab.set_rotation(45)
    for lab  in ax7b.get_xticklabels():
        lab.set_rotation(45)
    for lab  in ax7c.get_xticklabels():
        lab.set_rotation(45)
    plt.tight_layout()
    
    if len(figname)>0:
        plt.savefig(figname)

    plt.show()


    

# PLOT LEAF TRAITS
# Plot second figure that compares the litter trap obserations against observed litter accumulation.
def plot_residence_times(df,figname=''):
    sns.set_style("white")
    fig = plt.figure(8, facecolor='White',figsize=[10,3])
    # axis one - allocation to wood
    ax8a = plt.subplot2grid((1,3),(0,0))
    ax8a.annotate('a - Roots', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    sns.violinplot(x='plot',y='root_rt',data=df,inner=None,color='white',linewidth=0.5)
    sns.stripplot(x='plot',y='root_rt',data=df,jitter=0.05,alpha=0.005,linewidth=None)
    
    ax8a.set_ylabel('Residence time / yrs',fontsize = axis_size)
    
    # axis two - Narea
    ax8b = plt.subplot2grid((1,3),(0,1))
    ax8b.annotate('b - Foliage', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

    sns.violinplot(x='plot',y='fol_rt',data=df,inner=None,color='white',linewidth=0.5)
    sns.stripplot(x='plot',y='fol_rt',data=df,jitter=0.05,alpha=0.005,linewidth=None)
    
    ax8b.set_ylabel('Residence time / yrs',fontsize=10)
    
    # axis three - C:N ratio
    ax8c = plt.subplot2grid((1,3),(0,2))
    ax8c.annotate('c - Wood', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    sns.violinplot(x='plot',y='wood_rt',data=df,inner=None,color='white',linewidth=0.5)
    sns.stripplot(x='plot',y='wood_rt',data=df,jitter=0.05,alpha=0.005,linewidth=None)

    ax8c.set_ylabel('Residence time / yrs')

    ax8a.set_ylim(ymin=0)
    ax8b.set_ylim(ymin=0)
    ax8c.set_ylim(ymin=0)
    ax8a.set_xlabel('')
    ax8b.set_xlabel('')
    ax8c.set_xlabel('')
    for lab  in ax8a.get_xticklabels():
        lab.set_rotation(45)
    for lab  in ax8b.get_xticklabels():
        lab.set_rotation(45)
    for lab  in ax8c.get_xticklabels():
        lab.set_rotation(45)
        
    plt.tight_layout()
    
    if len(figname)>0:
        plt.savefig(figname)

    plt.show()


    
# PLOT ALLOCATION FRACTIONS FOR SINGLE SITE
def plot_allocation_fractions_single_site(df,plot,figname=''):
    sns.set_style("white")
    fig = plt.figure(9, facecolor='White',figsize=[4,4])
    # axis one - allocation to wood
    ax= plt.subplot2grid((1,1),(0,0))
    ax.annotate(plot, xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax.set_ylabel('allocation fraction',fontsize = axis_size)
    sns.violinplot(x='pool',y='allocation fraction',data=df,inner=None,linewidth=0.5,palette="muted")
    sns.stripplot(x='pool',y='allocation fraction',data=df,color='black',jitter=0.05,alpha=0.005,linewidth=None)

    ax.set_ylim((0,1))
    ax.set_xlabel('')
    
    for lab  in ax.get_xticklabels():
        lab.set_rotation(45)
    plt.tight_layout()
    
    if len(figname)>0:
        plt.savefig(figname)

    plt.show()


# PLOT RESIDENCE TIMES FOR SINGLE SITE
def plot_residence_times_single_site(df,plot,figname=''):
    sns.set_style("white")
    fig = plt.figure(10, facecolor='White',figsize=[4,4])
    # axis one - allocation to wood
    ax= plt.subplot2grid((1,1),(0,0))
    ax.annotate(plot, xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax.set_ylabel('allocation fraction',fontsize = axis_size)
    sns.violinplot(x='pool',y='residence time',data=df,inner=None,linewidth=0.5,scale="width",palette="muted")
    sns.stripplot(x='pool',y='residence time',data=df,color='black',jitter=0.05,alpha=0.005,linewidth=None)

    ax.set_xlabel('')
    ax.set_ylim(ymin=0.1)
    ax.set_yscale('log')
    for lab  in ax.get_xticklabels():
        lab.set_rotation(45)
    plt.tight_layout()
    
    if len(figname)>0:
        plt.savefig(figname)

    plt.show()
