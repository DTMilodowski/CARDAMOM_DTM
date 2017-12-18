#This set of functions prepares field data from GEM plots ready for assimilation into CARDAMOM
import numpy as np
import sys
import load_GEM_data as field

from scipy.interpolate import interp1d

# Get time series of woody biomass
def get_Cwood_ts(census_file,plot):
    census = field.collate_plot_level_census_data(census_file)
    plot_biomass = np.sum(census[plot]['C_wood'],axis=0)
    collection_date = np.max(census[plot]['CensusDate'],axis=0)
    plot_biomass[collection_date<np.datetime64('2000-01-01','D')]=np.nan
    return collection_date[np.isfinite(plot_biomass)], plot_biomass[np.isfinite(plot_biomass)]

# Get time series of fine root NPP
def get_root_NPP_ts(roots_file,plot,pad_ts=True):
    rootStocks,rootNPP = field.read_soil_stocks_and_npp(roots_file)
    N_cores,N_dates = rootNPP[plot]['AccumulationDays'].shape
    collection_dates = rootNPP[plot]['CollectionDate']
    previous_collection_dates = np.max(rootNPP[plot]['PreviousCollectionDate'],axis=0)

    interval = np.asarray(collection_dates-previous_collection_dates,dtype='float64')
    days = np.cumsum(interval)

    rootNPP_gapfilled = np.zeros((N_cores,N_dates))
    for ss in range(0,N_cores):
        # First check to see if there are gaps - if not, don't need to worry
        if (np.isnan(rootNPP[plot]['FineRootNPP'][ss,:])).sum()==0:
            rootNPP_gapfilled[ss,:]=rootNPP[plot]['FineRootNPP'][ss,:].copy()

        # We don't want to gapfill at the start or end of the time series
        # as we have no other constraints for the interpolation
        else:
            rootNPP_gapfilled[ss,:]=gapfill_field_data(rootNPP[plot]['FineRootNPP'][ss,:],days,pad_ts=pad_ts)

    rootNPP_gapfilled[rootNPP_gapfilled<0]=0

    rootNPP_ts = np.mean(rootNPP_gapfilled,axis=0)
    rootNPP_std = np.std(rootNPP_gapfilled,axis=0)
    return collection_dates, previous_collection_dates, rootNPP_ts, rootNPP_std

# Note that later root stocks surveys are pretty irregular - suggest only using the first survey as these are comprehensive.
def get_Croot(roots_file,plot):

    rootStocks,rootNPP = field.read_soil_stocks_and_npp(roots_file)
    N_core,N_dates = rootStocks[plot]['FineRootStocks'].shape

    Croot_plot = np.mean(rootStocks[plot]['FineRootStocks'][:,0])
    Croot_std = np.std(rootStocks[plot]['FineRootStocks'][:,0])
    collection_date = rootStocks[plot]['CollectionDate'][0]
    return collection_date, Croot_plot, Croot_std

# Get litterfall time series.  If pad_ts set to True (default), then nodata values at the ends of the time series for certain
# subplots will be padded with the first or last recorded value, so that average litter fluxes are taken across all subplots.
# This is an attempt to avoid biases in the averages, given that each 1 ha plot comprises only 25 subplots.
# Litter fluxes are returned as the total litter fall [g(C) m-2] collected within a specified time period [accumulation days]
# alongside collection dates
def get_litterfall_ts(litter_file,plot, pad_ts = True):
    
    litter = field.read_litterfall_data(litter_file)
    N_sp,N_dates = litter[plot]['mTotal'].shape
    
    acc_mass = litter[plot]['mTotal']/litter[plot]['TrapSize'] # convert from g(C) to g(C) m-2
    flux = litter[plot]['rTotal']*10.**6/10.**4/365.25 # convert flux from Mg(C)ha-1yr-1 to g(C)m-2d-1
    
    collection_dates = np.max(litter[plot]['CollectionDate'],axis=0)
    previous_collection_dates = np.max(litter[plot]['PreviousCollectionDate'],axis=0)

    accumulation_days = np.asarray(collection_dates-previous_collection_dates,dtype='float64')
    days = np.asarray(collection_dates-collection_dates[0],dtype='float64')
    
    litter_gapfilled = np.zeros((N_sp,N_dates))
    for ss in range(0,N_sp):
        # First check to see if there are gaps - if not, don't need to worry
        if (np.isnan(litter[plot]['mTotal'][ss,:])).sum()==0:
            litter_gapfilled[ss,:]=acc_mass.copy()

        # We don't want to gapfill at the start or end of the time series
        # as we have no other constraints for the interpolation
        else:
            gapfilled_fluxes=gapfill_field_data(litter[plot]['rTotal'][ss,:],days,pad_ts=pad_ts)
            litter_gapfilled[ss,:] = accumulation_days*gapfilled_fluxes # convert back to total accumulated C

    litter_gapfilled[litter_gapfilled<0]=0

    litter_fall_ts = np.mean(litter_gapfilled,axis=0)
    litter_fall_std = np.std(litter_gapfilled,axis=0)
    litter_fall_serr = litter_fall_std/np.sqrt(float(N_sp))
    
    return collection_dates, accumulation_days, litter_fall_ts, litter_fall_std, litter_fall_serr

"""
this needs updating
def get_subplot_litterfall_ts(litter_file,plot, pad_ts = True):
    
    litter = field.read_litterfall_data(litter_file)
    N_sp,N_dates = litter[plot]['rTotal'].shape

    collection_dates = np.max(litter[plot]['CollectionDate'],axis=0)
    previous_collection_dates = np.max(litter[plot]['PreviousCollectionDate'],axis=0)

    interval = np.asarray(collection_dates-previous_collection_dates,dtype='float64')
    days = np.cumsum(interval)
    
    litter_gapfilled = np.zeros((N_sp,N_dates))
    for ss in range(0,N_sp):
        # First check to see if there are gaps - if not, don't need to worry
        if (np.isnan(litter[plot]['rTotal'][ss,:])).sum()==0:
            litter_gapfilled[ss,:]=litter[plot]['rTotal'][ss,:].copy()

        # We don't want to gapfill at the start or end of the time series
        # as we have no other constraints for the interpolation
        else:
            litter_gapfilled[ss,:]=gapfill_field_data(litter[plot]['rTotal'][ss,:],days,pad_ts=pad_ts)

    litter_gapfilled[litter_gapfilled<0]=0
    litter_fall_ts = litter_gapfilled.copy()

    return collection_dates, accumulation_days, litter_fall_ts
"""

# Get time series of LAI using spline interpolation to fill the gaps
def get_LAI_ts(LAI_file,plot, pad_ts = True):
    LAI = field.load_LAI_time_series(LAI_file)
    N_sp, N_dates = LAI[plot]['LAI'].shape
    interval = np.zeros(N_dates,'timedelta64[D]')
    interval[1:] = LAI[plot]['date'][1:]-LAI[plot]['date'][:-1]
    days = np.cumsum(interval)
    indices = np.arange(0,N_dates,dtype='int')

    LAI_gapfilled = np.zeros((N_sp,N_dates))*np.nan
    for ss in range(0,N_sp):
        # First check to see if there are gaps - if not, don't need to worry
        if (np.isnan(LAI[plot]['LAI'][ss,:])).sum()==0:
            LAI_gapfilled[ss,:]=LAI[plot]['LAI'][ss,:].copy()

        # We don't want to gapfill at the start or end of the time series
        # as we have no other constraints for the interpolation
        else:
            LAI_gapfilled[ss,:]=gapfill_field_data(LAI[plot]['LAI'][ss,:],days,pad_ts=pad_ts)

    LAI_plot_ts = np.mean(LAI_gapfilled,axis=0)
    LAI_plot_std_ts = np.std(LAI_gapfilled,axis=0)
    LAI_plot_serr_ts = np.std(LAI_gapfilled,axis=0)/np.sqrt(float(N_sp))
    
    return  LAI[plot]['date'], LAI_plot_ts, LAI_plot_std_ts

def get_subplot_LAI_ts(LAI_file,plot, pad_ts = True):
    LAI = field.load_LAI_time_series(LAI_file)
    N_sp, N_dates = LAI[plot]['LAI'].shape
    interval = np.zeros(N_dates,'timedelta64[D]')
    interval[1:] = LAI[plot]['date'][1:]-LAI[plot]['date'][:-1]
    days = np.cumsum(interval)
    indices = np.arange(0,N_dates,dtype='int')

    LAI_gapfilled = np.zeros((N_sp,N_dates))*np.nan
    for ss in range(0,N_sp):
        # First check to see if there are gaps - if not, don't need to worry
        if (np.isnan(LAI[plot]['LAI'][ss,:])).sum()==0:
            LAI_gapfilled[ss,:]=LAI[plot]['LAI'][ss,:].copy()

        # We don't want to gapfill at the start or end of the time series
        # as we have no other constraints for the interpolation
        else:
            LAI_gapfilled[ss,:]=gapfill_field_data(LAI[plot]['LAI'][ss,:],days,pad_ts=pad_ts)

    LAI_plot_ts = LAI_gapfilled.copy()
    return  LAI[plot]['date'], LAI_plot_ts


# A generic gapfilling script. array is a 1D array with input data to be gapfilled
def gapfill_field_data(array,tsteps,pad_ts=True):

    indices = np.arange(0,array.size,dtype='int')
    # find first and last datapoint in time series
    first = indices[np.isfinite(array)][0]
    last = indices[np.isfinite(array)][-1]
    gapfilled = np.zeros(array.size)

    if (np.isnan(array[first:last+1])).sum()==0:
        gapfilled=array.copy()
    else:
        mask = np.isfinite(array[first:last+1])
        tsteps_nogaps = np.asarray(tsteps[first:last+1][mask],dtype='float')
        array_nogaps = array[first:last+1][mask]
        f = interp1d(tsteps_nogaps, array_nogaps, kind='linear')  # without specifying "kind", default is linear
        gapfilled[first:last+1] = f(tsteps[first:last+1])
        
    # for now, pad the time series with constant value where required so that plot average can be obtained
    if pad_ts == True:
        if first>0:
            gapfilled[:first] = gapfilled[first]
        if last<indices[-1]:
            gapfilled[last+1:] = gapfilled[last]

    return gapfilled
