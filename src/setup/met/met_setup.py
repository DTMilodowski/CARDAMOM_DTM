# This set of functions prepares met data for inclusion into DALEC runs
# Includes ERA interim and TRMM
import numpy as np

import sys
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/ERAinterim/")
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/TRMM/")
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/gapfilling/")
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/gapfilling/")

import weather_generator_BALI as wg
import TRMM_processing as TRMM
from metdata_processing import *

sys.path.append('/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/construct_drivers/')
import gapfill_station_metdata as gap

# generate daily met drivers from ERA-Interim and TRMM
# currently do not use shift in time zone, so min, max
# mean and cumulative variables represent GMT "days"
def generate_daily_met_drivers_ERAinterim_TRMM(ERA_file, TRMM_file, start_date, end_date):
    
   # read in ERA interim.
    ERA_start = '2011-01-01'
    ERA_start_date = np.datetime64(ERA_start,'D')# - np.timedelta64(time_zone_shift,'h')

    sw_eraint, tair_eraint, pptn_eraint, rh_eraint, sp_eraint = wg.readInEraInterim(ERA_file,1) 
    sw_era = sw_eraint[:sw_eraint.size/2].astype('float')  
    pptn_era = pptn_eraint[:sw_eraint.size/2].astype('float')
    airT_era = tair_eraint.astype('float')
    rh_era = rh_eraint.astype('float')
    sp_era = sp_eraint.astype('float')
    vpd_era=(0.6108*np.exp((17.27*airT_era)/(237.3+airT_era)))*(1-rh_era/100)

    n_days_eraint = tair_eraint.size/4
    
    mn2t = np.zeros(n_days_eraint)
    mx2t = np.zeros(n_days_eraint)
    ssrd = np.zeros(n_days_eraint)
    vpd = np.zeros(n_days_eraint)
    ERA_dates = np.zeros(n_days_eraint,dtype='datetime64[D]')
    iter_2 = 2
    iter_4 = 4
    for dd in range(0,n_days_eraint):
        # airT and rh both reported quarterly
        mn2t[dd]=np.min(airT_era[dd*iter_4:(dd+1)*iter_4])
        mx2t[dd]=np.max(airT_era[dd*iter_4:(dd+1)*iter_4])
        vpd[dd] = np.mean(vpd_era[dd*iter_4:(dd+1)*iter_4])

        # sw rad reported as 12hr totals
        ssrd[dd]=np.sum(sw_era[dd*iter_2:(dd+1)*iter_2])
    
        ERA_dates[dd] = ERA_start_date+np.timedelta64(1,'D')*dd

    ###############
    # Now load TRMM - TRMM is reported as three hourly totals
    print TRMM_file
    TRMM_dates_init, TRMM_pptn_init = TRMM.read_TRMM_file(TRMM_file)

    TRMM_dates_only =  TRMM_dates_init.astype('datetime64[D]')
    TRMM_dates = np.unique(TRMM_dates_only)
    N_TRMM = TRMM_dates.size
    TRMM_pptn = np.zeros(N_TRMM)
    for dd in range(0,N_TRMM):
        TRMM_pptn[dd] = np.sum(TRMM_pptn_init[TRMM_dates_only==TRMM_dates[dd]])
    # Now cut time series to period of interest
    start = np.datetime64(start_date,'D')
    end = np.datetime64(end_date,'D')
    if TRMM_dates[0]>start:
        print "Issue with TRMM - time series starts too late for period of interest"
    if TRMM_dates[-1]<end:
        print "Issue with TRMM - time series finishes too early for period of interest"
    if ERA_dates[0]>start:
        print "Issue with ERA-Interim - time series starts too late for period of interest"
    if ERA_dates[-1]<end:
        print "Issue with ERA-Interim - time series finishes too early for period of interest"
    TRMM_mask = np.all((TRMM_dates>=start,TRMM_dates<=end),axis=0)
    ERA_mask = np.all((ERA_dates>=start,ERA_dates<=end),axis=0)
    
    dates_out = ERA_dates[ERA_mask]
    ssrd_out = ssrd[ERA_mask]
    vpd_out = vpd[ERA_mask]
    mn2t_out = mn2t[ERA_mask]
    mx2t_out = mx2t[ERA_mask]
    TRMM_dates_out = TRMM_dates[TRMM_mask]
    pptn_out = TRMM_pptn[TRMM_mask]
    return dates_out,mn2t_out,mx2t_out,vpd_out,ssrd_out,TRMM_dates_out,pptn_out

# This downsamples a pre-existing half-hourly gapfilled dataset to a daily timestep
# the input file would usually be generated using the SPA prep scripts
def generate_daily_met_drivers_from_existing_halfhourly_time_series(met_data_file,start,end):
    dtype={'names':('time_in_days','time','airt','wind_spd','sw_rad','vpd','ppfd','precip','sfc_pressure','gap_airt','gap_sw','gap_vpd','gap_ppfd','gap_pptn','gap_sfc_pressure'),'formats':('f32','S32','f32','f32','f32','f32','f32','f32','f32','f32','f32','f32','f32')}
    met_data = np.genfromtxt(met_data_file,skiprows=1,delimiter=',',dtype=dtype)
    
    dates = np.datetime64(met_data['time']).astype('datetime64[D]')

    met_days = np.unique(dates)
    if met_days[0]>start:
        print "\tWARNING: time series starts after beginning full period of interest!"

    if met_days[-1]<end:
        print "\tWARNING: time series finishes before full period of interest!"

    if met_days[0]<start:
        print "\tWARNING: time series starts prior to full period of interest, so clipping"
        met_data = met_data[dates>=start]
        dates = np.datetime64(met_data['time']).astype('datetime64[D]')
        met_days = np.unique(dates)

    if met_days[-1]<end:
        print "\tWARNING: time series continues beyond full period of interest, so clipping"
        met_data = met_data[dates<=end]
        dates = np.datetime64(met_data['time']).astype('datetime64[D]')
        met_days = np.unique(dates)

    n_days=met_days.size

    daily_met_data = {}
    daily_met_data['date'] = met_days.copy()
    # temperature - get minimum and maximum temperatures
    daily_met_data['minT'] = np.min(met_data_clipped['airt'].reshape(n_days,48),axis=1)
    daily_met_data['maxT'] = np.max(met_data_clipped['airt'].reshape(n_days,48),axis=1)
    # precipitation - get total daily precipitation
    daily_met_data['pptn'] = np.sum(met_data_clipped['pptn'].reshape(n_days,48),axis=1)
    # radiation - get total daily radiation
    daily_met_data['ssrd'] = np.sum(met_data_clipped['swr'].reshape(n_days,48),axis=1)*10.**6 # convert J m2 d-1 to MJ m-2 d-1
    # vpd - get average vpd out
    daily_met_data['vpd'] = np.mean(met_data_clipped['vpd'].reshape(n_days,48),axis=1)/10.**3 # convert Pa to kPa

    return met_days,daily_met_data['minT'], daily_met_data['maxT'], daily_met_data['vpd'], daily_met_data['ssrd'],daily_met_data['pptn']


"""
# This is an old function that was used to read in and gapfill data directly
def generate_daily_met_drivers_from_existing_halfhourly_time_series(met_data):
    earliest_day = np.min(met_data['date'])
    latest_day = np.max(met_data['date'])

    # clip time series so that ends are trimmed to full days
    start_fullday=met_data['date'][0].astype('datetime64[D]')
    if met_data['date'][0]- met_data['date'][0].astype('datetime64[D]').astype('datetime64[m]') != np.timedelta64(0,'m'):
        start_fullday+=np.timedelta64(1,'D')
    
    end_fullday = met_data['date'][-1].astype('datetime64[D]')+np.timedelta64(1,'D')
    if met_data['date'][-1]- met_data['date'][-1].astype('datetime64[D]').astype('datetime64[m]') != np.timedelta64(23*60+30,'m'):
        end_fullday -= np.timedelta64(1,'D')

    mask = np.all((met_data['date']>=start_fullday.astype('datetime64[m]'),met_data['date']<end_fullday.astype('datetime64[m]')),axis=0)

    met_data_clipped = {}
    var=met_data.keys()
    for vv in range(0,len(var)):
        if var[vv]=='swr':
            met_data_clipped[var[vv]]=met_data[var[vv]][mask]*30.*60. # convert W/m2 to J/m2
            
        else:
            met_data_clipped[var[vv]]=met_data[var[vv]][mask]
        
    met_days = np.arange(start_fullday,end_fullday,dtype='datetime64[D]')
    n_days=met_days.size

    daily_met_data = {}
    daily_met_data['date'] = met_days
    # temperature - get minimum and maximum temperatures
    daily_met_data['minT'] = np.min(met_data_clipped['airT'].reshape(n_days,48),axis=1)
    daily_met_data['maxT'] = np.max(met_data_clipped['airT'].reshape(n_days,48),axis=1)
    # precipitation - get total daily precipitation
    daily_met_data['pptn'] = np.sum(met_data_clipped['pptn'].reshape(n_days,48),axis=1)
    # radiation - get total daily radiation
    daily_met_data['ssrd'] = np.sum(met_data_clipped['swr'].reshape(n_days,48),axis=1)
    # vpd - get average vpd out
    daily_met_data['vpd'] = np.mean(met_data_clipped['vpd'].reshape(n_days,48),axis=1)

    return met_days,daily_met_data['minT'], daily_met_data['maxT'], daily_met_data['vpd'], daily_met_data['ssrd'],daily_met_data['pptn']
"""

# rolling mean met data based on specified number of previous days (default 21). 
# edge effects dealt with by assuming fixed boundaries
def retro_looking_rolling_mean(array,window_width):
    convolve_window = np.concatenate((np.zeros(window_width-1),np.ones(window_width)),axis=0)/float(window_width)
    host = np.zeros(array.size+(window_width-1)*2)
    host[window_width-1:window_width-1+array.size]=array.copy()
    host[:window_width-1]=array[0]
    host[-window_width+1:]=array[-1]
    array_out = np.convolve(host,convolve_window,'valid')
    return array_out


# recast met data as 21 day moving averages
def retro_rolling_average_met_data(mn2t,mx2t,ssrd,vpd,pptn, t=21):
    mn2t_out = retro_looking_rolling_mean(mn2t,t)
    mx2t_out = retro_looking_rolling_mean(mx2t,t)
    ssrd_out = retro_looking_rolling_mean(ssrd,t)
    pptn_out = retro_looking_rolling_mean(pptn,t)
    vpd_out = retro_looking_rolling_mean(vpd,t)
    return mn2t_out,mx2t_out,ssrd_out, vpd_out, pptn_out
