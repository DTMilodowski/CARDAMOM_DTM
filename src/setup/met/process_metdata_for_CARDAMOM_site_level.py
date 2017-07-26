# This code sets up the met drivers for CARDAMOM
#
# CARDAMOM-DALEC is driven with remotely sensed meteorologoical data - from ERA-Interim Reanalysis and TRMM.
# At the plot level, I do not generally include fire and/or deforestation unless directly observed in the
# field campaigns - this is not observed at BALI sites

# import some libraries - update as required
import numpy as np
import sys
import met_setup as met

#---------------------------------------------------------------------------------------------------------------
# First of all, define the relevant input files
# MET DATA
#    -ERA Interim
ERA_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/ERAinterim/BALI_Met/BALI_ERA_Interim_Met.txt'
#    -TRMM
TRMM_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/TRMM/g4.areaAvgTimeSeries.TRMM_3B42_007_precipitation.20110101-20160429.117E_4N_117E_4N.csv'
# Gapfilled met data
met_file = '/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/construct_drivers/BALI_gapfilled_met_station_30mins_v2.csv'

# Second define output drivers. In this instance I am going to write drivers for both the EO data only, and gapfilled station
# data
# MET DATA
outdir = "/exports/csce/datastore/geos/users/dmilodow/BALI/CARDAMOM_BALI/met_data/"
outfile_EO = "BALI_ERAinterim_TRMM_daily_v1.csv"
outfile_station = "BALI_metstation_daily_v1.csv"

#---------------------------------------------------------------------------------------------------------------
# Now get some basic parameters for the run
start_date= '01/01/2011'
end_date= '01/03/2016'

# Initiate some arrays to host time series
d,m,y = start_date.split('/')
start = np.datetime64(y+'-'+m+'-'+d,'D')
d,m,y = end_date.split('/')
end = np.datetime64(y+'-'+m+'-'+d,'D')
date = np.arange(start,end+np.timedelta64(1,'D'), dtype = 'datetime64[D]')

N_t = date.size

mn2t_in = np.zeros(N_t)-9999.
mx2t_in = np.zeros(N_t)-9999.
vpd_in = np.zeros(N_t)-9999.
ssrd_in = np.zeros(N_t)-9999.
pptn_in = np.zeros(N_t)-9999.
mn2t21_in = np.zeros(N_t)-9999.
mx2t21_in = np.zeros(N_t)-9999.
vpd21_in = np.zeros(N_t)-9999.
ssrd21_in = np.zeros(N_t)-9999.
pptn21_in = np.zeros(N_t)-9999.

#---------------------------------------------------------------------------------------------------------------
# Process EO met data
#   -this function returns a dictionary with time series of meteorological variables to be assimilated into
#    CARDAMOM
#   - Variable keys: Time, airT, pptn, vpd, par, swr, sp
met_dates, mn2t, mx2t, vpd, ssrd, TRMM_dates, pptn= met.generate_daily_met_drivers_ERAinterim_TRMM(ERA_file, TRMM_file, start, end)
N_m = met_dates.size
for dd in range(0,N_m):
    mn2t_in[date == met_dates[dd]] = mn2t[dd]
    mx2t_in[date == met_dates[dd]] = mx2t[dd]
    ssrd_in[date == met_dates[dd]] = ssrd[dd]
    vpd_in[date == met_dates[dd]] = vpd[dd]
N_trmm = TRMM_dates.size
for dd in range(0,N_trmm):
    pptn_in[date == TRMM_dates[dd]] = pptn[dd]

# also get 21 day rolling average (retro-looking)
mn2t21,mx2t21,ssrd21,vpd21,pptn21 = met.retro_rolling_average_met_data(mn2t,mx2t,ssrd,vpd,pptn)
for dd in range(0,N_m):
    mn2t21_in[date == met_dates[dd]] = mn2t21[dd]
    mx2t21_in[date == met_dates[dd]] = mx2t21[dd]
    ssrd21_in[date == met_dates[dd]] = ssrd21[dd]
    vpd21_in[date == met_dates[dd]] = vpd21[dd]
for dd in range(0,N_trmm):
    pptn21_in[date == TRMM_dates[dd]] = pptn21[dd]

    mx2t_in[np.isnan(mx2t_in)]=-9999.
    mn2t_in[np.isnan(mn2t_in)]=-9999.
    ssrd_in[np.isnan(ssrd_in)]=-9999.
    vpd_in[np.isnan(vpd_in)]=-9999.
    pptn_in[np.isnan(pptn_in)]=-9999.
    mx2t21_in[np.isnan(mx2t_in)]=-9999.
    mn2t21_in[np.isnan(mn2t_in)]=-9999.
    ssrd21_in[np.isnan(ssrd_in)]=-9999.
    vpd21_in[np.isnan(vpd_in)]=-9999.
    pptn21_in[np.isnan(pptn_in)]=-9999.

# write met data to file
out_drivers = open(outdir+outfile_EO,'w')
out_drivers.write('timestep_days, date, mn2t, mx2t, vpd, ssrd, pptn, mn2t_21d, mx2t_21d, vpd_21d, ssrd_21d, pptn_21d\n')
for tt in range(0,N_t):
    out_drivers.write(str(tt) + ',' + str(date[tt]) + ', ' + str(mn2t_in[tt]) + ',' + str(mx2t_in[tt]) + ',' + str(vpd_in[tt]) + ',' + str(ssrd_in[tt]) + ',' + str(pptn_in[tt]) + ',' + str(mn2t21_in[tt]) + ',' + str(mx2t21_in[tt]) + ',' + str(vpd21_in[tt]) + ',' + str(ssrd21_in[tt]) + ',' + str(pptn21_in[tt]) + '\n')


#---------------------------------------------------------------------------------------------------------------
# Process station met data
#   -this function returns a dictionary with time series of meteorological variables to be assimilated into
#    CARDAMOM
#   - Variable keys: Time, airT, pptn, vpd, par, swr, sp
met_dates, mn2t, mx2t, vpd, ssrd, pptn = met.generate_daily_met_drivers_from_existing_halfhourly_time_series(met_file, start, end)

mn2t_in = np.zeros(N_t)-9999.
mx2t_in = np.zeros(N_t)-9999.
vpd_in = np.zeros(N_t)-9999.
ssrd_in = np.zeros(N_t)-9999.
pptn_in = np.zeros(N_t)-9999.
mn2t21_in = np.zeros(N_t)-9999.
mx2t21_in = np.zeros(N_t)-9999.
vpd21_in = np.zeros(N_t)-9999.
ssrd21_in = np.zeros(N_t)-9999.
pptn21_in = np.zeros(N_t)-9999.

#--------------------------------
N_m = met_dates.size
for dd in range(0,N_m):
    mn2t_in[date == met_dates[dd]] = mn2t[dd]
    mx2t_in[date == met_dates[dd]] = mx2t[dd]
    ssrd_in[date == met_dates[dd]] = ssrd[dd]
    vpd_in[date == met_dates[dd]] = vpd[dd]
    pptn_in[date == met_dates[dd]] = pptn[dd]

# also get 21 day rolling average (retro-looking)
mn2t21,mx2t21,ssrd21,vpd21,pptn21 = met.retro_rolling_average_met_data(mn2t,mx2t,ssrd,vpd,pptn)
for dd in range(0,N_m):
    mn2t21_in[date == met_dates[dd]] = mn2t21[dd]
    mx2t21_in[date == met_dates[dd]] = mx2t21[dd]
    ssrd21_in[date == met_dates[dd]] = ssrd21[dd]
    vpd21_in[date == met_dates[dd]] = vpd21[dd]
    pptn21_in[date == met_dates[dd]] = pptn21[dd]

    
    mx2t_in[np.isnan(mx2t_in)]=-9999.
    mn2t_in[np.isnan(mn2t_in)]=-9999.
    ssrd_in[np.isnan(ssrd_in)]=-9999.
    vpd_in[np.isnan(vpd_in)]=-9999.
    pptn_in[np.isnan(pptn_in)]=-9999.
    mx2t21_in[np.isnan(mx2t_in)]=-9999.
    mn2t21_in[np.isnan(mn2t_in)]=-9999.
    ssrd21_in[np.isnan(ssrd_in)]=-9999.
    vpd21_in[np.isnan(vpd_in)]=-9999.
    pptn21_in[np.isnan(pptn_in)]=-9999.

# write met data to file
outfile_drivers = "BALI_gapfilled_met_station_daily_v1.csv"
out_drivers = open(outfile_drivers,'w')
out_drivers.write('timestep_days, date, mn2t, mx2t, vpd, ssrd, pptn, mn2t_21d, mx2t_21d, vpd_21d, ssrd_21d, pptn_21d\n')
for tt in range(0,N_t):
    out_drivers.write(str(tt) + ',' + str(date[tt]) + ', ' + str(mn2t_in[tt]) + ',' + str(mx2t_in[tt]) + ',' + str(vpd_in[tt]) + ',' + str(ssrd_in[tt]) + ',' + str(pptn_in[tt]) + ',' + str(mn2t21_in[tt]) + ',' + str(mx2t21_in[tt]) + ',' + str(vpd21_in[tt]) + ',' + str(ssrd21_in[tt]) + ',' + str(pptn21_in[tt]) + '\n')

out_drivers.close()
