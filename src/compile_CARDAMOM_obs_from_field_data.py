import numpy as np
import sys
sys.path.append('/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/CARDAMOM/CARDAMOM_setup_drivers/')
import met_setup as met
import MODIS_setup as MODIS
import field_data_setup as field
sys.path.append('/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/construct_drivers/')
import gapfill_station_metdata as gap
# First of all, define the relevant input files
# MET DATA
#    -ERA Interim
ERA_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/ERAinterim/BALI_Met/BALI_ERA_Interim_Met.txt'
#    -TRMM
TRMM_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/TRMM/g4.areaAvgTimeSeries.TRMM_3B42_007_precipitation.20110101-20160429.117E_4N_117E_4N.csv'

# Gapfilled met data
met_file = '/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/construct_drivers/BALI_gapfilled_met_station_daily_v1.csv'

# FIELD DATA
#    -Plot census data
#    -Fine roots data
#    -Litter data
#    -LiDAR data (for LAI)
census_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_TreeCensus.csv'
roots_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_FineRoots_Stock_NPP_RawData.csv'
litter_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_Litterfall_RawData.csv'
LAI_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_LAI_fromHemisphericalPhotos_TimeSeries.csv'

#---------------------------------------------------------------------------------------------------------------
# Now get some basic parameters for the run
start_date= '01/01/2011'
end_date= '01/03/2016'
plot = ['Belian','LF','B North','B South', 'E', 'Seraya', 'DC1', 'DC2']
LAI_MH = [6.69,4.78,3.00,2.26,3.84,6.22,5.93,5.89]
LAI_rad = [8.30,5.76,4.87,3.73,5.70,9.01,8.25,9.35]
LAI_hemiphot = [4.46,3.76,3.65,3.44,3.93,4.27,4.40,4.05]
Csoil = [8295.66, 11275.18, 3934.03, 4916.91, 11925.08, 24347.79, 8144.94, -9999.]

# 


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
# Process met data
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
outfile_drivers = "BALI_ERAinterim_TRMM_daily_v1.csv"
out_drivers = open(outfile_drivers,'w')
out_drivers.write('timestep_days, date, mn2t, mx2t, vpd, ssrd, pptn, mn2t_21d, mx2t_21d, vpd_21d, ssrd_21d, pptn_21d\n')
for tt in range(0,N_t):
    out_drivers.write(str(tt) + ',' + str(date[tt]) + ', ' + str(mn2t_in[tt]) + ',' + str(mx2t_in[tt]) + ',' + str(vpd_in[tt]) + ',' + str(ssrd_in[tt]) + ',' + str(pptn_in[tt]) + ',' + str(mn2t21_in[tt]) + ',' + str(mx2t21_in[tt]) + ',' + str(vpd21_in[tt]) + ',' + str(ssrd21_in[tt]) + ',' + str(pptn21_in[tt]) + '\n')

#---------------------------------------------------------------------------------------------
# now produce an equivalent met data file that utilises station data where possible.
# Met station
met_file  = '/home/dmilodow/DataStore_DTM/BALI/SAFE_data/SAFE_FluxTower_AtmMet_data.csv'
soil_file  = '/home/dmilodow/DataStore_DTM/BALI/SAFE_data/SAFE_FluxTower_SoilMet_data.csv'
start_date= '01/01/2011 00:00'
end_date= '01/03/2016 00:00'
met_data_dict, soil_data_dict, RS_data_dict = gap.load_all_metdata(met_file, soil_file, ERA_file, TRMM_file, start_date, end_date)
# remove swr and PAR record from station prior to 22/09/2012 as the sensor was behaving oddly
mask = met_data_dict['date']<np.datetime64('2012-09-22 00:00','m')
met_data_dict['PAR'][mask]=np.nan
met_data_dict['swr'][mask]=np.nan 
minimum_pptn_rate = 0.5
STA_LTA_threshold = 4.
gaps = gap.locate_metdata_gaps_using_soil_moisture_time_series(met_data_dict, soil_data_dict, minimum_pptn_rate, STA_LTA_threshold)
gapfilled_met_data = gap.gapfill_metdata(met_data_dict,RS_data_dict,gaps)

met_dates, mn2t, mx2t, vpd, ssrd, pptn = met.generate_daily_met_drivers_from_existing_halfhourly_time_series(gapfilled_met_data)
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


# Now deal with the obs
for pp in range(0,len(plot)):
    print plot[pp]
    Cwood_in = np.zeros(N_t)*np.nan#-9999.
    Croot_in = np.zeros(N_t)*np.nan#-9999.
    Croot_std_in = np.zeros(N_t)*np.nan#-9999.
    Croot_npp_in = np.zeros(N_t)*np.nan
    Croot_npp_std_in = np.zeros(N_t)*np.nan
    Csoil_in = np.zeros(N_t)*np.nan#-9999.
    Litter_in = np.zeros(N_t)*np.nan#-9999.
    Litter_std_in = np.zeros(N_t)*np.nan#-9999.
    LAI_in = np.zeros(N_t)*np.nan#-9999.
    LAI_std_in = np.zeros(N_t)*np.nan#-9999.
    LAI_MH_in = np.zeros(N_t)*np.nan#-9999.
    LAI_rad_in = np.zeros(N_t)*np.nan#-9999.
    LAI_rad_std_in = np.zeros(N_t)*np.nan#-9999.
    LAI_MH_std_in = np.zeros(N_t)*np.nan#-9999.

    # LAI data
    LAI_date, LAI, LAI_std = field.get_LAI_ts(LAI_file,plot[pp])
    N_LAI = LAI_date.size
    for tt in range(0,N_LAI):
        LAI_in[date==LAI_date[tt]] = LAI[tt]
        LAI_rad_in[date==LAI_date[tt]] = LAI[tt]*LAI_rad[pp]/LAI_hemiphot[pp]#0.025*LAI[tt]**3.90
        LAI_MH_in[date==LAI_date[tt]] = LAI[tt]*LAI_MH[pp]/LAI_hemiphot[pp]#0.084*LAI[tt]**3.26
        LAI_std_in[date==LAI_date[tt]] = LAI_std[tt]
        LAI_rad_std_in[date==LAI_date[tt]] = LAI_std[tt]*LAI_rad[pp]/LAI_hemiphot[pp]#0.025*LAI_std[tt]**3.90
        LAI_MH_std_in[date==LAI_date[tt]] = LAI_std[tt]*LAI_MH[pp]/LAI_hemiphot[pp]#0.084*LAI_std[tt]**3.26

    # soil carbon reported in g/m2
    Csoil_in[0] = Csoil[pp]

    # Cwood reported in kg (for 1ha plot)
    census_date, Cwood = field.get_Cwood_ts(census_file,plot[pp])
    N_c = Cwood.size
    for dd in range(0,N_c):
        Cwood_in[date == census_date[dd]] = Cwood[dd] *1000./10.**4. # convert kg/ha to g/m^2
    # root stocks reported in Mg/ha
    root_stock_date, Croot, Croot_std = field.get_Croot(roots_file,plot[pp])
    root_npp_date,root_npp_previous_date,root_npp,root_npp_std = field.get_root_NPP_ts(roots_file,plot[pp])
    Croot_in[date == root_stock_date] = Croot*10**6/10.**4 # convert Mg/ha to g/m^2
    Croot_std_in[date == root_stock_date] = Croot_std*10**6/10.**4 # convert Mg/ha to g/m^2

    N_roo=root_npp.size
    for tt in range(0,N_roo):
        indices = np.all((date>=root_npp_previous_date[tt], date<root_npp_date[tt]),axis=0)
        root_npp_in[indices]= root_npp[tt] * (10.**6/10.**4/365.25) # convert Mg/ha/yr to g/m2/d
        root_npp_std_in[indices]= root_npp_std[tt] * (10.**6/10.**4/365.25) # convert Mg/ha/yr to g/m2/d

    # litter fluxes reported in Mg/ha/yr
    litter_collection_date, litter_previous_collection_date, litter_flux, litter_std = field.get_litterfall_ts(litter_file,plot[pp])

    # Initially assume average flux rates for litter between collection dates 
    N_lit=litter_flux.size
    for tt in range(0,N_lit):
        indices = np.all((date>=litter_previous_collection_date[tt], date<litter_collection_date[tt]),axis=0)
        Litter_in[indices]= litter_flux[tt] * (10.**6/10.**4/365.25) # convert Mg/ha/yr to g/m2/d
        Litter_std_in[indices]= litter_std[tt] * (10.**6/10.**4/365.25) # convert Mg/ha/yr to g/m2/d
        
    # Convert nodata to -9999
    Litter_in[np.isnan(Litter_in)]=-9999.
    Litter_std_in[np.isnan(Litter_std_in)]=-9999.
    root_npp_in[np.isnan(Litter_in)]=-9999.
    root_npp_std_in[np.isnan(Litter_std_in)]=-9999.
    LAI_MH_in[np.isnan(LAI_MH_in)]=-9999.
    LAI_MH_std_in[np.isnan(LAI_MH_std_in)]=-9999.
    LAI_rad_in[np.isnan(LAI_rad_in)]=-9999.
    LAI_rad_std_in[np.isnan(LAI_rad_std_in)]=-9999.
    LAI_in[np.isnan(LAI_in)]=-9999.
    LAI_std_in[np.isnan(LAI_std_in)]=-9999.
    Croot_in[np.isnan(Croot_in)]=-9999.
    Croot_std_in[np.isnan(Croot_in)]=-9999.
    Cwood_in[np.isnan(Cwood_in)]=-9999.
    Csoil_in[np.isnan(Csoil_in)]=-9999.



    # build an output matrix for the observations file
    obs_h = ['GPP','GPP_u','LAI','LAI_u','NEE','NEE_u','woo','woo_u','Reco','Reco_u','Cfol','Cfol_u','Cwoo','Cwoo_u','Croo','Croo_u','Clit','Clit_u','Csom','Csom_u','Cagb','Cagb_u','Cstem','Cstem_u','Cbranch','Cbranch_u','Ccroo','Ccroo_u','Cfol_max','Cfol_max_u','Evap','Evap_u','flit','flit_u','NPProo','NPProo_u']
    obs = np.zeros((len(obs_h),N_t))*-9999.
    
    # fill obs matrix with relevant data
    obs[obs_h=='Cwoo',:]=Cwood_in.copy()
    obs[obs_h=='Cwoo_u',:]=Cwood_in*0.20 # for now, assume 20% error on Cwood
    obs[obs_h=='Croo',:]=Croot_in.copy()
    obs[obs_h=='Croo_u',:]=Croot_std_in.copy()
    obs[obs_h=='flit',:]=Litter_in.copy()
    obs[obs_h=='flit_u',:]=Litter_std_in.copy()
    obs[obs_h=='NPProo',:]=root_npp_in.copy()
    obs[obs_h=='NPProo_u',:]=root_npp_std_in.copy()
    obs[obs_h=='LAI',:]=LAI_rad_in.copy()
    obs[obs_h=='LAI_u',:]=LAI_rad_std_in.copy()
    obs[obs_h=='Csom',:]=Csoil_in.copy()
    obs[obs_h=='Csom_u',:]=Csoil_in*0.50 # for now, assume 50% error on Csoil as only have one pit

    # write output to file
    outfile_obs = "CARDAMOM_obs_"+plot[pp]+".csv"
    out_obs = open(outfile_obs,'w')
    obs.write(obs_h[0])
    for vv in range(1,len(obs_h)):
        obs.write(','+obs_h[vv])
    obs.write(obs_h[vv]+'\n')
    # now write in data
    for tt in range(0,N_t):
        out_priors.write(str(tt) + ',' + str(date[tt]))
        for vv in range(0,len(obs_h)):
            obs.write(', ' + str(obs[vv,tt]))
            obs.write('\n')
    out_priors.close()
