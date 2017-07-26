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
end_date= '31/12/2015'
plot = ['Belian','LF','B North','B South', 'E', 'Seraya', 'DC1', 'DC2']
LAI_MH = [6.69,4.78,3.00,2.26,3.84,6.22,5.93,5.89]
LAI_rad = [8.30,5.76,4.87,3.73,5.70,9.01,8.25,9.35]
LAI_hemiphot = [4.46,3.76,3.65,3.44,3.93,4.27,4.40,4.05]
Csoil = [8295.66, 11275.18, 3934.03, 4916.91, 11925.08, 24347.79, 8144.94, -9999.]

# Initiate some arrays to host time series
d,m,y = start_date.split('/')
start = np.datetime64(y+'-'+m+'-'+d,'D')
d,m,y = end_date.split('/')
end = np.datetime64(y+'-'+m+'-'+d,'D')
date = np.arange(start,end+np.timedelta64(1,'D'), dtype = 'datetime64[D]')

N_t = date.size

# Now deal with the obs
for pp in range(0,len(plot)):
    print plot[pp]
    Cwood_in = np.zeros(N_t)*np.nan#-9999.
    Croot_in = np.zeros(N_t)*np.nan#-9999.
    Croot_std_in = np.zeros(N_t)*np.nan#-9999.
    root_npp_in = np.zeros(N_t)*np.nan
    root_npp_std_in = np.zeros(N_t)*np.nan
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
    root_npp_in[np.isnan(root_npp_in)]=-9999.
    root_npp_std_in[np.isnan(root_npp_std_in)]=-9999.
    LAI_MH_in[np.isnan(LAI_MH_in)]=-9999.
    LAI_MH_std_in[np.isnan(LAI_MH_std_in)]=-9999.
    LAI_rad_in[np.isnan(LAI_rad_in)]=-9999.
    LAI_rad_std_in[np.isnan(LAI_rad_std_in)]=-9999.
    LAI_in[np.isnan(LAI_in)]=-9999.
    LAI_std_in[np.isnan(LAI_std_in)]=-9999.
    Croot_in[np.isnan(Croot_in)]=-9999.
    Croot_std_in[np.isnan(Croot_std_in)]=-9999.

    Cwood_unc=Cwood_in*0.20 # for now, assume 20% error on Cwood
    Cwood_in[np.isnan(Cwood_in)]=-9999.
    Cwood_unc[np.isnan(Cwood_unc)]=-9999.

    Csoil_unc=Csoil_in*0.50 # for now, assume 50% error on Csoil as only have one pit
    Csoil_in[np.isnan(Csoil_in)]=-9999.
    Csoil_unc[np.isnan(Csoil_unc)]=-9999.

    # build an output matrix for the observations file
    obs_h = np.asarray(['GPP','GPP_u','LAI','LAI_u','NEE','NEE_u','woo','woo_u','Reco','Reco_u','Cfol','Cfol_u','Cwoo','Cwoo_u','Croo','Croo_u','Clit','Clit_u','Csom','Csom_u','Cagb','Cagb_u','Cstem','Cstem_u','Cbranch','Cbranch_u','Ccroo','Ccroo_u','Cfol_max','Cfol_max_u','Evap','Evap_u','flit','flit_u','NPProo','NPProo_u'])
    obs = np.zeros((len(obs_h),N_t))-9999.
    
    # fill obs matrix with relevant data
    obs[obs_h=='Cwoo',:]=Cwood_in.copy()
    obs[obs_h=='Cwoo_u',:]=Cwood_unc.copy()
    obs[obs_h=='Croo',:]=Croot_in.copy()
    obs[obs_h=='Croo_u',:]=Croot_std_in.copy()
    obs[obs_h=='flit',:]=Litter_in.copy()
    obs[obs_h=='flit_u',:]=Litter_std_in.copy()
    obs[obs_h=='NPProo',:]=root_npp_in.copy()
    obs[obs_h=='NPProo_u',:]=root_npp_std_in.copy()
    obs[obs_h=='LAI',:]=LAI_rad_in.copy()
    #obs[obs_h=='LAI_u',:]=LAI_rad_std_in.copy()
    obs[obs_h=='Csom',:]=Csoil_in.copy()
    obs[obs_h=='Csom_u',:]=Csoil_unc.copy()

    # write output to file
    outfile_obs = "CARDAMOM_obs_"+plot[pp]+".csv"
    out_obs = open(outfile_obs,'w')
    out_obs.write('tstep_days,date')
    for vv in range(0,len(obs_h)):
        out_obs.write(','+obs_h[vv])
    out_obs.write('\n')
    # now write in data
    for tt in range(0,N_t):
        out_obs.write(str(tt) + ',' + str(date[tt]))
        for vv in range(0,len(obs_h)):
            out_obs.write(', ' + str(obs[vv,tt]))
        out_obs.write('\n')
    out_obs.close()
