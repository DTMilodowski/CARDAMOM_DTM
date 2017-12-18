import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

def read_ICP_census_data(census_file):
    
    datatype = {'names': ('ForestType', 'Plot', 'Subplot', 'CensusDate1', 'Observers1', 'TagNumber1', 'D_POM1', 'H_POM1', 'Height1', 'RAINFOR_flag1', 'Alive_flag1', 'C_stem1','C_coarse_root1', 'Comment1', 'CommentData1', 'CensusDate2', 'Observers2', 'TagNumber2', 'D_POM2', 'H_POM2', 'Height2', 'RAINFOR_flag2', 'Alive_flag2', 'C_stem2','C_coarse_root2', 'Comment2', 'CommentData2', 'CensusDate3', 'Observers3', 'TagNumber3', 'D_POM3', 'H_POM3', 'Height3', 'RAINFOR_flag3', 'Alive_flag3', 'C_stem3','C_coarse_root3', 'Comment3', 'CommentData3', 'SubplotX', 'SubplotY', 'SpeciesID', 'WoodDensity','Source','Quality'),'formats': ('S16','S16','i8','S10','S64','f8','f16','f16','f16','S8','i8','f32','f32','S16','S64','S10','S64','f8','f16','f16','f16','S8','i8','f32','f32','S16','S64','S10','S64','f8','f16','f16','f16','S8','i8','f32','f32','S16','S64','i8','i8','S32','f16','S100','S32')}
    census_data = np.genfromtxt(census_file, skiprows = 1, delimiter = ',',dtype=datatype)
    
    nrows = census_data['Plot'].size
    BALI_plot = []
    Subplot = np.zeros(nrows)
    CensusDates = np.zeros((nrows,3),dtype = 'datetime64[D]')
    TreeTag = np.zeros(nrows)
    AltTag = np.zeros(nrows)
    DPOM = np.zeros((nrows,3))
    HPOM = np.zeros((nrows,3))
    Height = np.zeros((nrows,3))
    C_stem = np.zeros((nrows,3))
    C_coarse_root = np.zeros((nrows,3))
    RAINFOR_flag = []
    Alive_flag = np.ones((nrows,3)) # assume alive unless specified
    Species = []
    SubplotCoords = np.ones((nrows,2))
    WoodDensity = np.zeros(nrows)
    AltTag[:] = np.nan

    for i in range(0,nrows):

        BALI_plot.append(census_data['Plot'][i])
        Subplot[i] = census_data['Subplot'][i]
        Species.append(census_data['SpeciesID'][i])
        SubplotCoords[i,0] = census_data['SubplotX'][i]
        SubplotCoords[i,1] = census_data['SubplotY'][i]
        WoodDensity[i] = census_data['WoodDensity'][i]
        RAINFOR_flag.append([])

        # Census #1
        prior_census = 0
        if np.isfinite(census_data['TagNumber1'][i]):
            prior_census = 1
            if len(census_data['CensusDate1'][i])<3:
                CensusDates[i,0] = np.datetime64(0,'D')#np.datetime64('1990-01-01')
            elif census_data['CensusDate1'][i][2]=="/":
                day,month,year = census_data['CensusDate1'][i].split("/")
                CensusDates[i,0] = np.datetime64(year+'-'+month+'-'+day)
            elif census_data['CensusDate1'][i][2]=="-":
                month_string = np.asarray(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
                month_num = np.asarray(['01','02','03','04','05','06','07','08','09','10','11','12'])
                day,month,year = census_data['CensusDate1'][i].split("-")
                year_str = str(2000+int(year))
                CensusDates[i,0] = np.datetime64(year_str+'-'+month_num[month_string==month][0]+'-'+day)
            else:
                print "row", i, "\t incorrect date format in CensusDate1"
                CensusDates[i,0] =np.datetime64(0,'D')# np.datetime64('1990-01-01')

            TreeTag[i] = census_data['TagNumber1'][i]
            DPOM[i,0] = census_data['D_POM1'][i]
            HPOM[i,0] = census_data['H_POM1'][i]
            Height[i,0] = census_data['Height1'][i]
            C_stem[i,0] = census_data['C_stem1'][i]
            C_coarse_root[i,0] = census_data['C_coarse_root1'][i]
            RAINFOR_flag[i].append(census_data['RAINFOR_flag1'][i])
            Alive_flag[i,0] = census_data['Alive_flag1'][i]
        else:
            CensusDates[i,0] = np.datetime64(0,'D')#np.datetime64('1990-01-01')
            TreeTag[i] = np.nan
            DPOM[i,0] = np.nan
            HPOM[i,0] = np.nan
            Height[i,0] = np.nan
            C_stem[i,0] = 0#np.nan
            C_coarse_root[i,0] = 0#np.nan
            RAINFOR_flag[i].append('NaN')
            Alive_flag[i,0] = np.nan

        # Subsequent censuses
        for j in range(1,3):
            if prior_census==0: # New recruits
                if np.isfinite(census_data['TagNumber'+str(j+1)][i]):
                    prior_census = 1
                    if len(census_data['CensusDate'+str(j+1)][i])<3:
                        CensusDates[i,j] = np.datetime64(0,'D')#np.datetime64('1990-01-01')
                    elif census_data['CensusDate'+str(j+1)][i][2]=="/":
                        day,month,year = census_data['CensusDate'+str(j+1)][i].split("/")
                        CensusDates[i,j] = np.datetime64(year+'-'+month+'-'+day)
                    elif census_data['CensusDate'+str(j+1)][i][2]=="-":
                        month_string = np.asarray(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
                        month_num = np.asarray(['01','02','03','04','05','06','07','08','09','10','11','12'])
                        day,month,year = census_data['CensusDate'+str(j+1)][i].split("-")
                        year_str = str(2000+int(year))
                        CensusDates[i,j] = np.datetime64(year_str+'-'+month_num[month_string==month][0]+'-'+day)
                    else:
                        print "row", i, "\t incorrect date format in CensusDate"+str(j+1)
                        CensusDates[i,j] =np.datetime64(0,'D')# np.datetime64('1990-01-01')

                    #day,month,year = census_data['CensusDate'+str(j+1)][i].split("/")
                    #CensusDates[i,j] = np.datetime64(year+'-'+month+'-'+day)
                    TreeTag[i] = census_data['TagNumber'+str(j+1)][i]
                    DPOM[i,j] = census_data['D_POM'+str(j+1)][i]
                    HPOM[i,j] = census_data['H_POM'+str(j+1)][i]
                    Height[i,j] = census_data['Height'+str(j+1)][i]
                    C_stem[i,j] = census_data['C_stem'+str(j+1)][i]
                    C_coarse_root[i,j] = census_data['C_coarse_root'+str(j+1)][i]
                    RAINFOR_flag[i].append(census_data['RAINFOR_flag'+str(j+1)][i])
                    Alive_flag[i,j] = census_data['Alive_flag'+str(j+1)][i]
                else:
                    CensusDates[i,j] = np.datetime64(0,'D')#np.datetime64('1990-01-01')
                    TreeTag[i] = np.nan
                    DPOM[i,j] = np.nan
                    HPOM[i,j] = np.nan
                    Height[i,j] = np.nan
                    C_stem[i,j] = 0
                    C_coarse_root[i,j] = 0
                    RAINFOR_flag[i].append('NaN')
                    Alive_flag[i,j] = np.nan

            else: # trees present in prior census
                if np.isfinite(census_data['TagNumber'+str(j+1)][i]):

                    if len(census_data['CensusDate'+str(j+1)][i])<3:
                        CensusDates[i,j] = np.datetime64(0,'D')#np.datetime64('1990-01-01')
                    elif census_data['CensusDate'+str(j+1)][i][2]=="/":
                        day,month,year = census_data['CensusDate'+str(j+1)][i].split("/")
                        CensusDates[i,j] = np.datetime64(year+'-'+month+'-'+day)
                    elif census_data['CensusDate'+str(j+1)][i][2]=="-":
                        month_string = np.asarray(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
                        month_num = np.asarray(['01','02','03','04','05','06','07','08','09','10','11','12'])
                        day,month,year = census_data['CensusDate'+str(j+1)][i].split("-")
                        year_str = str(2000+int(year))
                        CensusDates[i,j] = np.datetime64(year_str+'-'+month_num[month_string==month][0]+'-'+day)
                    else:
                        print "row", i, "\t incorrect date format in CensusDate"+str(j+1)
                        CensusDates[i,j] = np.datetime64(0,'D')#np.datetime64('1990-01-01')

                    DPOM[i,j] = census_data['D_POM'+str(j+1)][i]
                    HPOM[i,j] = census_data['H_POM'+str(j+1)][i]
                    Height[i,j] = census_data['Height'+str(j+1)][i]
                    C_stem[i,j] = census_data['C_stem'+str(j+1)][i]
                    C_coarse_root[i,j] = census_data['C_coarse_root'+str(j+1)][i]
                    RAINFOR_flag[i].append(census_data['RAINFOR_flag'+str(j+1)][i])
                    Alive_flag[i,j] = census_data['Alive_flag'+str(j+1)][i]
                    if census_data['TagNumber'+str(j+1)][i]!=TreeTag[i]:
                        AltTag[i]= census_data['TagNumber'+str(j+1)][i]
                else:
                    CensusDates[i,j] =np.datetime64(0,'D')# np.datetime64('1990-01-01')
                    #TreeTag[i] = np.nan
                    DPOM[i,j] = np.nan
                    HPOM[i,j] = np.nan
                    Height[i,j] = np.nan
                    C_stem[i,j] = 0#np.nan
                    C_coarse_root[i,j] = 0#np.nan
                    RAINFOR_flag[i].append('NaN')
                    Alive_flag[i,j] = np.nan
    
    plot_names=np.asarray(BALI_plot)
    RAINFOR = np.asarray(RAINFOR_flag)
    spp = np.asarray(Species)
    C_stem[np.isnan(C_stem)]=0
    C_coarse_root[np.isnan(C_coarse_root)]=0
    return plot_names, Subplot, CensusDates, TreeTag, AltTag, DPOM, HPOM, Height, C_stem, C_coarse_root, RAINFOR, Alive_flag, spp, SubplotCoords, WoodDensity

def collate_plot_level_census_data(census_file):
    BALI_plot, Subplot, CensusDates, TreeTag_census, AltTag, DPOM, HPOM, Height, C_stem, C_coarse_root, RAINFOR_flag, Alive_flag, Species, SubplotCoords, WoodDensity = read_ICP_census_data(census_file)

    plot_names = np.unique(BALI_plot)
    N_plots = plot_names.size
    
    CensusDict={}
    # Interesting properties: CanopyHeight, C_stem, C_coarse_roots, CensusDate
    for i in range(0,N_plots):
        plot_indices = BALI_plot==plot_names[i]
        #Set up arrays to save
        subplot_ids = np.unique(Subplot)
        n_subplots = subplot_ids.size
        dates = np.zeros((n_subplots,3),dtype = 'datetime64[D]')
        CanHt = np.zeros((n_subplots,3))
        Cstem = np.zeros((n_subplots,3))
        Croot = np.zeros((n_subplots,3))
        BasalArea = np.zeros((n_subplots,3))
        # note that for growth, mortality and recruitment, first year of census will be nan values because no there are no previous surveys!
        Growth = np.zeros((n_subplots,3))*np.nan
        Mortality = np.zeros((n_subplots,3))*np.nan
        Recruitment = np.zeros((n_subplots,3))*np.nan

        for s in range(0,n_subplots):
            subplot_indices = plot_indices * Subplot==subplot_ids[s]
            for y in range(0,3):
                datetemp = CensusDates[subplot_indices,y]
                dates[s,y]= np.max(datetemp)

                Cstemtemp = C_stem[subplot_indices,y]
                if np.isfinite(Cstemtemp).sum()>0:
                    Cstem[s,y]= np.sum(Cstemtemp[np.isfinite(Cstemtemp)])
                else:
                    Cstem[s,y]=np.nan

                Croottemp = C_coarse_root[subplot_indices,y]
                if np.isfinite(Croottemp).sum()>0:
                    Croot[s,y]= np.sum(Croottemp[np.isfinite(Croottemp)])  
                else:
                    Croot[s,y]= np.nan

                httemp = Height[subplot_indices,y]
                if np.isfinite(httemp).sum()>0:
                    CanHt[s,y]= np.mean(httemp[np.isfinite(httemp)])   
                else:
                    CanHt[s,y]=np.nan

                DBHtemp = DPOM[subplot_indices,y]
                if np.isfinite(DBHtemp).sum()>0:
                    BasalArea[s,y]= np.pi*np.sum((DBHtemp[np.isfinite(DBHtemp)]/2)**2)   
                else:
                    BasalArea[s,y]=np.nan
        
        # need to do catch for where there are no trees in subplot that has been surveyed!
        for s in range(0,n_subplots):
            for y in range(0,3):
                if dates[s,y]==np.datetime64('1970-01-01','D'):
                    if np.max(dates[:,y])>np.datetime64('1970-01-01','D'):
                        dates[s,y]=np.max(dates[:,y])
                        Croot[s,y]=0.
                        Cstem[s,y]=0.
        # now lets do the growth, mortality and recruitment
        for s in range(0,n_subplots):
            subplot_indices = plot_indices * Subplot==subplot_ids[s]
            Cwood_temp = C_stem[subplot_indices]+C_coarse_root[subplot_indices]            
            for y in range(1,3):
                """
                growth_indices = np.all((np.isfinite(Cwood_temp[:,y-1]),np.isfinite(Cwood_temp[:,y])),axis=0)
                recruit_indices = np.all((np.isfinite(Cwood_temp[:,y]),~np.isfinite(Cwood_temp[:,y-1])),axis=0)
                mortality_indices = np.all((np.isfinite(Cwood_temp[:,y-1]),~np.isfinite(Cwood_temp[:,y])),axis=0)
                """
                growth_indices = np.all((Cwood_temp[:,y-1]>0,Cwood_temp[:,y]>0),axis=0)
                recruit_indices = np.all((Cwood_temp[:,y]>0,Cwood_temp[:,y-1]==0),axis=0)
                mortality_indices = np.all((Cwood_temp[:,y-1]>0,Cwood_temp[:,y]==0),axis=0)
                if np.isfinite(Cwood_temp).sum()>0:
                    Growth[s,y] = np.sum(Cwood_temp[:,y][growth_indices]-Cwood_temp[:,y-1][growth_indices])
                    Recruitment[s,y] = np.sum(Cwood_temp[:,y][recruit_indices])
                    Mortality[s,y] = np.sum(Cwood_temp[:,y-1][mortality_indices])
                else:
                     Growth[s,y] = np.nan
                     Recruitment[s,y] = np.nan
                     Mortality[s,y] = np.nan

        PlotDict = {}
        PlotDict['n_subplots']=n_subplots
        PlotDict['CanopyHeight']=CanHt
        PlotDict['C_stem']=Cstem
        PlotDict['C_coarse_roots']=Croot
        PlotDict['C_wood']=Croot+Cstem
        PlotDict['CensusDate']=dates
        PlotDict['BasalArea']=BasalArea
        PlotDict['Growth']= Growth
        PlotDict['Recruitment']= Recruitment
        PlotDict['Mortality']=Mortality
        CensusDict[plot_names[i]]=PlotDict
    return CensusDict


# Read litterfall from GEM plot database file (converted to csv).  Mass collections (denoted with prefix m) are in g accumulated.  Otherwise given as flux in Mg C ha-1 yr-1
def read_litterfall_data(litter_file):
   
    datatype = {'names': ('ForestType', 'Plot', 'CollectionDate', 'PreviousCollectionDate', 'AccumulationDays', 'Trap', 'TrapSize', 'mLeaves', 'mTwigs', 'mFruit', 'mFlowers', 'mSeeds', 'mMisc','Comments','Leaves', 'Twigs', 'Fruit', 'Flowers', 'Seeds', 'Misc', 'Reproductive','Total'), 'formats': ('S16','S16','S10','S10','f16','i8','f16','f16','f16','f16','f16','f16','f16','S64','f16','f16','f16','f16','f16','f16','f16','f16')}
    litter_data = np.genfromtxt(litter_file, skiprows = 1, delimiter = ',',dtype=datatype)

    plot_names = np.unique(litter_data['Plot'])
    N_plots = plot_names.size

    LitterDict = {}
    
    for i in range(0,N_plots):
        plot_data = litter_data[litter_data['Plot']==plot_names[i]]
        N_sub = np.max(np.unique(plot_data['Trap']))
        PlotDict = {}
        CollectionDates = []
        PreviousDates = []
        AccumulationDays = []
        TrapSize = []
        mLeaves = []
        mTwigs = []
        mFruit = []
        mFlowers = []
        mSeeds = []
        mMisc = []
        rLeaves = []
        rTwigs = []
        rFruit = []
        rFlowers = []
        rSeeds = []
        rRepro = []
        rTotal = []
        rMisc = []

        for j in range(0,N_sub):
            subplot_data = plot_data[plot_data['Trap']==j+1]
            N_collections = subplot_data.size
            
            cDates = np.zeros(N_collections, dtype = 'datetime64[D]')
            pDates = np.zeros(N_collections, dtype = 'datetime64[D]')
            
            # mass collected by component - note that masses need converting from g per trap to g m-2
            mLeaves.append(subplot_data['mLeaves']/subplot_data['TrapSize'])
            mTwigs.append(subplot_data['mTwigs']/subplot_data['TrapSize'])
            mFruit.append(subplot_data['mFruit']/subplot_data['TrapSize'])
            mFlowers.append(subplot_data['mFlowers']/subplot_data['TrapSize'])
            mSeeds.append(subplot_data['mSeeds']/subplot_data['TrapSize'])
            mMisc.append(subplot_data['mMisc']/subplot_data['TrapSize'])
            # Rates of litter fall by component
            rLeaves.append(subplot_data['Leaves'])
            rTwigs.append(subplot_data['Twigs'])
            rFruit.append(subplot_data['Fruit'])
            rFlowers.append(subplot_data['Flowers'])
            rSeeds.append(subplot_data['Seeds'])
            rMisc.append(subplot_data['Misc'])
            rRepro.append(subplot_data['Reproductive'])
            rTotal.append(subplot_data['Total'])
            AccumulationDays.append(subplot_data['AccumulationDays'])
            for k in range(0, N_collections):
                if len(subplot_data['CollectionDate'][k])>0:
                    day,month,year = subplot_data['CollectionDate'][k].split("/")
                    cDates[k]= np.datetime64(year+'-'+month+'-'+day)
                else:
                    cDates[k]= np.datetime64(0,'D')#

                if len(subplot_data['PreviousCollectionDate'][k])>0:
                    day,month,year = subplot_data['PreviousCollectionDate'][k].split("/")
                    pDates[k]= np.datetime64(year+'-'+month+'-'+day)
                else:
                    pDates[k]= np.datetime64(0,'D')
                #accDays[k]= cDates[k]-pDates[k]

            #AccumulationDays.append(accDays)
            PreviousDates.append(pDates)
            CollectionDates.append(cDates)
            TrapSize.append(subplot_data['TrapSize'][0])
        
        # Load plot data into plot dictionary
        PlotDict['N_Subplots']=N_sub
        PlotDict['CollectionDate']=np.asarray(CollectionDates)
        PlotDict['PreviousCollectionDate']=np.asarray(PreviousDates)
        PlotDict['AccumulationDays']=np.asarray(AccumulationDays)
        PlotDict['mLeaves']=np.asarray(mLeaves)
        PlotDict['mTwigs']=np.asarray(mTwigs)
        PlotDict['mFruit']=np.asarray(mFruit)
        PlotDict['mFlowers']=np.asarray(mFlowers)
        PlotDict['mSeeds']=np.asarray(mSeeds)
        PlotDict['mMisc']=np.asarray(mMisc)
        PlotDict['mTotal']=PlotDict['mLeaves']+PlotDict['mTwigs']+PlotDict['mFruit']+PlotDict['mFlowers']+PlotDict['mSeeds']+PlotDict['mMisc']
        PlotDict['rLeaves']=np.asarray(rLeaves)
        PlotDict['rTwigs']=np.asarray(rTwigs)
        PlotDict['rFruit']=np.asarray(rFruit)
        PlotDict['rFlowers']=np.asarray(rFlowers)
        PlotDict['rSeeds']=np.asarray(rSeeds)
        PlotDict['rMisc']=np.asarray(rMisc)
        PlotDict['rReproductive']=np.asarray(rRepro)
        PlotDict['rTotal']=np.asarray(rTotal)
        PlotDict['TrapSize']=np.asarray(TrapSize)

        LitterDict[plot_names[i]]=PlotDict

    return LitterDict


def read_soil_respiration_data(soil_resp_file):
    datatype = {'names': ('ForestType', 'Plot', 'Date', 'Observer', 'Collar', 'SoilMoisture', 'SoilT', 'AirT', 'RespirationFlux', 'Remarks', 'QualityFlag'), 'formats': ('S16','S16','S10','S32','i8','f16','f16','f16','f16','S64','i8')}
    resp_data = np.genfromtxt(soil_resp_file, skiprows = 1, delimiter = ',',dtype=datatype)

    plot_names = np.unique(resp_data['Plot'])
    N_plots = plot_names.size

    SoilRespDict = {}
    
    for i in range(0,N_plots):
        plot_data = resp_data[resp_data['Plot']==plot_names[i]]
        N_collars = np.max(np.unique(resp_data['Collar']))
        PlotDict = {}
        Dates = []
        SoilMoisture = []
        SoilT = []
        AirT = []
        Flux = []
        Flag = []

        for j in range(0,N_collars):
            collar_data = plot_data[plot_data['Collar']==j+1]
            N_collections = collar_data.size
            cDates = np.zeros(N_collections, dtype = 'datetime64[D]')

            SoilMoisture.append(collar_data['SoilMoisture'])
            SoilT.append(collar_data['SoilT'])
            AirT.append(collar_data['AirT'])
            Flux.append(collar_data['RespirationFlux'])
            Flag.append(collar_data['QualityFlag'])
            
            for k in range(0, N_collections):
                day,month,year = collar_data['Date'][k].split("/")
                cDates[k]= np.datetime64(year+'-'+month+'-'+day)

            Dates.append(cDates)
        
        # Load plot data into plot dictionary
        PlotDict['N_Collars']=N_collars
        PlotDict['CollectionDate']=np.asarray(Dates)
        PlotDict['SoilMoisture']=np.asarray(SoilMoisture)
        PlotDict['SoilT']=np.asarray(SoilT)
        PlotDict['AirT']=np.asarray(AirT)
        PlotDict['RespirationFlux']=np.asarray(Flux)
        PlotDict['QualityFlag']=np.asarray(Flag)

        SoilRespDict[plot_names[i]]=PlotDict

    return SoilRespDict

def read_soil_stocks_and_npp(roots_file):
    datatype = {'names': ('ForestType', 'Plot', 'Core', 'DataType', 'CollectionDate', 'PreviousCollectionDate', 'AccumulationDays', 'CoarseRoots', 'FineRoots', 'Remarks'), 'formats': ('S16','S16','i8','S5','S10','S10','i8','f16','f16','S64')}
    roots_data = np.genfromtxt(roots_file, skiprows = 1, delimiter = ',',dtype=datatype)

    plot_names = np.unique(roots_data['Plot'])
    N_plots = plot_names.size

    RootStocksDict = {}
    RootNPPDict = {}
    
    for i in range(0,N_plots):
        plot_data = roots_data[roots_data['Plot']==plot_names[i]]
        NPP_data = plot_data[plot_data['DataType']=='NPP']
        Stocks_data = plot_data[plot_data['DataType']=='Stock']
        N_cores = np.max(np.unique(roots_data['Core']))
        StocksDict = {}
        NPPDict = {}

        # First do stocks
        cDates = np.unique(Stocks_data['CollectionDate'])
        N_stocks = cDates.size
        StocksDates = np.zeros(N_stocks,dtype='datetime64[D]')
        for dd in range(0,N_stocks):
            day,month,year = cDates[dd].split("/")
            StocksDates[dd]= np.datetime64(year+'-'+month+'-'+day)

        StocksDates = np.sort(StocksDates)

        FineRootsStocks = np.zeros((N_cores,N_stocks))*np.nan
        CoarseRootsStocks =np.zeros((N_cores,N_stocks))*np.nan

        for cc in range(0,N_cores):
            core_data = Stocks_data[Stocks_data['Core']==cc+1]
            N_collections = core_data.size

            for kk in range(0, N_collections):
                day,month,year = core_data['CollectionDate'][kk].split("/") 
                dd = (StocksDates == np.datetime64(year+'-'+month+'-'+day))
                CoarseRootsStocks[cc,dd] = core_data['CoarseRoots'][kk]
                FineRootsStocks[cc,dd] = core_data['FineRoots'][kk]

        StocksDict['N_Cores']=N_cores
        StocksDict['CollectionDate']=StocksDates
        StocksDict['CoarseRootStocks']=CoarseRootsStocks
        StocksDict['FineRootStocks']=FineRootsStocks
        StocksDict['Mask'] = np.isfinite(FineRootsStocks)


        # Next do NPP
        cDates = np.unique(NPP_data['CollectionDate'])
        N_NPP = cDates.size
        CollectionDates = np.zeros(N_NPP,dtype='datetime64[D]')
        for dd in range(0,N_NPP):
            day,month,year = cDates[dd].split("/")
            CollectionDates[dd]= np.datetime64(year+'-'+month+'-'+day)

        CollectionDates = np.sort(CollectionDates)
        PreviousDates = np.zeros((N_cores,N_NPP),dtype='datetime64[D]')
        AccumulationDays = np.zeros((N_cores,N_NPP))*np.nan
        FineRootsNPP = np.zeros((N_cores,N_NPP))*np.nan
        CoarseRootsNPP = np.zeros((N_cores,N_NPP))*np.nan

        for cc in range(0,N_cores):
            core_data = NPP_data[NPP_data['Core']==cc+1]
            N_collections = core_data.size

            for kk in range(0, N_collections):
                day,month,year = core_data['CollectionDate'][kk].split("/")
                dd = (CollectionDates == np.datetime64(year+'-'+month+'-'+day))
                
                if len(core_data['PreviousCollectionDate'][kk])==0:
                    PreviousDates[cc,dd]= np.datetime64(0,'D')
                else:
                    pday,pmonth,pyear = core_data['PreviousCollectionDate'][kk].split("/")            
                    PreviousDates[cc,dd]= np.datetime64(pyear+'-'+pmonth+'-'+pday)
                    AccumulationDays[cc,dd]=core_data['AccumulationDays'][kk]
                    CoarseRootsNPP[cc,dd]=core_data['CoarseRoots'][kk]
                    FineRootsNPP[cc,dd]=core_data['FineRoots'][kk]

        NPPDict['N_Cores']=N_cores
        NPPDict['CollectionDate']=CollectionDates
        NPPDict['PreviousCollectionDate']=PreviousDates
        NPPDict['AccumulationDays']=AccumulationDays
        NPPDict['CoarseRootNPP']=CoarseRootsNPP
        NPPDict['FineRootNPP']=FineRootsNPP
        NPPDict['Mask'] = np.isfinite(FineRootsNPP)

        RootStocksDict[plot_names[i]]=StocksDict
        RootNPPDict[plot_names[i]]=NPPDict

    return RootStocksDict,RootNPPDict



# Load spreadsheet of LAI derived from hemispherical photographs (courtesy of Terhi Riutta at Oxford).  LAI estimated using Hemisfer.
def load_field_LAI(LAI_file):
    datatype = {'names': ('ForestType','Plot', 'Subplot', 'LAI'), 'formats': ('S32','S32'
,'i8','f16')}
    hemisfer_LAI = np.genfromtxt(LAI_file, skiprows = 1, delimiter = ',',dtype=datatype)

    return hemisfer_LAI

def load_LAI_time_series(LAI_file):
    datatype = {'names': ('ForestType', 'Plot', 'Date', 'Subplot', 'SkyConditions', 'Exposure', 'LAI', 'Remarks', 'Qflag', 'Qreason'), 'formats': ('S16','S16','S10','i8','S16','i8','f8','S64','i8','S64')}
    LAI_data = np.genfromtxt(LAI_file, skiprows = 1, delimiter = ',',dtype=datatype)

    plot_names = np.unique(LAI_data['Plot'])
    N_plots = plot_names.size
    plot_dict = {}
    for pp in range(0,N_plots):
        LAI_dict = {}
        plot_data = LAI_data[LAI_data['Plot']==plot_names[pp]]
        N_subplots = np.unique(plot_data['Subplot']).size
        # set up some arrays for data output
        date_str = np.unique(plot_data['Date'])
        N_dates = date_str.size
        dates = np.zeros(N_dates,dtype = 'datetime64[D]')
        for dd in range(0,N_dates):
            day,month,year = date_str[dd].split("/")
            dates[dd] = np.datetime64(year+'-'+month+'-'+day)
        dates=np.sort(dates)

        # Need to loop through the dates column in the input array, and produce equivalent but in np.datetime64 so that cross referencing is easier
        date_ref = np.zeros(plot_data['Date'].size,dtype = 'datetime64[D]')
        for dd in range(0,plot_data['Date'].size):
            day,month,year = plot_data['Date'][dd].split("/")            
            date_ref[dd] =  np.datetime64(year+'-'+month+'-'+day)

        LAI = np.zeros((N_subplots,N_dates))
        for dd in range(0,N_dates):
            for ss in range(0,N_subplots):
                index=np.all((date_ref==dates[dd],plot_data['Subplot']==ss+1,plot_data['Exposure']==2,plot_data['Qflag']==1),axis=0)
                if np.sum(index)>0:
                    LAI[ss,dd] = plot_data['LAI'][index][0]
                else:
                    LAI[ss,dd] = np.nan

        LAI_dict['date']=dates.copy()
        LAI_dict['LAI']=LAI.copy()
        plot_dict[plot_names[pp]]=LAI_dict    
    return plot_dict
