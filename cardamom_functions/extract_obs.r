
###
## Function to extract obs needed for CARDAMOM
###

extract_obs<-function(latlon_wanted,lai_all,Csom_all,forest_all,Cwood_all,sand_clay_all,crop_man_all
		     ,burnt_all
		     ,ctessel_pft,site_name,start_year,end_year,timestep_days,spatial_type,resolution
		     ,grid_type,modelname) {
		      
    ###
    ## Get some LAI information
    ###
#    print("checking lai")
    if (lai_source == "MODIS") {
	# get lai
	lai=extract_modis_lai(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,lai_all,as.numeric(start_year):as.numeric(end_year))
	# assume default uncertainty (log scale)
	lai_unc=rep(1.5,length.out=length(lai))
    } else if (lai_source == "site_specific") {
	# read from .csv or netcdf
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	lai=read_site_specific_obs("LAI",infile)
	lai_unc=rep(1.5,length.out=length(lai))
    } else {
	lai=-9999
	lai_unc=-9999
    }

    ###
    ## Get some Csom information
    ###
#    print("checking Csom initial")
    if (Csom_source == "HWSD") {
	# could add other variables such as SOM (gC.m-2)
	som=extract_hwsd_Csom(spatial_type,resolution,grid_type,latlon_wanted,Csom_all)
    } else if (Csom_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
	som=read_site_specific_obs("Csom_initial",infile)
    } else {
	# assume no data available
	som=-9999
    }

    ###
    ## Get some sand / clay information
    ###

    if (sand_clay_source == "HWSD") {
	# could add other variables such as SOM (gC.m-2)
	sand_clay=extract_hwsd_sand_clay(spatial_type,resolution,grid_type,latlon_wanted,sand_clay_all)
	top_sand = sand_clay$top_sand ; bot_sand = sand_clay$bot_sand
	top_clay = sand_clay$bot_sand ; bot_clay = sand_clay$bot_clay
    } else if (sand_clay_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
	top_sand = read_site_specific_obs("top_sand_initial",infile)
	bot_sand = read_site_specific_obs("bot_sand_initial",infile)
	top_clay = read_site_specific_obs("top_clay_initial",infile)
	bot_clay = read_site_specific_obs("bot_clay_initial",infile)
    } else {
	# assume no data available
	top_sand = 40 ; bot_sand = 40
	top_clay = 15 ; bot_clay = 15
    }

    ###
    ## Get some crop management information
    ###

    if (ctessel_pft == 1 & crop_management_source == "sacks_crop_calendar") {
	# could add other variables such as SOM (gC.m-2)
	crop_dates=extract_sacks_crop_info(spatial_type,resolution,grid_type,latlon_wanted,crop_man_all)
	plant = crop_dates$plant ; plant_range = crop_dates$plant_range
	harvest = crop_dates$harvest ; harvest_range = crop_dates$harvest_range
    } else if (crop_management_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
	plant = read_site_specific_obs("plant_initial",infile)
	plant_range = read_site_specific_obs("plant_range_initial",infile)
	harvest = read_site_specific_obs("harvest_initial",infile)
	harvest_range = read_site_specific_obs("harvest_range_initial",infile)
    } else {
	# assume no data available
	plant = 304 ; plant_range = -9999 # 15
	harvest = 208 ; harvest_range = -9999 # 15
    }

    ###
    ## Get some Wood increment information (time series)
    ###
#    print("checking woodinc")
    if (woodinc_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	woodinc=read_site_specific_obs("woodinc",infile)
	woodinc_unc=rep(2,length.out=length(woodinc))
    } else {
	# assume no data available
	woodinc=-9999
	woodinc_unc=-9999
    }

    ###
    ## Get some GPP information (time series)
    ###
#    print("checking gpp")
    if (GPP_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")}
#	if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater.csv",sep="")}
	if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_copy.csv",sep="")}
	GPP = read_site_specific_obs("GPP",infile)
	GPP_unc = pmax(0.20,GPP * 0.20) #rep(1,length.out=length(GPP))
    } else {
	# assume no data available
	GPP = -9999
	GPP_unc = -9999
    }

    ###
    ## Get some Evapotranspiration information (time series)
    ###
#    print("checking gpp")
    if (Evap_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")}
#	if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater.csv",sep="")}
	if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_copy.csv",sep="")}
	Evap = read_site_specific_obs("Evap",infile)
	Evap_unc = pmax(0.2,Evap * 0.20) #rep(1,length.out=length(Evap))
	# borrow woody increment for soil evaporation in ACM_ET recalibration
	woodinc = read_site_specific_obs("soilevap",infile)
	woodinc_unc = pmax(0.2,woodinc * 0.20) #rep(1,length.out=length(woodinc))
    } else {
	# assume no data available
	Evap = -9999
	Evap_unc = -9999
    }

    ###
    ## Get some Reco information (time series)
    ###
#    print("checking Reco")
    if (Reco_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	Reco=read_site_specific_obs("Reco",infile)
	Reco_unc=rep(2,length.out=length(Reco))
    } else {
	# assume no data available
	Reco=-9999
	Reco_unc=-9999
    }

    ###
    ## Get some NEE information (time series)
    ###
#    print("checking nee")
    if (NEE_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	NEE=read_site_specific_obs("NEE",infile)
	NEE_unc=rep(2,length.out=length(NEE))
    } else {
	# assume no data available
	NEE=-9999
	NEE_unc=-9999
    }

    ###
    ## Get some Cfoliage information (initial conditions)
    ###
#    print("checking Cfoliage initial")
    if (Cfol_initial_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
	Cfol_initial=read_site_specific_obs("Cfol_initial",infile)
    } else {
	# assume no data available
	Cfol_initial=-9999
    }

    ###
    ## Get some Cwood information (initial conditions)
    ###
#    print("checking Cwood initial")
    if (Cwood_initial_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
	Cwood_initial=read_site_specific_obs("Cwood_initial",infile)
    } else if (Cwood_initial_source == "Avitabile") {
	# get Cwood
	output=extract_avitabile_biomass(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,Cwood_all)
	Cwood_initial=output$Cwood_stock
	Cwood_initial_unc=output$Cwood_stock_unc/output$Cwood_stock
    } else if (Cwood_initial_source == "mpi_biomass") {
	# get Cwood
	output=extract_mpi_biomass(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,Cwood_all)
	# locate time time data 
	tmp = which(output$Cwood_stock == max(output$Cwood_stock,na.rm=TRUE))[1]
	if (length(tmp) > 0) {
	    Cwood_initial=output$Cwood_stock[tmp]
	    Cwood_initial_unc=output$Cwood_stock_unc[tmp]/output$Cwood_stock[tmp]
	} else {
	    Cwood_initial=-9999
	    Cwood_initial_unc=-9999
	}
    } else {
	# assume no data available
	Cwood_initial=-9999 ; Cwood_initial_unc=-9999
    } 

    ###
    ## Get some Croots information (initial conditions)
    ###
#    print("checking Croots initial")
    if (Croots_initial_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
	Croots_initial=read_site_specific_obs("Croots_initial",infile)
    } else {
	# assume no data available
	Croots_initial=-9999
    }

    ###
    ## Get some Clitter information (initial conditions)
    ###
#    print("checking Clitter initial")
    if (Clit_initial_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
	Clit_initial=read_site_specific_obs("Clit_initial",infile)
    } else {
 	# assume no data available
	Clit_initial=-9999
    }
    ###
    ## Get some Cfoliage information (stock)
    ###
#    print("checking Cfoliage")
    if (Cfol_stock_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	Cfol_stock=read_site_specific_obs("Cfol_stock",infile)
	Cfol_stock_unc=read_site_specific_obs("Cfol_stock_unc",infile)
	if (length(which(Cfol_stock_unc != -9999)) > 0) {
	    # if we have uncertainty data convert to fraction to be compatable with analysis
	    tmp=which(Cfol_stock_unc == -9999)
	    Cfol_stock_unc=(Cfol_stock_unc/Cfol_stock)
	    Cfol_stock_unc[tmp]=-9999
	} else {
	    # on the other hand if not then we have no uncertainty info, so use default
	    Cfol_stock_unc=rep(0.38,length.out=length(Cfol_stock))
	}
    } else {
	# assume no data available
	Cfol_stock=-9999 ; Cfol_stock_unc=-9999
    }

    ###
    ## Get some Cwood information (stock)
    ###
#    print("checking Cwood")
    if (Cwood_stock_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	Cwood_stock=read_site_specific_obs("Cwood_stock",infile)
	Cwood_stock_unc=read_site_specific_obs("Cwood_stock_unc",infile)
	if (length(which(Cwood_stock_unc != -9999)) > 0) {
	    # if we have uncertainty data convert to fraction to be compatable with analysis
	    tmp=which(Cwood_stock_unc == -9999)
	    Cwood_stock_unc=(Cwood_stock_unc/Cwood_stock)
	    Cwood_stock_unc[tmp]=-9999
	} else {
	    # on the other hand if not then we have no uncertainty info, so use default
	    Cwood_stock_unc=rep(0.25,length.out=length(Cwood_stock))
	}
    } else if (Cwood_stock_source == "mpi_biomass") {
	# get Cwood
	output=extract_mpi_biomass(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,Cwood_all)
	Cwood_stock=output$Cwood_stock
	Cwood_stock_unc=output$Cwood_stock_unc/Cwood_stock
    } else {
	# assume no data available
	Cwood_stock=-9999
	Cwood_stock_unc=-9999
    } 

    ###
    ## Get some Cagb information (stock)
    ###
#    print("checking Cagb")
    if (Cagb_stock_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	Cagb_stock=read_site_specific_obs("Cagb_stock",infile)
	Cagb_stock_unc=read_site_specific_obs("Cagb_stock_unc",infile)
	if (length(which(Cagb_stock_unc != -9999)) > 0) {
	    # if we have uncertainty data convert to fraction to be compatable with analysis
	    tmp=which(Cagb_stock_unc == -9999)
	    Cagb_stock_unc=(Cagb_stock_unc/Cagb_stock)
	    Cagb_stock_unc[tmp]=-9999
	} else {
	    # on the other hand if not then we have no uncertainty info, so use default
	    Cagb_stock_unc=rep(0.25,length.out=length(Cagb_stock))
	}
    } else {
	# assume no data available
	Cagb_stock=-9999
	Cagb_stock_unc=-9999
    } 
    
    ###
    ## Get some Croots information (stock)
    ###
#    print("checking Croots")
    if (Croots_stock_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	Croots_stock=read_site_specific_obs("Croots_stock",infile)
	Croots_stock_unc=read_site_specific_obs("Croots_stock_unc",infile)
	if (length(which(Croots_stock_unc != -9999)) > 0) {
	    # if we have uncertainty data convert to fraction to be compatable with analysis
	    tmp=which(Croots_stock_unc == -9999)
	    Croots_stock_unc=(Croots_stock_unc/Croots_stock)
	    Croots_stock_unc[tmp]=-9999
	} else {
	    # on the other hand if not then we have no uncertainty info, so use default
	    Croots_stock_unc=rep(0.44,length.out=length(Croots_stock))
	}
    } else {
	# assume no data available
	Croots_stock=-9999
	Croots_stock_unc=-9999
    }

    ###
    ## Get some Clitter information (stock)
    ###
#    print("checking Clitter")
    if (Clit_stock_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	Clit_stock=read_site_specific_obs("Clit_stock",infile)
	Clit_stock_unc=read_site_specific_obs("Clit_stock_unc",infile)
	if (length(which(Clit_stock_unc != -9999)) > 0) {
	    # if we have uncertainty data convert to fraction to be compatable with analysis
	    tmp=which(Clit_stock_unc == -9999)
	    Clit_stock_unc=(Clit_stock_unc/Clit_stock)
	    Clit_stock_unc[tmp]=-9999
	} else {
	    # on the other hand if not then we have no uncertainty info, so use default
	    Clit_stock_unc=rep(0.38,length.out=length(Clit_stock))
	}
    } else {
 	# assume no data available
	Clit_stock=-9999 ; Clit_stock_unc=-9999
    }

    ###
    ## Get some Csom information (stock)
    ###
#    print("checking Csom")
    if (Csom_stock_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	Csom_stock=read_site_specific_obs("Csom_stock",infile)
	Csom_stock_unc=read_site_specific_obs("Csom_stock_unc",infile)
	if (length(which(Csom_stock_unc != -9999)) > 0) {
	    # if we have uncertainty data convert to fraction to be compatable with analysis
	    tmp=which(Csom_stock_unc == -9999)
	    Csom_stock_unc=(Csom_stock_unc/Csom_stock)
	    Csom_stock_unc[tmp]=-9999
	} else {
	    # on the other hand if not then we have no uncertainty info, so use default
	    Csom_stock_unc=rep(0.24,length.out=length(Csom_stock))
	}
    } else {
 	# assume no data available
	Csom_stock=-9999 ; Csom_stock_unc=-9999
    }

    ###
    ## Get some Cstem information (stock)
    ###
#    print("checking Cstem")
    if (Cstem_stock_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	Cstem_stock=read_site_specific_obs("Cstem_stock",infile)
	Cstem_stock_unc=read_site_specific_obs("Cstem_stock_unc",infile)
	if (length(which(Cstem_stock_unc != -9999)) > 0) {
	    # if we have uncertainty data convert to fraction to be compatable with analysis
	    tmp=which(Cstem_stock_unc == -9999)
	    Cstem_stock_unc=(Cstem_stock_unc/Cstem_stock)
	    Cstem_stock_unc[tmp]=-9999
	} else {
	    # on the other hand if not then we have no uncertainty info, so use default
	    Cstem_stock_unc=rep(0.24,length.out=length(Cstem_stock))
	}
    } else {
 	# assume no data available
	Cstem_stock=-9999 ; Cstem_stock_unc=-9999
    }

    ###
    ## Get some Cbranch information (stock)
    ###
#    print("checking Cbranch")
    if (Cbranch_stock_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	Cbranch_stock=read_site_specific_obs("Cbranch_stock",infile)
	Cbranch_stock_unc=read_site_specific_obs("Cbranch_stock_unc",infile)
	if (length(which(Cbranch_stock_unc != -9999)) > 0) {
	    # if we have uncertainty data convert to fraction to be compatable with analysis
	    tmp=which(Cbranch_stock_unc == -9999)
	    Cbranch_stock_unc=(Cbranch_stock_unc/Cbranch_stock)
	    Cbranch_stock_unc[tmp]=-9999
	} else {
	    # on the other hand if not then we have no uncertainty info, so use default
	    Cbranch_stock_unc=rep(0.24,length.out=length(Cbranch_stock))
	}
    } else {
 	# assume no data available
	Cbranch_stock=-9999 ; Cbranch_stock_unc=-9999
    }

    ###
    ## Get some Ccoarseroot information (stock)
    ###
#    print("checking Ccoarseroot")
    if (Ccoarseroot_stock_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	Ccoarseroot_stock=read_site_specific_obs("Ccoarseroot_stock",infile)
	Ccoarseroot_stock_unc=read_site_specific_obs("Ccoarseroot_stock_unc",infile)
	if (length(which(Ccoarseroot_stock_unc != -9999)) > 0) {
	    # if we have uncertainty data convert to fraction to be compatable with analysis
	    tmp=which(Ccoarseroot_stock_unc == -9999)
	    Ccoarseroot_stock_unc=(Ccoarseroot_stock_unc/Ccoarseroot_stock)
	    Ccoarseroot_stock_unc[tmp]=-9999
	} else {
	    # on the other hand if not then we have no uncertainty info, so use default
	    Ccoarseroot_stock_unc=rep(0.24,length.out=length(Ccoarseroot_stock))
	}
    } else {
 	# assume no data available
	Ccoarseroot_stock=-9999 ; Ccoarseroot_stock_unc=-9999
    }

    ###
    ## Get some Cfolmax information (stock)
    ###
#    print("checking Cfolmax")
    if (Cfolmax_stock_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	Cfolmax_stock=read_site_specific_obs("Cfolmax_stock",infile)
	Cfolmax_stock_unc=read_site_specific_obs("Cfolmax_stock_unc",infile)
	if (length(which(Cfolmax_stock_unc != -9999)) > 0) {
	    # if we have uncertainty data convert to fraction to be compatable with analysis
	    tmp=which(Cfolmax_stock_unc == -9999)
	    Cfolmax_stock_unc=(Cfolmax_stock_unc/Cfolmax_stock)
	    Cfolmax_stock_unc[tmp]=-9999
	} else {
	    # on the other hand if not then we have no uncertainty info, so use default
	    Cfolmax_stock_unc=rep(0.24,length.out=length(Cfolmax_stock))
	}
    } else {
 	# assume no data available
	Cfolmax_stock=-9999 ; Cfolmax_stock_unc=-9999
    }

    ###
    ## Get some deforestation information (fraction time series)
    ###
#    print("checking deforestation")
    if (deforestation_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	deforestation=read_site_specific_obs("deforestation",infile)
	forest_management=read_site_specific_obs("management",infile)
	yield_class=0 #read_site_specific_obs("yield_class",infile)
	age=read_site_specific_obs("age",infile)
	age=age[1] # we only want the age at the beginning of the simulation
    } else if (deforestation_source == "combined_dataset" | deforestation_source == "GFW") {
	output=extract_forestry_information(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,forest_all,start_year,end_year,ctessel_pft)
	ctessel_pft=output$ctessel_pft
	deforestation=output$deforestation
	yield_class=output$yield_class
	age=output$age #; if (age > 1) {age = -9999}
	forest_management=2
    } else {
 	# assume no data available
	deforestation=0
	forest_management=1
	yield_class=0
	age=-9999
#	if (modelname == "DALEC_GSI_FR" | modelname == "DALECN_GSI_FR") {age=30}
    }

    ###
    ## Get some burnt area information (fraction time series)
    ###
#    print("checking burnt_area")
    if (burnt_area_source == "site_specific") {
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
	burnt_area=read_site_specific_obs("burnt_area",infile)
    } else if (burnt_area_source == " "){
 	# assume no data available
	burnt_area=0
    } else {
	burnt_area=extract_burnt_area_information(latlon_wanted,timestep_days,spatial_type,grid_type,resolution,start_year,end_year,burnt_all)
    }

    # return output now
    return(list(LAT=latlon_wanted[1],LAI=lai,SOM=som,GPP=GPP,Evap=Evap,NEE=NEE,Reco=Reco,woodinc=woodinc
		,Cfol_stock=Cfol_stock,Cwood_stock=Cwood_stock,Cagb_stock=Cagb_stock,Croots_stock=Croots_stock,Clit_stock=Clit_stock
		,Cfol_initial=Cfol_initial,Cwood_initial=Cwood_initial,Cwood_initial_unc=Cwood_initial_unc,Croots_initial=Croots_initial,Clit_initial=Clit_initial
		,Csom_stock=Csom_stock,deforestation=deforestation,burnt_area=burnt_area,ctessel_pft=ctessel_pft,yield_class=yield_class,age=age
		,forest_management=forest_management,LAI_unc=lai_unc,GPP_unc=GPP_unc,Evap_unc=Evap_unc,NEE_unc=NEE_unc,Reco_unc=Reco_unc,woodinc_unc=woodinc_unc
		,Cfol_stock_unc=Cfol_stock_unc,Cwood_stock_unc=Cwood_stock_unc,Cagb_stock_unc=Cagb_stock_unc,Croots_stock_unc=Croots_stock_unc
		,Clit_stock_unc=Clit_stock_unc,Csom_stock_unc=Csom_stock_unc,Cstem_stock=Cstem_stock,Cstem_stock_unc=Cstem_stock_unc
		,Cbranch_stock=Cbranch_stock,Cbranch_stock_unc=Cbranch_stock_unc,Ccoarseroot_stock=Ccoarseroot_stock
		,Ccoarseroot_stock_unc=Ccoarseroot_stock_unc,Cfolmax_stock=Cfolmax_stock,Cfolmax_stock_unc=Cfolmax_stock_unc
		,top_sand=top_sand,bot_sand=bot_sand,top_clay=top_clay,bot_clay=bot_clay,plant=plant,plant_range=plant_range
		,harvest=harvest,harvest_range=harvest_range))

}

