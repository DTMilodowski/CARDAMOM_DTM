
###
## Function to load met data from global field ECMWF data
## subsequently extracted in extract_met_drivers.txt
###

load_met_fields_for_extraction<-function(latlon_in,met_source,modelname,startyear,endyear) {
    # let the user know this might take some time
    print("Loading global met fields for subsequent sub-setting ...")
    # declare timing variables
    t_grid=0

    if (met_source == "site_specific") {

	# contruct output
	met_all=list(site_specific=TRUE)
  
    } else {

	if (met_source == "ECMWF") {
	    # declare variable ids needed to select files / infile variables
	    varid=c("SWdown","Tair","Rainf","Qair","PSurf","Wind_E","Wind_N") ; infile_varid=c("SSRD","T","LSP","Q","LNSP","U","V")
	    # open first ecmwf file to extract needed information
	    input_file_1=paste(path_to_met_source,varid[1],"_ei",startyear,"01.nc",sep="") 
	    data1=nc_open(input_file_1)

	    # get timing variable
	    hour_since1=ncvar_get(data1,"time")
	    # step size in hours of the ecmwf files
	    input_step_size=hour_since1[2]-hour_since1[1] 
	    steps_in_day=24/input_step_size
		    
	    # extract location variables
	    lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "lon")
	    # expand the one directional values here into 2 directional
	    lat_dim=length(lat) ; long_dim=length(long)
	    long=array(long,dim=c(long_dim,lat_dim))
	    lat=array(lat,dim=c(lat_dim,long_dim)) ; lat=t(lat)

	    # If we are using the GSI model we need the 21 days (or month) before the start date of the simulation, so we need to check if we have this information
	    # NOTE that this section of code is duplicated for each of the available datasets because of differences in storage and file name
	    extra_year=FALSE
    #	if (grepl("GSI",modelname)) {
		present=0
		for (lag in seq(1,length(varid))) {
		    input_file_1=paste(path_to_met_source,varid[1],"_ei",as.character(as.numeric(startyear)-1),"01.nc",sep="") 
		    # check whether the files exist or not
		    if (file.exists(input_file_1)) {present=present+1}
		}
		# if all the files exist then we will use them
		if (present == length(varid)) {extra_year=TRUE}
    #	}

	} else if (met_source == "ERA") {    
	    # declare variable ids needed to select files / infile variables
	    varid=c("sw_radiation_daily_mean","airt_daily_max","precipitation_daily_mean","vpd_daily_mean","sfc_pressure_daily_mean","airt_daily_min","wind_spd_daily_mean") 
	    infile_varid=c("daily_swrad","airt_max","daily_precip","vpd_mean","sfc_pressure_mean","airt_min","wind_spd")

	    # open first ecmwf file to extract needed information
	    input_file_1=paste(path_to_met_source,varid[1],"_",startyear,"01.nc",sep="") 
	    data1=nc_open(input_file_1)

	    # get timing variable
	    input_step_size=24
	    steps_in_day=1
		    
	    # extract location variables
	    lat=ncvar_get(data1, "Latitude") ; long=ncvar_get(data1, "Longitude")
	    # expand the one directional values here into 2 directional
	    lat_dim=length(lat) ; long_dim=length(long)
	    long=array(long,dim=c(long_dim,lat_dim))
	    lat=array(lat,dim=c(lat_dim,long_dim)) ; lat=t(lat)

	    # If we are using the GSI model we need the 21 days (or month) before the start date of the simulation, so we need to check if we have this information
	    # NOTE that this section of code is duplicated for each of the available datasets because of differences in storage and file name
	    extra_year=FALSE
    #	if (grepl("GSI",modelname)) {
		present=0
		for (lag in seq(1,length(varid))) {
		    input_file_1=paste(path_to_met_source,varid[1],"_",as.character(as.numeric(startyear)-1),"01.nc",sep="") 
		    # check whether the files exist or not
		    if (file.exists(input_file_1)) {present=present+1}
		}
		# if all the files exist then we will use them
		if (present == length(varid)) {extra_year=TRUE}
    #	}

	} else if (met_source == "PRINCETON") {
	    # declare variable ids needed to select files / infile variables
	    varid=c("dswrf","tas","prcp","shum","pres","wind") ; infile_varid=c("dswrf","tas","prcp","shum","pres","wind")
	    # open first ecmwf file to extract needed information
	    input_file_1=paste(path_to_met_source,varid[1],"_3hourly_",startyear,"-",startyear,".nc",sep="") 
	    data1=nc_open(input_file_1)

	    # get timing variable
	    hour_since1=ncvar_get(data1,"time")
	    # step size in hours of the ecmwf files
	    input_step_size=hour_since1[2]-hour_since1[1] 
	    steps_in_day=24/input_step_size
		    
	    # extract location variables
	    lat=ncvar_get(data1, "latitude") ; long=ncvar_get(data1, "longitude")
	    # expand the one directional values here into 2 directional
	    lat_dim=length(lat) ; long_dim=length(long)
	    long=array(long,dim=c(long_dim,lat_dim))
	    lat=array(lat,dim=c(lat_dim,long_dim)) ; lat=t(lat)

	    # If we are using the GSI model we need the 21 days (or month) before the start date of the simulation, so we need to check if we have this information
	    # NOTE that this section of code is duplicated for each of the available datasets because of differences in storage and file name
	    extra_year=FALSE
    #	if (grepl("GSI",modelname)) {
		present=0
		for (lag in seq(1,length(varid))) {
		    input_file_1=paste(path_to_met_source,varid[lag],"_3hourly_",as.character(as.numeric(startyear)-1),"-",as.character(as.numeric(startyear)-1),".nc",sep="") 
		    # check whether the files exist or not
		    if (file.exists(input_file_1)) {present=present+1}
		}
		# if all the files exist then we will use them
		if (present == length(varid)) {extra_year=TRUE}	
    #	}

	} else if (met_source == "CHESS") {
	    # declare variable ids needed to select files / infile variables
	    varid=c("rsds","tas","precip","huss","psurf","dtr","sfcWind") ; infile_varid=c("rsds","tas","precip","huss","psurf","dtr","sfcWind")
	    # open first ecmwf file to extract needed information
	    input_file_1=paste(path_to_met_source,"chess_",varid[1],"_",startyear,"01.nc",sep="") 
	    data1=nc_open(input_file_1)

	    # get timing variable
	    input_step_size=24
	    steps_in_day=1
		    
	    # extract location variables
	    lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "lon")

	    # If we are using the GSI model we need the 21 days (or month) before the start date of the simulation, so we need to check if we have this information
	    # NOTE that this section of code is duplicated for each of the available datasets because of differences in storage and file name
	    extra_year=FALSE
    #	if (grepl("GSI",modelname)) {
		present=0
		for (lag in seq(1,length(varid))) {
		    input_file_1=paste(path_to_met_source,"chess_",varid[lag],"_",as.character(as.numeric(startyear)-1),"12.nc",sep="") 
		    # check whether the files exist or not
		    if (file.exists(input_file_1)) {present=present+1}
		}
		# if all the files exist then we will use them
		if (present == length(varid)) {extra_year=TRUE}	
    #	}
	}

	# now find out which pixels we will filter through
	# if lat / long is 2 dimensional rather than vector we need to adjust for this
	lat_dim=dim(lat)[2] ; long_dim=dim(long)[1]
	# convert input data long to conform to what we need
	check1=which(long > 180) ; if (length(check1) > 0) { long[check1]=long[check1]-360 }
	# which locations are within the desired zone
	remove_lat=intersect(which(lat < (max(latlon_in[,1])+1.0)),which(lat > (min(latlon_in[,1])-1.0)))
	remove_long=intersect(which(long < (max(latlon_in[,2])+1.0)),which(long > (min(latlon_in[,2])-1.0)))
	if (length(check1) > 0) { long[check1]=long[check1]+360 }
	# now find common where out in both contexts
	remove_lat=intersect(remove_lat,remove_long)
	# update both variables because of common matrix
	remove_long=remove_lat
	# adjust for matrix rather than vector arrangement
	remove_lat=remove_lat/dim(lat)[1]
	remove_long=(remove_long-(floor(remove_lat)*dim(lat)[1]))+1
	remove_lat=ceiling(remove_lat)
	# update new dimensions
	lat_dim=length(min(remove_lat):max(remove_lat)) ; long_dim=length(min(remove_long):max(remove_long))
	lat=lat[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)] ; long=long[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
	# next deal with extracting the wheat from the chaff
	# this should always be the swrad variable as we can easily put a lower limit
	tmp1=ncvar_get(data1,infile_varid[1])
        tmp1[which(is.na(tmp1) == TRUE)] = -9999
	# now because Burn during its update period has conflicting R packages and this means that it does not always 
	# correctly identify the missing data. So we use the Short Wave Radiation as a condition for searching for this as we know that
	# swrad cannot be less than 0
	tmp1[which(tmp1 < 0)] = NA
	# this section is key as near land sea borders it is possible for the nearest lat/long location to actually be a sea pixel
	wheat_from_chaff=which(is.na(tmp1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),1]) == FALSE)

	# convert input data long to conform to what we need
	check1=which(long > 180) ; if (length(check1) > 0) { long[check1]=long[check1]-360 }
	# now filter through the reduced dataset for the specific locations
	# not this selection by met_in$wheat vectorises the lat and long variables therefore we need to use closest2 option 1 (I think...)
	output=lapply(1:dim(latlon_in)[1],FUN=closest2d,lat=lat[wheat_from_chaff],long=long[wheat_from_chaff],lat_in=latlon_in[,1],long_in=latlon_in[,2],nos_dim=1) 
	var1_out=unlist(output) ; rm(output)
	# select the correct location values for the original vector form the wheat_from_chaff
	wheat_from_chaff=wheat_from_chaff[var1_out]
	# return long to 0-360
	if (length(check1) > 0) { long[check1]=long[check1]+360 }

	# tidy
	nc_close(data1)

	# define timing variables
	years_to_load=as.numeric(startyear):as.numeric(endyear)
	if (extra_year) {
	    load_years=c(as.character(as.numeric(startyear)-1),as.numeric(startyear):as.numeric(endyear))
	} else {
	    load_years=as.numeric(startyear):as.numeric(endyear)
	}
	print("have finished determining met extraction parameters, now for the main load")
	# extract the data into lists which we will contruct.
	if (use_parallel) {
	    cl <- makeCluster(min(length(load_years),numWorkers), type = "PSOCK")    
	    # load R libraries in cluster
	    clusterExport(cl,"load_r_libraries") ; clusterEvalQ(cl, load_r_libraries())
	    if (met_source == "CHESS") {
		print("...bonus CHESS diurnal temperature range and wind spd to be read in") 
		var6_out_list=parLapply(cl,load_years,fun=load_met_function,varid=varid[6],infile_varid=infile_varid[6],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff) 
		wind_out_list=parLapply(cl,load_years,fun=load_met_function,varid=varid[7],infile_varid=infile_varid[7],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    } else if (met_source == "ERA") {
		print("...bonus ERA-Interim minimum temperature and wind spd to be read in") 
		var6_out_list=parLapply(cl,load_years,fun=load_met_function,varid=varid[6],infile_varid=infile_varid[6],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff) 
		wind_out_list=parLapply(cl,load_years,fun=load_met_function,varid=varid[7],infile_varid=infile_varid[7],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    } else if (met_source == "ECMWF") {
		print("...bonus ECMWF wind speed north and east to be read")
		wind_out_list=parLapply(cl,load_years,fun=load_met_function,varid=varid[6],infile_varid=infile_varid[6],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)            
		windN_out_list=parLapply(cl,load_years,fun=load_met_function,varid=varid[7],infile_varid=infile_varid[7],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    } else if (met_source == "PRINCETON") {
		print("...bonus wind spd to be read")
		wind_out_list=parLapply(cl,load_years,fun=load_met_function,varid=varid[6],infile_varid=infile_varid[6],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    }
	    var1_out_list=parLapply(cl,load_years,fun=load_met_function,varid=varid[1],infile_varid=infile_varid[1],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff) 
	    print("...met load 20 %")
	    var2_out_list=parLapply(cl,load_years,fun=load_met_function,varid=varid[2],infile_varid=infile_varid[2],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    print("...met load 40 %")
	    var3_out_list=parLapply(cl,load_years,fun=load_met_function,varid=varid[3],infile_varid=infile_varid[3],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)  
	    print("...met load 60 %")
	    var4_out_list=parLapply(cl,load_years,fun=load_met_function,varid=varid[4],infile_varid=infile_varid[4],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    print("...met load 80 %")
	    var5_out_list=parLapply(cl,load_years,fun=load_met_function,varid=varid[5],infile_varid=infile_varid[5],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)    
	    print("...met load 100 %")
	    stopCluster(cl)	
	} else {
	    # or use serial
	    if (met_source == "CHESS") {
		print("...bonus CHESS daily temperature range and wind speed file to be read in") 
		var6_out_list=lapply(load_years,FUN=load_met_function,varid=varid[6],infile_varid=infile_varid[6],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
		wind_out_list=lapply(load_years,FUN=load_met_function,varid=varid[7],infile_varid=infile_varid[7],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaf)
	    } else if (met_source == "ERA") {
		print("...bonus ERA-Interim minimum temperature and wind speed to be to be read in") 
		var6_out_list=lapply(load_years,FUN=load_met_function,varid=varid[6],infile_varid=infile_varid[6],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
		wind_out_list=lapply(load_years,FUN=load_met_function,varid=varid[7],infile_varid=infile_varid[7],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    } else if (met_source == "ECMWF") {
		print("...bonus wind speed to be read in")
		wind_out_list=lapply(load_years,FUN=load_met_function,varid=varid[6],infile_varid=infile_varid[6],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
		windN_out_list=lapply(load_years,FUN=load_met_function,varid=varid[7],infile_varid=infile_varid[7],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    } else if (met_source == "PRINCETON"){
		print("...bonus wind speed to be read in")
		wind_out_list=lapply(load_years,FUN=load_met_function,varid=varid[6],infile_varid=infile_varid[6],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    }
	    var1_out_list=lapply(load_years,FUN=load_met_function,varid=varid[1],infile_varid=infile_varid[1],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    print("...met load 20 %")
	    var2_out_list=lapply(load_years,FUN=load_met_function,varid=varid[2],infile_varid=infile_varid[2],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    print("...met load 40 %")
	    var3_out_list=lapply(load_years,FUN=load_met_function,varid=varid[3],infile_varid=infile_varid[3],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    print("...met load 60 %")
	    var4_out_list=lapply(load_years,FUN=load_met_function,varid=varid[4],infile_varid=infile_varid[4],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    print("...met load 80 %")
	    var5_out_list=lapply(load_years,FUN=load_met_function,varid=varid[5],infile_varid=infile_varid[5],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
	    print("...met load 100 %")
	} # parallel option

    # hack not to be left alone
    if (met_source == "PRINCETON") {
    #    copy_to=readline("employ hack fix for duke meteorology (y/n)")
    #    if (copy_to == "y" ) {
	    # 1983 has an error in the short wave radiation so we will loop over it
	    this_one = which(load_years == "1983")
	    if (length(this_one) > 0) {var1_out_list[[this_one]]$var_out=var1_out_list[[(this_one+1)]]$var_out}
    #    }
    }

	# user update
	print("...beginning restructuring of meteorological datasets")
	var1_out=0 ; var2_out=0 ; var3_out=0 ; var4_out=0 ; var5_out=0 ; var6_out=0 ; wind_out=0 ; tmp_out=0; t_grid=0
	for (i in seq(1, length(var1_out_list))) {var1_out=append(var1_out,var1_out_list[[i]]$var_out) ; t_grid=append(t_grid,var1_out_list[[i]]$t_grid)}
	rm(var1_out_list) 
	for (i in seq(1, length(var2_out_list))) {var2_out=append(var2_out,var2_out_list[[i]]$var_out)}
	rm(var2_out_list) 
	for (i in seq(1, length(var3_out_list))) {var3_out=append(var3_out,var3_out_list[[i]]$var_out)}
	rm(var3_out_list) 
	for (i in seq(1, length(var4_out_list))) {var4_out=append(var4_out,var4_out_list[[i]]$var_out)}
	rm(var4_out_list) 
	for (i in seq(1, length(var5_out_list))) {var5_out=append(var5_out,var5_out_list[[i]]$var_out)}
	rm(var5_out_list) 
	if (met_source == "CHESS" | met_source == "ERA") {
	    for (i in seq(1, length(var6_out_list))) {var6_out=append(var6_out,var6_out_list[[i]]$var_out)}
	    rm(var6_out_list)
	    for (i in seq(1, length(wind_out_list))) {wind_out=append(wind_out,wind_out_list[[i]]$var_out)}
	    rm(wind_out_list)
	} else if (met_source == "PRINCETON" | met_source == "CHESS"){
	    for (i in seq(1, length(wind_out_list))) {wind_out=append(wind_out,wind_out_list[[i]]$var_out)}
	    rm(wind_out_list)
	} else if (met_source == "ECMWF") {
	    for (i in seq(1, length(wind_out_list))) {wind_out=append(wind_out,wind_out_list[[i]]$var_out)}
	    rm(wind_out_list)
	    for (i in seq(1, length(windN_out_list))) {tmp_out=append(tmp_out,windN_out_list[[i]]$var_out)}
	    rm(windN_out_list)
	    # convert wind components into single magnitude 
	    wind_out = sqrt((wind_out**2) + (tmp_out**2))
	    # and remove excess 
	    rm(tmp_out)
	}

	# remove initial value
	var1_out=var1_out[-1] ; var2_out=var2_out[-1]
	var3_out=var3_out[-1] ; var4_out=var4_out[-1]
	var5_out=var5_out[-1] ; var6_out=var6_out[-1]
	wind_out=wind_out[-1] ; gc()

	if (met_source == "ECMWF") {
	    # if ECMWF then convert log(Pa) -> Pa
	    var5_out=exp(var5_out)
	}

	# we want total t_grid only
	t_grid=sum(t_grid)

	# generate some additional timing and data
	# NOTE that if we are using an extra of data for the purposes of a GSI model run then the number of days should be the number intended for simulation 
	# not to match with the met files just at the moment.
	# This will be corrected later on when the location specific data are extracted
	print("...generating day of year variables")
	for (yr in seq(1,length(years_to_load))) {
	    # if this includes the extra year add it on to the beginning
	    if (extra_year) {extra_nos_days=nos_days_in_year(load_years[1])}
	    # is current year a leap or not
	    nos_days=nos_days_in_year(years_to_load[yr])

	    # reconstruct arrays before doing site level
	    if (yr == 1) {
		doy=seq(1,nos_days)
	  #	    co2=rep(c(450,445,440,440,440,430,420,400,390,380,360,350,340,330,330,340,355,370,390,400,420,440,445,450),length.out=(nos_days*steps_in_day), each=floor((nos_days*steps_in_day)/24))
		# mass of dry air = 28.97(g/mol) ; mass of co2 = 44.01 (g/mol); *1e-6 scale from umol -> mol
		co2=rep(380,length.out=(nos_days*steps_in_day)) # default ppm 570 / 370 face
		if (extra_year) {co2=append(rep(380,length.out=(extra_nos_days*steps_in_day)),co2)}
	    } else {
		doy=append(doy,seq(1,nos_days))
	  #	    co2=append(co2,rep(c(450,445,440,440,440,430,420,400,390,380,360,350,340,330,330,340,355,370,390,400,420,440,445,450),length.out=(nos_days*steps_in_day), each=floor((nos_days*steps_in_day)/24))+(2*yr))
		co2=append(co2,rep(380,length.out=(nos_days*steps_in_day)))
	    }
	} # end of years loop

	# create day of run variable
	run_day=seq(1:length(doy))

	# user update
	print("...preparing final met output variables")
	# restructure to output variables; currently wheat from chaff only used for CHESS but actually might be better to use for all cases later on when recoding occurs...
	var1_out=array(var1_out, dim=c(length(wheat_from_chaff),t_grid))
	var2_out=array(var2_out, dim=c(length(wheat_from_chaff),t_grid))
	var3_out=array(var3_out, dim=c(length(wheat_from_chaff),t_grid))
	var4_out=array(var4_out, dim=c(length(wheat_from_chaff),t_grid))
	var5_out=array(var5_out, dim=c(length(wheat_from_chaff),t_grid))
	var6_out=array(var6_out, dim=c(length(wheat_from_chaff),t_grid))
	wind_out=array(wind_out, dim=c(length(wheat_from_chaff),t_grid))

	# output variables
	met_all=list(input_step_size=input_step_size,steps_in_day=steps_in_day,lat=lat,long=long,run_day=run_day,wheat=wheat_from_chaff,t_grid=t_grid,airt=var2_out,swrad=var1_out,co2=co2,doy=doy,precip=var3_out,sp_humidity=var4_out,pressure=var5_out,extra_year=extra_year,diurnal_trange=var6_out,wind_spd=wind_out)

	# quick sanity check
	if(min(as.vector(met_all$swrad) < -5)) {print(paste("SW_RAD summary: ",summary(as.vector(met_all$swrad)),sep="")) } #; stop('SW Radiation has values less than zero - oh dear ')}
	met_all$swrad[which(met_all$swrad < 0)] = 0

	# clean up loose memory
	rm(lat_dim,long_dim,remove_lat,remove_long,years_to_load,load_years,tmp1,check1,nos_days)
	rm(input_step_size,steps_in_day,lat,long,run_day,wheat_from_chaff,t_grid,var2_out,var1_out,co2,doy,var3_out,var4_out,var5_out,extra_year,var6_out,wind_out)

    }
    gc(reset=TRUE,verbose=FALSE)
    return(met_all)

} # function end
## Use byte compile
load_met_fields_for_extraction<-cmpfun(load_met_fields_for_extraction)
