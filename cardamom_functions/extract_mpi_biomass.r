

extract_mpi_biomass<- function(timestep_days,spatial_type,resolution,grid_type,latlon_in,Cwood_all) {

	# find the nearest location
	output=closest2d(1,Cwood_all$lat,Cwood_all$long,latlon_in[1],latlon_in[2],3)
	#i1=unlist(output)[2] ; j1=unlist(output)[1]
	i1=unlist(output)[1] ; j1=unlist(output)[2]
	print(paste("MPI Biomass data extracted for current location ",Sys.time(),sep=""))

#	# work out number of pixels to average over
#	if (spatial_type == "grid") {
#	    if (grid_type == "wgs84") {
#		# calculate pixel area and convert from m2 -> km2
#		area=calc_pixel_area(latlon_in[1],latlon_in[2],resolution)*1e-6
#		radius=max(0,ceiling(sqrt(area)*0.5))
#	    } else if (grid_type == "UK") {
#		radius=max(0,floor(1*resolution*1e-3*0.5))
#	    } else {
#		stop("have not specified the grid used in this analysis")
#	    }
#	} else {
#	    radius=0
#	}

	# work out number of pixels to average over
	if (spatial_type == "grid") {
	    if (grid_type == "wgs84") {
		# resolution of the product
		product_res = (Cwood_all$lat[1,2]-Cwood_all$lat[1,1])+(Cwood_all$long[2,1]-Cwood_all$long[1,1])
		product_res = product_res * 0.5
		# radius is ceiling of the ratio of the product vs analysis ratio
		radius = ceiling(resolution / product_res)
	    } else if (grid_type == "UK") {
		radius=max(0,floor(1*resolution*1e-3*0.5))
	    } else {
		stop("have not specified the grid used in this analysis")
	    }
	} else {
	    radius=0
	}

	# work out average areas
	average_i=(i1-radius):(i1+radius) ; average_j=(j1-radius):(j1+radius)
	# carry out averaging
	Cwood=array(NA, dim=c(dim(Cwood_all$biomass_all)[3]))
	Cwood_unc=array(NA, dim=c(dim(Cwood_all$biomass_all_unc)[3]))

	for (n in seq(1, dim(Cwood_all$biomass_all)[3])) {
	    Cwood[n]=mean(Cwood_all$biomass_all[average_i,average_j,n], na.rm=TRUE)
	    Cwood_unc[n]=mean(Cwood_all$biomass_all_unc[average_i,average_j,n], na.rm=TRUE)
	} #  for loop

	# convert missing data back to -9999
	Cwood[which(is.na(Cwood))]=-9999.0 ; Cwood_unc[which(is.na(Cwood_unc))]=-9999.0
	# next work out how many days we should have in the year
	doy_out=0
	for (i in seq(1, length(years_to_do))) {
	    # is current year a leap or not
	    nos_days = 365
	    mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/4))*4
	    if (mod == 0) {
		nos_days = 366
		mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/100))*100
		if (mod == 0) {
		    nos_days  = 365
		    mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/400))*400
		    if (mod == 0) {
			nos_days  = 366
		    }
		}
	    }
	    # count up days needed
	    doy_out=append(doy_out,1:nos_days)
	}
	doy_out=doy_out[-1]

	# just incase there is no missing data we best make sure there is a value which can be assessed
	if (length(Cwood_all$missing_years) == 0) { Cwood_all$missing_years=1066 }

	# declare output variable
	Cwood_out=array(-9999, dim=length(doy_out))
	Cwood_unc_out=array(-9999, dim=length(doy_out))
	# now line up the obs days with all days
	b=1 ; i=1 ; a=1 ; start_year=as.numeric(years_to_do[1]) 
	while (b <= length(Cwood_all$doy_obs)) {
	    # if we are in a year which is missing then we do not allow consideration of DOY
	    if (start_year != Cwood_all$missing_years[a]) {
		if (doy_out[i] == Cwood_all$doy_obs[b]) {
		    Cwood_out[i]=Cwood[b] ; Cwood_unc_out[i]=Cwood_unc[b] ; b=b+1
		} # end if doy matches
	    } # end if missing year
	    # but we do keep counting through the total vector length which we expect
	    i=i+1
	    # each time we come back to doy_out[i]==1 we need to count on the year
	    if (doy_out[i] == 1) {
		# and if we have just been in a missing year we need to count on the missing years vector to
		if (start_year == Cwood_all$missing_years[a]) {a=min(length(Cwood_all$missing_years),a+1)}
		start_year=start_year+1 
	    } # end if doy_out[i] == 1
	} # end while condition

	if (length(timestep_days) == 1 & timestep_days[1] == 1) {
	    # well actually we do nothing
	} else {
	    # generally this now deals with time steps which are not daily.
	    # However if not monthly special case
	    if (length(timestep_days) == 1) {
		run_day_selector=seq(1,length(Cwood_out),timestep_days)
		timestep_days=rep(timestep_days, length.out=length(Cwood_out))
	    }
	    print("...calculating monthly averages for Cwood stocks")
	    # determine the actual daily positions
	    run_day_selector=cumsum(timestep_days)
	    # create needed variables
	    Cwood_agg=array(-9999,dim=length(run_day_selector))
	    Cwood_unc_agg=array(-9999,dim=length(run_day_selector))
	    for (y in seq(1,length(run_day_selector))) {
		pick=Cwood_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]]
		pick2=Cwood_unc_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]]
		Cwood_agg[y]=mean(pick[which(pick != -9999)],na.rm=TRUE)
		Cwood_unc_agg[y]=mean(pick2[which(pick2 != -9999)],na.rm=TRUE)
	    }
	    # now convert the missing values back to -9999
	    Cwood_agg[which(is.na(Cwood_agg))]=-9999
	    Cwood_unc_agg[which(is.na(Cwood_unc_agg))]=-9999
	    # update with new output information
	    Cwood_out=Cwood_agg ; Cwood_unc_out=Cwood_unc_agg
	    # clean up
	    rm(Cwood_agg,Cwood_unc_agg) ; gc()
	} # monthly aggregation etc

	# pass the information back
	return(list(Cwood_stock=Cwood_out,Cwood_stock_unc=Cwood_unc_out))

} # end of function

