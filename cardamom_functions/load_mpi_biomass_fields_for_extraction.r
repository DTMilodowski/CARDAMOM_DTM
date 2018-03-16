
###
## Function to load met data from global field ECMWF data
## subsequently extracted in extract_met_drivers.txt
###

load_mpi_biomass_fields_for_extraction<-function(latlon_in,Cwood_source,Cwood_initial_source) {

    if (Cwood_source == "mpi_biomass" | Cwood_initial_source == "mpi_biomass") {

	# let the user know this might take some time
	print("Loading processed 0.01 degree MPI-Biomass estimates for subsequent sub-setting ...")

	lat_done = FALSE ; missing_years=0 ; keepers=0 ; yrs=1
	# loop for year here
	for (yr in seq(1, length(years_to_do))) {
	    print(paste("... ",round((yr/length(years_to_do))*100,0),"% completed ",Sys.time(),sep=""))

	    if (yr == 1) {
		# first check how many files we have
		for (yrr in seq(1, length(years_to_do))) {
		    if (years_to_do[yrr] == "2010") {keepers=keepers+1} else {missing_years=append(missing_years,years_to_do[yrr])}
		}
	    }
	    # open processed modis files
	    if (years_to_do[yr] == "2010") {
		input_file_1=paste(path_to_biomass,"2014121116258biomass_v3_total.nc",sep="") 
	    } else {
		# this is currently applicable for this year only so if it aint here then we don't care
		input_file_1=" "
	    } 

	    # check to see if file exists if it does then we read it in, if not then we assume its a year we don't have data for and move on
	    if (file.exists(input_file_1) == TRUE) {
		# open the file
		data1=nc_open(input_file_1)

		# as this is a single time stamp we will assume day 1 for the moment
		doy_in=1

		# extract location variables
		if (lat_done == FALSE) {lat=ncvar_get(data1, "latitude") ; long=ncvar_get(data1, "longitude")}

		# read the biomass estimates and uncertainty
		var1=ncvar_get(data1, "biomass_total") ; var2=ncvar_get(data1, "uncertainty_biomass_total")
		# set actual missing data to -9999
		var1[which(is.na(as.vector(var1)))] = -9999 ; var2[which(is.na(as.vector(var2)))] = -9999

		# close files after use
		nc_close(data1)

		# remove additional spatial information
		if (lat_done == FALSE) {
		    if (length(dim(var1)) > 2) {est_obs_in_year=dim(var1)[3]*2} else {est_obs_in_year=1}
		    biomass_hold=array(NA, dim=c(length(long)*length(lat)*est_obs_in_year,keepers))
		    unc_biomass_hold=array(NA, dim=c(length(long)*length(lat)*est_obs_in_year,keepers))
		    biomass_hold[1:length(as.vector(var1)),yrs]=as.vector(var1)
		    unc_biomass_hold[1:length(as.vector(var2)),yrs]=as.vector(var2)
		    doy_out=doy_in
		} else { 
		    biomass_hold[1:length(as.vector(var1)),yrs]=as.vector(var1)
		    unc_biomass_hold[1:length(as.vector(var2)),yrs]=as.vector(var2)
		    doy_out=append(doy_out,doy_in)
		}

		# update flag for lat / long load
		if (lat_done == FALSE) {lat_done = TRUE}
		# keep track of years actually ran
		yrs=yrs+1
	    } # end of does file exist

	} # year loop

	# remove initial value
	missing_years=missing_years[-1]

	# clean up variables
	rm(var1,var2,doy_in,yr,yrr,lat_done,est_obs_in_year) ; gc(reset=TRUE,verbose=FALSE)

	# now remove the ones that are actual missing data
	biomass_hold[which(as.vector(biomass_hold) == -9999)]=NA
	unc_biomass_hold[which(as.vector(unc_biomass_hold) == -9999)]=NA
	# return spatial structure to data
	biomass_out=array(as.vector(biomass_hold), dim=c(length(long),length(lat),length(doy_out)))
	unc_biomass_out=array(as.vector(unc_biomass_hold), dim=c(length(long),length(lat),length(doy_out)))
	max_lat=max(latlon_in[,1])+1.0 ; max_long=max(latlon_in[,2])+1.0
	min_lat=min(latlon_in[,1])-1.0 ; min_long=min(latlon_in[,2])-1.0

	keep_lat=which(lat > min_lat & lat < max_lat)
	keep_long=which(long > min_long & long < max_long)
	unc_biomass_out=unc_biomass_out[,keep_lat,] 
	unc_biomass_out=unc_biomass_out[keep_long,] 
	biomass_out=biomass_out[keep_long,,] 
	biomass_out=biomass_out[,keep_lat]
	# restructure
	biomass_out=array(as.vector(biomass_out), dim=c(dim(biomass_out),length(doy_out)))
	unc_biomass_out=array(as.vector(unc_biomass_out), dim=c(dim(unc_biomass_out),length(doy_out)))

	lat=lat[keep_lat] ; long=long[keep_long]

	# clean up variables
	rm(biomass_hold,unc_biomass_hold,max_lat,min_lat,max_long,min_long,keep_lat,keep_long) ; gc(reset=TRUE,verbose=FALSE)
	
	# output variables; adjust from kgC.m-2 to gC.m-2 
	return(list(biomass_all=(biomass_out*1e3),biomass_all_unc=(unc_biomass_out*1e3),doy_obs=doy_out,lat=lat,long=long,missing_years=missing_years))

    } else if (Cwood_initial_source == "Avitabile"){

	# let the user know this might take some time
	print("Loading processed 0.25 degree Avitabile biomass priors for subsequent sub-setting ...")

	input_file_1=paste(path_to_biomass,"Biomass_stocks_with_lat_long.nc",sep="") 
	# open the file
	data1=nc_open(input_file_1)
	# extract location variables
	lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "long")

	# read the biomass estimates and uncertainty
	biomass=ncvar_get(data1, "Biomass") ; biomass_uncertainty=ncvar_get(data1, "Biomass_Uncertainty")
	# set actual missing data to -9999
	biomass[which(is.na(as.vector(biomass)))] = -9999 ; biomass_uncertainty[which(is.na(as.vector(biomass_uncertainty)))] = -9999

	# close files after use
	nc_close(data1)

	# clean up variables
	gc(reset=TRUE,verbose=FALSE)

	# now remove the ones that are actual missing data
	biomass[which(as.vector(biomass) == -9999)]=NA
	biomass_uncertainty[which(as.vector(biomass_uncertainty) == -9999)]=NA
	# filter around target area
	max_lat=max(latlon_in[,1])+1.0 ; max_long=max(latlon_in[,2])+1.0
	min_lat=min(latlon_in[,1])-1.0 ; min_long=min(latlon_in[,2])-1.0
	keep_lat_min=min(which(lat[1,] > min_lat))
	keep_lat_max=max(which(lat[1,] < max_lat))
	keep_long_min=min(which(long[,1] > min_long))
	keep_long_max=max(which(long[,1] < max_long))
	# remove data outside of target area
	biomass=biomass[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max] 
	biomass_uncertainty=biomass_uncertainty[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max] 
	lat=lat[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max] ; long=long[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]

	# clean up variables
	gc(reset=TRUE,verbose=FALSE)
	
	# output variables; convert MgCha-> gCm-2
	return(list(biomass_all=biomass*1e2,biomass_uncertainty=biomass_uncertainty*1e2,lat=lat,long=long))

    } # if MPI biomass

} # function end

