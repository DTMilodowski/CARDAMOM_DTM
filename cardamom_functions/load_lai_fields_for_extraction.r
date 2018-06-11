
###
## Function to load met data from global field ECMWF data
## subsequently extracted in extract_met_drivers.txt
###

load_lai_fields_for_extraction<-function(latlon_in,lai_source,years_to_load) {

    if (lai_source == "MODIS") {

	# let the user know this might take some time
	print("Loading processed lai fields for subsequent sub-setting ...")

	lat_done = FALSE ; missing_years=0 ; keepers=0 ; yrs=1
	# loop for year here
	for (yr in seq(1, length(years_to_load))) {
	    print(paste("... ",round((yr/length(years_to_load))*100,0),"% completed ",Sys.time(),sep=""))

	    if (yr == 1) {
		# first check how many files we have
		for (yrr in seq(1, length(years_to_load))) {
		    # open processed modis files
		    input_file_1=paste(path_to_lai,"modis_lai_with_lat_long_",years_to_load[yrr],".nc",sep="") 
		    if (file.exists(input_file_1) == TRUE) {keepers=keepers+1} else {missing_years=append(missing_years,years_to_load[yrr])}
		}
	    }
	    # open processed modis files
	    input_file_1=paste(path_to_lai,"modis_lai_with_lat_long_",years_to_load[yr],".nc",sep="")

	    # check to see if file exists if it does then we read it in, if not then we assume its a year we don't have data for and move on
	    if (file.exists(input_file_1) == TRUE) {
		# open the file
		data1=nc_open(input_file_1)

		# get timing variable
		doy_in=ncvar_get(data1,"doy")

		# extract location variables
		if (lat_done == FALSE) {
		    lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "long")
		    # restrict the spatial extent based on latlong ranges provided
		    remove_lat=intersect(which(lat < (max(latlon_in[,1])+2)),which(lat > (min(latlon_in[,1])-2)))
		    remove_long=intersect(which(long < (max(latlon_in[,2])+2)),which(long > (min(latlon_in[,2])-2)))
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
		}

		# read the modis lai drivers
		var1=ncvar_get(data1, "modis_lai")
		# reduce spatial cover to the desired area only
		var1=var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
		# set actual missing data to -9999
		var1[which(is.na(as.vector(var1)))] = -9999

		# close files after use
		nc_close(data1)

		# remove additional spatial information
		if (lat_done == FALSE) {
		    lai_hold=array(NA, dim=c(long_dim*lat_dim*48,keepers))
		    lai_hold[1:length(as.vector(var1)),yrs]=as.vector(var1)
		    doy_out=doy_in
		} else { 
		    lai_hold[1:length(as.vector(var1)),yrs]=as.vector(var1)
		    doy_out=append(doy_out,doy_in)
		}

		# update flag for lat / long load
		if (lat_done == FALSE) {lat_done = TRUE}
		# keep track of years actually ran
		yrs=yrs+1
		# clean up allocated memeory
		rm(var1) ; gc()
	    } # end of does file exist

	} # year loop

        # Sanity check for LAI
        if (lat_done == FALSE) {stop('No LAI information could be found...')}

	# remove initial value
	missing_years=missing_years[-1]

	# check which ones are NA because I made them up
	not_na=is.na(as.vector(lai_hold))
        not_na=which(not_na == FALSE)

	# now remove the ones that are actual missing data
	lai_hold[which(as.vector(lai_hold) == -9999)]=NA
	# return spatial structure to data
	lai_out=array(as.vector(lai_hold)[not_na], dim=c(long_dim,lat_dim,length(doy_out)))

	# output variables
	lai_all=list(lai_all=lai_out,doy_obs=doy_out,lat=lat,long=long,missing_years=missing_years)
	# clean up variables
	rm(doy_in,lai_hold,not_na,lai_out,doy_out,lat,long,missing_years) ; gc(reset=TRUE,verbose=FALSE)
	return(lai_all)

    } # if MODIS

} # function end
## Use byte compile
load_lai_fields_for_extraction<-cmpfun(load_lai_fields_for_extraction)
