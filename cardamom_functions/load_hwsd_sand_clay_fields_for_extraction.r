
###
## Function to load met data from global field ECMWF data
## subsequently extracted in extract_met_drivers.txt
###

load_hwsd_sand_clay_fields_for_extraction<-function(latlon_in,sand_clay_source) {
    
    if (sand_clay_source == "HWSD") {
    
	# let the user know this might take some time
	print("Loading processed HWSD sand clay fields for subsequent sub-setting ...")

	# open processed modis files
	input_file_1=paste(path_to_hwsd_Csom,"/HWSD_sand_clay_with_lat_long.nc",sep="") 
	data1=nc_open(input_file_1)

	# extract location variables
	lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "long")
	# read the modis lai drivers
	hwsd_top_sand=ncvar_get(data1, "HWSD_top_sand") ; hwsd_top_clay=ncvar_get(data1, "HWSD_top_clay")
	hwsd_bot_sand=ncvar_get(data1, "HWSD_bot_sand") ; hwsd_bot_clay=ncvar_get(data1, "HWSD_bot_clay")

	if (length(dim(latlon_in)) > 1) {
	    max_lat=max(latlon_in[,1])+0.5 ; max_long=max(latlon_in[,2])+0.5
	    min_lat=min(latlon_in[,1])-0.5 ; min_long=min(latlon_in[,2])-0.5
	} else {
	    max_lat=max(latlon_in[1])+0.5 ; max_long=max(latlon_in[2])+0.5
	    min_lat=min(latlon_in[1])-0.5 ; min_long=min(latlon_in[2])-0.5	    
	}
	keep_lat=which(lat[1,] > min_lat & lat[1,] < max_lat)
	keep_long=which(long[,1] > min_long & long[,1] < max_long)
	hwsd_top_sand=hwsd_top_sand[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	hwsd_bot_sand=hwsd_bot_sand[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	hwsd_top_clay=hwsd_top_clay[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	hwsd_bot_clay=hwsd_bot_clay[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]    
	lat=lat[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)] 
	long=long[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	# close files after use
	nc_close(data1)

	# clean
	rm(keep_lat,keep_long,max_lat,max_long,min_lat,min_long) ; gc(reset=TRUE,verbose=FALSE)

	# output variables
	return(list(top_sand=hwsd_top_sand,top_clay=hwsd_top_clay,bot_sand=hwsd_bot_sand,bot_clay=hwsd_bot_clay,lat=lat,long=long))

    } else {
	# output variables
	return(list(top_sand=40,top_clay=15,bot_sand=40,bot_clay=15,lat=-9999,long=-9999))
    }

} # function end

