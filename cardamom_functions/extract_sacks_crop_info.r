####
## Function to extract the HWSD soil carbon estimate for initial conditions
####

extract_sacks_crop_info<- function(spatial_type,resolution,grid_type,latlon_in,crop_man_all) {

	# convert input data long to conform to what we need
	check1=which(crop_man_all$long > 180) ; if (length(check1) > 0) { crop_man_all$long[check1]=crop_man_all$long[check1]-360 }
	
	# find the nearest location
	output=closest2d(1,crop_man_all$lat,crop_man_all$long,latlon_in[1],latlon_in[2],3)
	#i1=unlist(output)[2] ; j1=unlist(output)[1]
	i1=unlist(output)[1] ; j1=unlist(output)[2]
	print(paste("Crop management data extracted for current location ",Sys.time(),sep=""))

	# return long to 0-360
	if (length(check1) > 0) { crop_man_all$long[check1]=crop_man_all$long[check1]+360 }

	# work out number of pixels to average over
#	print("NOTE all Csom values are a minimum average of 9 pixels (i.e centre+1 )")
	if (spatial_type == "grid") {
	    if (grid_type == "wgs84") {
		# calculate pixel area and convert from m2 -> km2
		area=calc_pixel_area(latlon_in[1],latlon_in[2],resolution)*1e-6
		radius=max(0,ceiling(sqrt(area)*0.5))
	    } else if (grid_type == "UK") {
		radius=max(0,floor(1*resolution*1e-3*0.5))
	    } else {
		stop("have not specified the grid used in this analysis")
	    }
	} else {
	    radius=0
	}
	answer=NA ; plant = -9999 ; harvest = -9999 ; plant_range = -9999 ; harvest_range = -9999
	while (is.na(answer) == TRUE) {
	    # work out average areas
	    average_i=max(1,i1-radius):min(length(crop_man_all$long),i1+radius) ; average_j=max(1,j1-radius):min(length(crop_man_all$lat),j1+radius)
	    # carry out averaging
	    tmp=crop_man_all$plant[average_i,average_j] ; tmp[which(tmp == -9999)]=NA
	    plant=mean(tmp, na.rm=TRUE) 
	    tmp=crop_man_all$harvest[average_i,average_j] ; tmp[which(tmp == -9999)]=NA
	    harvest=mean(tmp, na.rm=TRUE) 
#	    tmp=crop_man_all$plant_range[average_i,average_j] ; tmp[which(tmp == -9999)]=NA
#	    plant_range=mean(tmp, na.rm=TRUE) 
#	    tmp=crop_man_all$harvest_range[average_i,average_j] ; tmp[which(tmp == -9999)]=NA
#	    harvest_range=mean(tmp, na.rm=TRUE) 
	    # error checking
	    if (is.na(plant) | plant == 0 | harvest == 0 | is.na(harvest)) {radius=radius+1 ; answer=NA} else {answer=0}
	}
	print(paste("NOTE Sowing and Harvest dates averaged over a pixel radius (i.e. centre + radius) of ",radius," points",sep=""))
	# pass the information back
	return(list(plant=plant,plant_range=plant_range,harvest=harvest,harvest_range=harvest_range))

} # end of function