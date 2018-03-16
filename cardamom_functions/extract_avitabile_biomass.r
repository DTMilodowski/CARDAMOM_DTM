

extract_avitabile_biomass<- function(timestep_days,spatial_type,resolution,grid_type,latlon_in,Cwood_all) {

	# find the nearest location
	output=closest2d(1,Cwood_all$lat,Cwood_all$long,latlon_in[1],latlon_in[2],2)
	i1=unlist(output)[1] ; j1=unlist(output)[2]

	print(paste("Avitabile Biomass data extracted for current location ",Sys.time(),sep=""))

	# work out number of pixels to average over
	if (spatial_type == "grid") {
	    if (grid_type == "wgs84") {
		# resolution of the product
		product_res = (Cwood_all$lat[1,2]-Cwood_all$lat[1,1])+(Cwood_all$long[2,1]-Cwood_all$long[1,1])
		product_res = product_res * 0.5
		# radius is ceiling of the ratio of the product vs analysis ratio
		radius = ceiling(resolution / product_res)
	    } else if (grid_type == "UK") {
		stop("In correct grid selection for this product")
	    } else {
		stop("have not specified the grid used in this analysis")
	    }
	} else {
	    radius=0
	}

	# work out average areas
	average_i=max(1,(i1-radius)):min(dim(Cwood_all$biomass_all)[1],(i1+radius)) ; average_j=max(1,(j1-radius)):min(dim(Cwood_all$biomass_all)[2],(j1+radius))
	# carry out averaging
	Cwood=mean(Cwood_all$biomass_all[average_i,average_j], na.rm=TRUE)
	Cwood_unc=mean(Cwood_all$biomass_all_unc[average_i,average_j], na.rm=TRUE)

	# convert missing data back to -9999
	Cwood[which(is.na(Cwood))]=-9999.0 ; Cwood_unc[which(is.na(Cwood_unc))]=-9999.0

	# pass the information back
	return(list(Cwood_stock=Cwood,Cwood_stock_unc=Cwood_unc))

} # end of function

