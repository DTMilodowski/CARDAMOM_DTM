
###
## Function which determines how many grid cells 
## are within the defined box
###

find_pft<- function (lat,long) {

    # check input data
    if (length(which(long > 180)) > 0) {stop("Long should be -180 to +180")}

    # now work out how many of these are land points
    # determine whether using LCM 2007 or default ECMWF land cover map
    if (use_lcm == "LCM2007") {
	data2=nc_open("/home/lsmallma/WORK/GREENHOUSE/LCM2007/LCM2007_with_lat_long.nc")
	lcm=ncvar_get(data2,"LCM2007")
    } else if (use_lcm == "CORINE2006") {
	data2=nc_open("/home/lsmallma/WORK/GREENHOUSE/Corine_lcm/Corine2006_at250m_with_lat_long.nc")
	lcm=ncvar_get(data2,"Corine2006")
    } else if (use_lcm == "CORINE2006_1km") {
	data2=nc_open("/home/lsmallma/WORK/GREENHOUSE/Corine_lcm/Corine2006_at1km_with_lat_long.nc")
	lcm=ncvar_get(data2,"Corine2006")
    } else if (use_lcm == "ECMWF") {
	# load global surfclim file and info file for surfclim
	data2=nc_open("/home/lsmallma/gcel/ECMWF_ERA_MET_2000_2012/surfclim_all.nc")
	# extract high vegetation cover fraction
	hi_veg_frac=ncvar_get(data2, "CVH")
	# extract low vegetation cover fraction
	low_veg_frac=ncvar_get(data2, "CVL")
	# extract high vegetation type
	hi_veg_type=ncvar_get(data2, "TVH")
	# extract low vegetation type
	low_veg_type=ncvar_get(data2, "TVL")
    } else {
	stop("bugger no land cover option found / set")
    }

    if (use_lcm == "ECMWF") {
	# download location data
	lat_lcm=ncvar_get(data2,"lat")
	long_lcm=ncvar_get(data2,"lon")
    } else {
	# download location data
	lat_lcm=ncvar_get(data2,"lat")
	long_lcm=ncvar_get(data2,"long")
    }

    # convert input data long to conform to what we need
    check1=which(long_lcm > 180) ; if (length(check1) > 0) { long_lcm[check1]=long_lcm[check1]-360 }

    # house keeping
    nc_close(data2)

    # find out the pft now
    pft_keep=0
    # find locations
    if (use_parallel) {
	if (use_lcm == "ECMWF") {
	    cl <- makeCluster(min(length(lat),numWorkers), type = "PSOCK")    
	    output=parLapply(cl,1:length(lat),fun=closest2d,lat=lat_lcm,long=long_lcm,lat_in=lat,long_in=long,nos_dim=1) 
	    stopCluster(cl)
	    output_i=unlist(output)
	} else {
	    cl <- makeCluster(min(length(lat),numWorkers), type = "PSOCK")    
	    output=parLapply(cl,1:length(lat),fun=closest2d,lat=lat_lcm,long=long_lcm,lat_in=lat,long_in=long,nos_dim=2) 
	    stopCluster(cl)
	    # extract the i,j values seperately
	    output_i=unlist(output)[which((1:length(unlist(output))*0.5) != floor(1:length(unlist(output))*0.5))]
	    output_j=unlist(output)[which((1:length(unlist(output))*0.5) == floor(1:length(unlist(output))*0.5))]
	}
     } else {
	if (use_lcm == "ECMWF") {
	    # search for nearest matches
	    output=lapply(1:length(lat),FUN=closest2d,lat=lat_lcm,long=long_lcm,lat_in=lat,long_in=long,nos_dim=1) 
	    output_i=unlist(output)[1]
	} else {
	    output=lapply(1:length(lat),FUN=closest2d,lat=lat_lcm,long=long_lcm,lat_in=lat,long_in=long,nos_dim=2) 
	    # extract the i,j values seperately
	    output_i=unlist(output)[which((1:length(unlist(output))*0.5) != floor(1:length(unlist(output))*0.5))]
	    output_j=unlist(output)[which((1:length(unlist(output))*0.5) == floor(1:length(unlist(output))*0.5))]
	}
    }

    # return the longitude back to normal
    if (length(check1) > 0) { long_lcm[check1]=long_lcm[check1]+360 }

    # now iterate throught the list
    for (pft in seq(1, length(lat))) {
	if (use_lcm == "LCM2007") {
	    new_pft=lcm2007_to_ctessel(lcm[output_i[pft],output_j[pft]])
	} else if (use_lcm == "CORINE2006") {
	    new_pft=corine2006_to_ctessel(lcm[output_i[pft],output_j[pft]])
	} else if (use_lcm == "CORINE2006_1km") {
	    new_pft=corine2006_to_ctessel(lcm[output_i[pft],output_j[pft]])
	} else if (use_lcm == "ECMWF") {
	    # if using ECMWF then no change needed, pass direct
	    # while selecting for dominant value
	    if (low_veg_frac[output_i[pft]] > hi_veg_frac[output_i[pft]]) {
		new_pft=low_veg_type[output_i[pft]]
	    } else {
		new_pft=hi_veg_type[output_i[pft]]
	    }
	}
	# now exclude if not a land site
	if (new_pft == 0 | new_pft == 14 | new_pft == 15) {
	    stop(paste('lat/long provides has been identified as a water body... position ',pft," of ",length(lat),sep=""))
	} else {
	    pft_keep=append(pft_keep,new_pft)
	}
    }
    # remove initial value
    pft_keep=pft_keep[-1]

    # product is the number of points
    return(list(ctessel_pft=pft_keep))
    # clean up
    gc(reset=TRUE, verbose=FALSE)

} # end of find_pft

## Use byte compile
find_pft<-cmpfun(find_pft)