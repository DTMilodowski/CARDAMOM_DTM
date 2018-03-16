
###
## Function which determines how many grid cells 
## are within the defined box
###

## lcm2007_to_ctessel, a function which matches the dominant classifications of the lcm2007 to the appropriate C/D-TESSEL PFT

lcm2007_to_ctessel<- function(input_pft) {
    ## LCM2007 class in order 1-23 = C/D-TESSEL equivalent
    # 1) "Broadleaf forest"        = 5
    # 2) "Needleleaf forest"       = 3
    # 3) "Arable"                  = 1
    # 4) "Improved Grassland"      = 2
    # 5) "Rough Grassland"         = 2
    # 6) "Natural Graassland"      = 2
    # 7) "Calcareous Grassland"    = 2
    # 8) "Acid Grassland"          = 2
    # 9) "Fen Marsh Swamp"         = 13
    #10) "Heather"                 = 2
    #11) "Heather Grassland"       = 2
    #12) "Bog"                     = 13
    #13) "Montane"                 = 9
    #14) "Inland rock"             = 8
    #15) "Saltwater"               = 15
    #16) "Freshwater"              = 14
    #17) "Supra-littoral rock"     = 20
    #18) "Supra-littoral sediment" = 20
    #19) "Littoral rock"           = 20
    #20) "Littoral sediment"       = 20
    #21) "Saltmarsh"               = 20
    #22) "Urban"                   = 19
    #23) "Suburban"                = 19

    # vector of corresponding C/D-TESSEL PFTs in order of the LCM2007 types
    tessel_types=c(5,3,1,2,2,2,2,2,13,2,2,13,9,8,15,14,20,20,20,20,20,19,19)
    # use input LCM2007 cover type to select and return the ctessel PFT
    lcm2007_to_ctessel=tessel_types[input_pft]
    # if location does not have a pft in the lcm2007 make 0 and this will use default values from ECMWF
    if (input_pft==0) {lcm2007_to_ctessel=0}
    # now return needef value
    return(lcm2007_to_ctessel)
}

## Use byte compile
lcm2007_to_ctessel<-cmpfun(lcm2007_to_ctessel)

# corine2006_to_ctessel, a function which matches the dominant classifications of the corine2006 to the appropriate C/D-TESSEL PFT
corine2006_to_ctessel<- function(input_pft) {
    ## Corine2006 class in order 1-43 = C/D-TESSEL equivalent
    # 1) "Continuous urban fabric"         = 0
    # 2) "Discontinuous urban fabric"      = 19
    # 3) "Industrial or urban"             = 0
    # 4) "Road or Rail + associated"       = 0
    # 5) "Port Areas"                      = 0
    # 6) "Airports"                        = 0
    # 7) "Mineral extraction"              = 0
    # 8) "Dump sites"                      = 0
    # 9) "Construction sites"              = 0
    #10) "Green urban areas"               = 19
    #11) "Sport / Leisure facilities"      = 0
    #12) "Non-irrigated arable"            = 1
    #13) "Irrigated arable"                = 10
    #14) "Rice fields"                     = 10
    #15) "Vineyards"                       = 17
    #16) "Fruit tree/berry plantation"     = 18
    #17) "Olive groves"                    = 18
    #18) "Pastures"                        = 2
    #19) "Annual crops + fixed associated" = 1
    #20) "Complex cultivation patterns"    = 1
    #21) "Agriculture with signif natural" = 1
    #22) "Agro-forest areas"               = 3
    #23) "Broadleaf forest"                = 5
    #24) "Coniferous forest"               = 3
    #25) "Mixed forest"                    = 18
    #26) "Natural grassland"               = 2
    #27) "Moors and heathland"             = 2
    #28) "Sclerophyllous veg"              = 17
    #29) "Transitional wood-shrub"         = 19
    #30) "Beaches, Dune, sands"            = 20
    #31) "Bare Rock"                       = 8
    #32) "Sparsely vegetated"              = 11
    #33) "Burnt areas"                     = 0
    #34) "Glaciers and snow"               = 12
    #35) "Inland marsh"                    = 13
    #36) "Peat bog"                        = 13
    #37) "Salt marshes"                    = 20
    #38) "Salines"                         = 20
    #39) "Inter-tidal flats"               = 20
    #40) "Water courses"                   = 14
    #41) "Water bodies"                    = 14
    #42) "Coastal lagoons"                 = 14
    #43) "Estaries"                        = 20

    # vector of corresponding C/D-TESSEL PFTs in order of the Corine2006 types
    tessel_types=c(0,19,0,0,0,0,0,0,0,19,0,1,10,10,17,18,18,2,1,1,1,3,5,3,18,2,2,17,19,20,8,11,0,152,13,13,20,20,20,14,14,14,20)
    # use input Corine2006 cover type to select and return the ctessel PFT
    corine2006_to_ctessel=tessel_types[input_pft]
    # if location does not have a pft in the corine make 0 and this will use default values from ECMWF
    if (input_pft==0) {corine2006_to_ctessel=0}
    # now return needef value
    return(corine2006_to_ctessel)
}
## Use byte compile
corine2006_to_ctessel<-cmpfun(corine2006_to_ctessel)

how_many_points<- function (lat,long,resolution,grid_type) {

    # check input data
    if (length(which(long > 180)) > 0) {stop("Long should be -180 to +180")}

    # generate UK or WGS-84 lat long grid 
    if (grid_type == "UK") {
	output=generate_uk_grid(lat,long,resolution)
    } else if (grid_type=="wgs84") {
	output=generate_wgs84_grid(lat,long,resolution)
    } else {
	stop('have selected invalid grid type, the valid options are "UK" and "wgs84"')
    }
    lat=output$lat ; long=output$long
    lat_dim=output$lat_dim ; long_dim=output$long_dim

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
    } else if (use_lcm == "forestry_commission") {
	data2=nc_open("/home/lsmallma/data_store/UK_forest_information/UK_forestry_planting_public.nc")
	# read pft information
	primary_pft=ncvar_get(data2, "primary_pft")
	secondary_pft=ncvar_get(data2, "secondary_pft")
	tertiary_pft=ncvar_get(data2, "tertiary_pft")
	dims=dim(primary_pft) ; planting_pft=array(-9999, dim=dims)
	# read cover to determine the max cover
	primary_cover=ncvar_get(data2, "primary_cover")
	secondary_cover=ncvar_get(data2, "secondary_cover")
	tertiary_cover=ncvar_get(data2, "tertiary_cover")
	# first make areas with no pft information lost to the analysis
	primary_cover[which(is.na(as.vector(primary_pft)))]=-9999
	secondary_cover[which(is.na(as.vector(secondary_pft)))]=-9999
	tertiary_cover[which(is.na(as.vector(tertiary_pft)))]=-9999
	# might not be the most efficient way of doing things but here we go
	keep_primary=which(as.vector(primary_cover) > 0 & as.vector(primary_cover) > as.vector(secondary_cover) & as.vector(primary_cover) > as.vector(tertiary_cover))
	keep_secondary=which(as.vector(secondary_cover) > 0 & as.vector(secondary_cover) > as.vector(primary_cover) & as.vector(secondary_cover) > as.vector(tertiary_cover))
	keep_tertiary=which(as.vector(tertiary_cover) > 0 & as.vector(tertiary_cover) > as.vector(primary_cover) & as.vector(tertiary_cover) > as.vector(secondary_cover))
	planting_pft[keep_primary]=as.vector(primary_pft)[keep_primary]
	planting_pft[keep_secondary]=as.vector(secondary_pft)[keep_secondary]
	planting_pft[keep_tertiary]=as.vector(tertiary_pft)[keep_tertiary]
	# reconstruct the data
	lcm=array(planting_pft, dim=dims)
	rm(primary_cover,secondary_cover,tertiary_cover,dims,planting_pft,keep_primary,keep_secondary,keep_tertiary)
    } else if (use_lcm == "forestry_commission_LCM2007") {
	data2=nc_open("/home/lsmallma/data_store/UK_forest_information/UK_forestry_planting_public_and_private.nc")
	# read pft information
	primary_pft=ncvar_get(data2, "primary_pft")
	secondary_pft=ncvar_get(data2, "secondary_pft")
	tertiary_pft=ncvar_get(data2, "tertiary_pft")
	dims=dim(primary_pft) ; planting_pft=array(-9999, dim=dims)
	# read cover to determine the max cover
	primary_cover=ncvar_get(data2, "primary_cover")
	secondary_cover=ncvar_get(data2, "secondary_cover")
	tertiary_cover=ncvar_get(data2, "tertiary_cover")
	# first make areas with no pft information lost to the analysis
	primary_cover[which(is.na(as.vector(primary_pft)))]=-9999
	secondary_cover[which(is.na(as.vector(secondary_pft)))]=-9999
	tertiary_cover[which(is.na(as.vector(tertiary_pft)))]=-9999
	# might not be the most efficient way of doing things but here we go
	keep_primary=which(as.vector(primary_cover) > 0 & as.vector(primary_cover) > as.vector(secondary_cover) & as.vector(primary_cover) > as.vector(tertiary_cover))
	keep_secondary=which(as.vector(secondary_cover) > 0 & as.vector(secondary_cover) > as.vector(primary_cover) & as.vector(secondary_cover) > as.vector(tertiary_cover))
	keep_tertiary=which(as.vector(tertiary_cover) > 0 & as.vector(tertiary_cover) > as.vector(primary_cover) & as.vector(tertiary_cover) > as.vector(secondary_cover))
	planting_pft[keep_primary]=as.vector(primary_pft)[keep_primary]
	planting_pft[keep_secondary]=as.vector(secondary_pft)[keep_secondary]
	planting_pft[keep_tertiary]=as.vector(tertiary_pft)[keep_tertiary]
	# reconstruct the data
	lcm=array(planting_pft, dim=dims)
	rm(primary_cover,secondary_cover,tertiary_cover,dims,planting_pft,keep_primary,keep_secondary,keep_tertiary)
    } else if (use_lcm == "forestry_commission_public_private") {
	data2=nc_open("/home/lsmallma/data_store/UK_forest_information/UK_forestry_planting_FC_public_and_private.nc")
	# read pft information
	primary_pft=ncvar_get(data2, "primary_pft")
	secondary_pft=ncvar_get(data2, "secondary_pft")
	tertiary_pft=ncvar_get(data2, "tertiary_pft")
	dims=dim(primary_pft) ; planting_pft=array(-9999, dim=dims)
	# read cover to determine the max cover
	primary_cover=ncvar_get(data2, "primary_cover")
	secondary_cover=ncvar_get(data2, "secondary_cover")
	tertiary_cover=ncvar_get(data2, "tertiary_cover")
	# first make areas with no pft information lost to the analysis
	primary_cover[which(is.na(as.vector(primary_pft)))]=-9999
	secondary_cover[which(is.na(as.vector(secondary_pft)))]=-9999
	tertiary_cover[which(is.na(as.vector(tertiary_pft)))]=-9999
	# might not be the most efficient way of doing things but here we go
	keep_primary=which(as.vector(primary_cover) > 0 & as.vector(primary_cover) > as.vector(secondary_cover) & as.vector(primary_cover) > as.vector(tertiary_cover))
	keep_secondary=which(as.vector(secondary_cover) > 0 & as.vector(secondary_cover) > as.vector(primary_cover) & as.vector(secondary_cover) > as.vector(tertiary_cover))
	keep_tertiary=which(as.vector(tertiary_cover) > 0 & as.vector(tertiary_cover) > as.vector(primary_cover) & as.vector(tertiary_cover) > as.vector(secondary_cover))
	planting_pft[keep_primary]=as.vector(primary_pft)[keep_primary]
	planting_pft[keep_secondary]=as.vector(secondary_pft)[keep_secondary]
	planting_pft[keep_tertiary]=as.vector(tertiary_pft)[keep_tertiary]
	# reconstruct the data
	lcm=array(planting_pft, dim=dims)
	rm(primary_cover,secondary_cover,tertiary_cover,dims,planting_pft,keep_primary,keep_secondary,keep_tertiary)
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
	hi_veg_frac=as.vector(hi_veg_frac) ; low_veg_frac=as.vector(low_veg_frac)
	hi_veg_type=as.vector(hi_veg_type) ; low_veg_type=as.vector(low_veg_type)
	lcm=hi_veg_type ; lcm[which(low_veg_frac > hi_veg_frac)]=low_veg_type[which(low_veg_frac > hi_veg_frac)]
	lat_lcm=ncvar_get(data2, "lat")
	long_lcm=ncvar_get(data2, "lon")
	long_lcm[which(long_lcm > 180)] = long_lcm[which(long_lcm > 180)]-360
    } else {
	stop("bugger no land cover option found / set")
    }
    # download location data
    if (use_lcm != "ECMWF") {
	lat_lcm=ncvar_get(data2,"lat")
	long_lcm=ncvar_get(data2,"long")
    }

    # house keeping
    nc_close(data2)

    # raw total pixels
    print(paste("Raw pixel total is ",length(lat)," next filter for water bodies"))

    # find locations
    if (use_parallel) {
	cl <- makeCluster(numWorkers, type = "PSOCK")    
	# load R libraries in cluster
	clusterExport(cl,"load_r_libraries") ; clusterEvalQ(cl, load_r_libraries())
	if (use_lcm == "ECMWF") {
	    output=parLapply(cl,1:length(lat),fun=closest2d,lat=lat_lcm,long=long_lcm,lat_in=lat,long_in=long,nos_dim=1) 
	    stopCluster(cl)
	    # extract the i,j values seperately
	    output_i=unlist(output)
	} else {
	    output=parLapply(cl,1:length(lat),fun=closest2d,lat=lat_lcm,long=long_lcm,lat_in=lat,long_in=long,nos_dim=2) 
	    stopCluster(cl)
	    # extract the i,j values seperately
	    output_i=unlist(output)[which((1:length(unlist(output))*0.5) != floor(1:length(unlist(output))*0.5))]
	    output_j=unlist(output)[which((1:length(unlist(output))*0.5) == floor(1:length(unlist(output))*0.5))]
	}

     } else {
	if (use_lcm == "ECMWF") {
	    output=lapply(1:length(lat),FUN=closest2d,lat=lat_lcm,long=long_lcm,lat_in=lat,long_in=long,nos_dim=1) 
	    # extract the i,j values seperately
	    output_i=unlist(output)
	} else {
	    output=lapply(1:length(lat),FUN=closest2d,lat=lat_lcm,long=long_lcm,lat_in=lat,long_in=long,nos_dim=2) 
	    # extract the i,j values seperately
	    output_i=unlist(output)[which((1:length(unlist(output))*0.5) != floor(1:length(unlist(output))*0.5))]
	    output_j=unlist(output)[which((1:length(unlist(output))*0.5) == floor(1:length(unlist(output))*0.5))]
	}
    }

    print("Generating land sea mask")

    # load global shape file for land sea mask
    landmask=shapefile("./cardamom_functions/global_map/ne_10m_admin_0_countries.shx")
    # just to be sure enforce the projection to WGS-84
    landmask=spTransform(landmask,CRS("+init=epsg:4326"))
    # extract lat/long information for subsetting
    landsea_long=coordinates(landmask)[,1] ; landsea_lat=coordinates(landmask)[,2]
    # filter to constrain the target area
    filter=which(landsea_lat < max(lat)+3) 
    landsea_long=landsea_long[filter] ; landsea_lat=landsea_lat[filter] ; landmask=landmask[filter,]
    filter=which(landsea_lat > min(lat)-3) 
    landsea_long=landsea_long[filter] ; landsea_lat=landsea_lat[filter] ; landmask=landmask[filter,]
    filter=which(landsea_long < max(long)+3) 
    landsea_long=landsea_long[filter] ; landsea_lat=landsea_lat[filter] ; landmask=landmask[filter,]
    filter=which(landsea_long > min(long)-3) 
    landsea_long=landsea_long[filter] ; landsea_lat=landsea_lat[filter] ; landmask=landmask[filter,]
    # generate mask at the resolution of the analysis - ok this wont be an exact match 
    # because of the slightly larger area being selected here but it will be close enough for maintaining resolution
    landsea=raster(ncols=lat_dim,nrows=long_dim)
    # update the extent inforation for the raster
    extent(landsea)<-extent(landmask)
    landsea=rasterize(landmask,landsea,coordinates(landmask))
    # extract lat/long information for the raster now
    landsea_long=coordinates(landsea)[,1] ; landsea_lat=coordinates(landsea)[,2]
    # now convert the raster into a readable array
    landsea=as.vector(landsea) ; landsea[which(is.na(landsea) == FALSE)]=1 ; landsea[which(is.na(landsea))]=0

#    # load the land sea mask 20 km resolution
#    load("./cardamom_functions/landmask20km.rda")
#    # 1 = land, 0 = sea
#    #image(landmask20km[1]) ; length(which(landmask20km[[1]] == 1)) ; length(which(landmask20km[[1]] == 0))
#    # force a change from the Robinson projection to the WGS-84 lat/long
#    landmask20km=spTransform(landmask20km,CRS("+init=epsg:4326"))
#    # extract the lat long coordinate systems
#    landsea_long=coordinates(landmask20km)[,1] ; landsea_lat=coordinates(landmask20km)[,2]
#    # extract the land sea mask information 
#    landsea=landmask20km[[1]]
#    # filter to constrain the target area
#    filter=which(landsea_lat < max(lat)+0.5) 
#    landsea_long=landsea_long[filter] ; landsea_lat=landsea_lat[filter] ; landsea=landsea[filter]
#    filter=which(landsea_lat > min(lat)-0.5) 
#    landsea_long=landsea_long[filter] ; landsea_lat=landsea_lat[filter] ; landsea=landsea[filter]
#    filter=which(landsea_long < max(long)+0.5) 
#    landsea_long=landsea_long[filter] ; landsea_lat=landsea_lat[filter] ; landsea=landsea[filter]
#    filter=which(landsea_long > min(long)-0.5) 
#    landsea_long=landsea_long[filter] ; landsea_lat=landsea_lat[filter] ; landsea=landsea[filter]

    # find locations
    if (use_parallel) {
	cl <- makeCluster(numWorkers, type = "PSOCK")    
	# load R libraries in cluster
	clusterExport(cl,"load_r_libraries") ; clusterEvalQ(cl, load_r_libraries())
	output=parLapply(cl,1:length(lat),fun=closest2d,lat=landsea_lat,long=landsea_long,lat_in=lat,long_in=long,nos_dim=1) 
	stopCluster(cl)
	# extract the i,j values seperately
	output_k=unlist(output)
     } else {
	output=lapply(1:length(lat),FUN=closest2d,lat=landsea_lat,long=landsea_long,lat_in=lat,long_in=long,nos_dim=1) 
	# extract the i,j values seperately
	output_k=unlist(output)
    }

    # now select only the bits of the vector we want
    landsea = landsea[output_k]

    # find out the pft now
    remove=0 ; pft_keep=0
    # now iterate through the sites
    for (pft in seq(1, length(lat))) {
	# update the user, but only sometimes
	if (pft%%2000 == 0 | pft < 500) {print(paste("Ocean filter ",round((pft/length(lat))*100,0),"% complete",sep=""))}
	# convert incoming pft to common values (in this case CTESSEL)
	if (use_lcm == "LCM2007") {
	    new_pft=lcm2007_to_ctessel(lcm[output_i[pft],output_j[pft]])
	} else if (use_lcm == "CORINE2006") {
	    new_pft=corine2006_to_ctessel(lcm[output_i[pft],output_j[pft]])
	} else if (use_lcm == "CORINE2006_1km") {
	    new_pft=corine2006_to_ctessel(lcm[output_i[pft],output_j[pft]])
	} else if (use_lcm == "forestry_commission" | use_lcm == "forestry_commission_LCM2007" | use_lcm == "forestry_commission_public_private") {
	    new_pft=lcm[output_i[pft],output_j[pft]]
	    if (new_pft < 0 | length(new_pft) == 0) {new_pft = 0}
	} else if (use_lcm == "ECMWF") {
	    new_pft = lcm[output_i[pft]]
	}
	# now exclude if not a land site
	if (new_pft == 0 | new_pft == 14 | new_pft == 15) {
	    remove=append(remove,pft)
	} else {
	    pft_keep=append(pft_keep,new_pft)
	}
    } # sites for loop

    # remove initial value
    pft_keep=pft_keep[-1] ; remove=remove[-1] ; rm(lat_lcm,long_lcm)

    # now remove water pixels from analysis
    lat=lat[-remove] ; long=long[-remove]

    # generate the site names
    b=1 ; sites=rep("NA",times=(length(lat)+length(remove)))
    for (n in seq(1, (length(lat)+length(remove)))) {
	if (n%%2000 == 0 | n < 500) {print(paste("Have generated ",round((b/(length(lat)+length(remove)))*100,0),"% of site IDs" ,sep=""))}
	# we want the numbers to match their location within the domain not of the land pixels.
	# this is needed for easy reconstruction later
	sites[b]=sprintf('%05i',n) ; b=b+1
    }
    # now sub selecting for the ones we want
    sites=sites[-remove]

    # product is the number of points
    return(list(nosites=length(lat),waterpixels=remove,landsea=landsea,ctessel_pft=pft_keep,lat_dim=lat_dim,long_dim=long_dim,sites=sites))
    # clean up
    gc(reset=TRUE, verbose=FALSE)

}

## Use byte compile
how_many_points<-cmpfun(how_many_points)