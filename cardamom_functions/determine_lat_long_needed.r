
###
## Function determines all lat long coordinates needed for grid run mode
###

determine_lat_long_needed<- function(lat,long,resolution,grid_type,remove) {

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
    lat=output$lat ; long=output$long ; rm(output)

    if (use_lcm == "LCM2007") {
	data2=nc_open("/home/lsmallma/WORK/GREENHOUSE/LCM2007/LCM2007_with_lat_long.nc")
    } else if (use_lcm == "CORINE2006") {
	data2=nc_open("/home/lsmallma/WORK/GREENHOUSE/Corine_lcm/Corine2006_at250m_with_lat_long.nc")
    } else if (use_lcm == "CORINE2006_1km") {
	data2=nc_open("/home/lsmallma/WORK/GREENHOUSE/Corine_lcm/Corine2006_at1km_with_lat_long.nc")
    } else if (use_lcm == "forestry_commission") {
	data2=nc_open("/home/lsmallma/data_store/UK_forest_information/UK_forestry_planting_public.nc")
    } else if (use_lcm == "forestry_commission_LCM2007") {
	data2=nc_open("/home/lsmallma/data_store/UK_forest_information/UK_forestry_planting_public_and_private.nc")
    } else if (use_lcm == "forestry_commission_public_private") {
	data2=nc_open("/home/lsmallma/data_store/UK_forest_information/UK_forestry_planting_FC_public_and_private.nc")
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

    # remove the values we don't want
    lat=lat[-remove] ; long=long[-remove]
   
    # output the result
    output=list(lat=lat,long=long)
    rm(lat_lcm,long_lcm,lat,long) ; gc(reset=TRUE, verbose=FALSE)
    return(output)

}
## Use byte compile
determine_lat_long_needed<-cmpfun(determine_lat_long_needed)