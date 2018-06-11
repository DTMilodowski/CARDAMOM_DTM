

generate_wgs84_grid<-function(lat,long,resolution) {
    # load locations into the dataframe
    pt=data.frame(x=long,y=lat)
    # create coordinates
    coordinates(pt)=~x+y
    # impose lat/long system (WGS-84)
    proj4string(pt)=CRS("+init=epsg:4326")
    # transform to ORGB system
    # extract box edges in (m)
    bottom_left_x=coordinates(pt)[1,1]
    bottom_left_y=coordinates(pt)[1,2]
    top_right_x=coordinates(pt)[2,1]
    top_right_y=coordinates(pt)[2,2]

    # inputs for northings and eastings (m)
    N_to_do=seq(bottom_left_y,top_right_y, by=resolution) 
    E_to_do=seq(bottom_left_x,top_right_x, by=resolution) 

    # loop to fill in the blacks
    lat_dim=length(N_to_do) ; long_dim=length(E_to_do)
    N_to_do=rep(N_to_do, each=long_dim)
    E_to_do=rep(E_to_do, times=lat_dim)

    # convert these back to the lat long grid
    pt=data.frame(x=E_to_do,y=N_to_do)
    coordinates(pt)=~x+y
    # impose WGS-84
    proj4string(pt)=CRS("+init=epsg:4326")

    # extract needed values
    lat=coordinates(pt)[,2] ; long=coordinates(pt)[,1]

    # clean up

    # return
    output=list(lat=lat,long=long,lat_dim=lat_dim,long_dim=long_dim)
    rm(lat,long,lat_dim,long_dim,pt,E_to_do,N_to_do,bottom_left_x,bottom_left_y,top_right_x,top_right_y) ; gc(reset=TRUE,verbose=FALSE)	
    return(output)

} #end function

## Use byte compile
generate_wgs84_grid<-cmpfun(generate_wgs84_grid)