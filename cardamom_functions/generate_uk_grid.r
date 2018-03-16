

generate_uk_grid<-function(lat,long,resolution) {
    # load locations into the dataframe
    pt=data.frame(x=long,y=lat)
    # create coordinates
    coordinates(pt)=~x+y
    # impose lat/long system (WGS-84)
    proj4string(pt)=CRS("+init=epsg:4326")
    # transform to ORGB system
    # +proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs 
    pt=spTransform(pt,CRS("+init=epsg:27700"))
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
    # impose ORGB
    proj4string(pt)=CRS("+init=epsg:27700")
    # reproject to WGS-84
    pt=spTransform(pt,CRS("+init=epsg:4326"))
    # extract needed values
    lat=coordinates(pt)[,2] ; long=coordinates(pt)[,1]

    # return
    output=list(lat=lat,long=long,lat_dim=lat_dim,long_dim=long_dim)
    rm(lat,long,lat_dim,long_dim,pt,E_to_do,N_to_do,bottom_left_x,bottom_left_y,top_right_x,top_right_y) ; gc(reset=TRUE,verbose=FALSE)	
    return(output)

} #end function

## Use byte compile
generate_uk_grid<-cmpfun(generate_uk_grid)