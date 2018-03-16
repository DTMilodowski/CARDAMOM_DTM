

###
## Function to determine the area of a pixel at any location in meter squared
###

## function to convert degrees to radians
deg2rad<-function(degree_in) {return(degree_in*(pi/180))}

calc_pixel_area<-function(latitude,longitude,resolution) {

    # resolution in degrees
    # latitude (-90,90) degrees
    # longitude (-180,180) degrees

    # mean earth radius (m)
    R = 6371e3 

    pixel_area=R**2 * ( deg2rad(longitude+resolution/2.)-deg2rad(longitude-resolution/2.) ) * (sin(deg2rad(latitude+resolution/2.))-sin(deg2rad(latitude-resolution/2.)))

    # return to user in meters
    return(pixel_area)

}
