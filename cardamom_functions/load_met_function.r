
###
## Load met function
###

# This function assumed ERA-Interim data files are being used

load_met_function<- function (year_to_do,varid,infile_varid,remove_lat,remove_long,path_to_met_source,met_source,wheat) {

      if (met_source == "ECMWF") {
	  # open first file in the sequence
	  m=1
	  input_file_1=paste(path_to_met_source,varid[1],"_ei",year_to_do,"0",m,".nc",sep="")
	  # open netcdf files
	  data1=nc_open(input_file_1) 

	  # read the met drivers
	  var1=ncvar_get(data1, infile_varid[1]) ; var1=var1[,,1:(dim(var1)[3]-1)]

	  # close files after use
	  nc_close(data1)

	  # keep count of time steps
	  t_grid=dim(var1)[3]

	  # filter spatial extent
	  var1=var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
	  # move through time removing the "chaff"
	  for (i in seq(1, t_grid)) {if (i == 1) {tmp=as.vector(var1[,,i])[wheat]} else {tmp=append(tmp,as.vector(var1[,,i])[wheat])} }
	  # now append to the output variable
	  var1_out=tmp

	  # loop through months in the year
	  for (m in seq(2, 12)) {
	      # select correct file names
	      if (m > 9) {
		  input_file_1=paste(path_to_met_source,varid[1],"_ei",year_to_do,m,".nc",sep="")
	      } else {
		  input_file_1=paste(path_to_met_source,varid[1],"_ei",year_to_do,"0",m,".nc",sep="")
	      }
	      # open netcdf files
	      data1=nc_open(input_file_1)
	      # read the met drivers
	      var1=ncvar_get(data1, infile_varid) ; var1=var1[,,1:(dim(var1)[3]-1)]
	      t_grid=t_grid+dim(var1)[3] ; tmp_t=dim(var1)[3]
	      # close files after use
	      nc_close(data1)

	      # update to the correct arrays
	      # filter spatial extent
	      var1=var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
	      # move through time removing the "chaff"
	      for (i in seq(1, tmp_t)) {if (i == 1) {tmp=as.vector(var1[,,i])[wheat]} else {tmp=append(tmp,as.vector(var1[,,i])[wheat])} }
	      # now append to the output variable
	      var1_out=append(var1_out,tmp)

	  } # looping through dataset

	  # clean up
	  rm(var1,tmp,i,m,tmp_t) ; gc(reset=TRUE,verbose=FALSE)

      } else if (met_source == "PRINCETON") {
	  # open first ecmwf file to extract needed information
	  input_file_1=paste(path_to_met_source,varid[1],"_3hourly_",year_to_do,"-",year_to_do,".nc",sep="") 
	  data1=nc_open(input_file_1)
	  # read the met drivers
	  var1=ncvar_get(data1, infile_varid[1]) ; var1=var1[,,1:(dim(var1)[3])]

	  # close files after use
	  nc_close(data1)

	  # keep count of time steps
	  t_grid=dim(var1)[3]

	  # filter spatial extent
	  var1=var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
	  # move through time removing the "chaff"
	  for (i in seq(1, t_grid)) {if (i == 1) {tmp=as.vector(var1[,,i])[wheat]} else {tmp=append(tmp,as.vector(var1[,,i])[wheat])} }
	  # now append to the output variable
	  var1_out=tmp

	  # clean up
	  rm(var1,tmp) ; gc(reset=TRUE,verbose=FALSE)

     } else if (met_source == "CHESS") {
	  # open first file in the sequence
	  m=1
	  input_file_1=paste(path_to_met_source,"chess_",varid[1],"_",year_to_do,"0",m,".nc",sep="") 
	  # open netcdf files
	  data1=nc_open(input_file_1) 

	  # read the met drivers
	  var1=ncvar_get(data1, infile_varid[1]) ; var1=var1[,,1:(dim(var1)[3])]

	  # close files after use
	  nc_close(data1)

	  # keep count of time steps
	  t_grid=dim(var1)[3]
	  # filter spatial extent
	  var1=var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
	  # move through time removing the "chaff"
	  for (i in seq(1, t_grid)) {if (i == 1) {var1_out=as.vector(var1[,,i])[wheat]} else {var1_out=append(var1_out,as.vector(var1[,,i])[wheat])} }

	  # loop through months in the year
	  for (m in seq(2, 12)) {
	      # select correct file names
	      if (m > 9) {
		  input_file_1=paste(path_to_met_source,"chess_",varid[1],"_",year_to_do,m,".nc",sep="") 
	      } else {
		  input_file_1=paste(path_to_met_source,"chess_",varid[1],"_",year_to_do,"0",m,".nc",sep="") 
	      }
	      # open netcdf files
	      data1=nc_open(input_file_1)
	      # read the met drivers
	      var1=ncvar_get(data1, infile_varid) ; var1=var1[,,1:(dim(var1)[3])]
	      t_grid=t_grid+dim(var1)[3] ; tmp_t = dim(var1)[3]
	      # close files after use
	      nc_close(data1)

	      # filter spatial extent
	      var1=var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
	      # move through time removing the "chaff"
	      for (i in seq(1, tmp_t)) {if (i == 1) {tmp=as.vector(var1[,,i])[wheat]} else {tmp=append(tmp,as.vector(var1[,,i])[wheat])} }
	      # now append to the output variable
	      var1_out=append(var1_out,tmp)

	  } # looping through dataset

	  # clean up
	  rm(var1,tmp,i,m,tmp_t) ; gc(reset=TRUE,verbose=FALSE)

      } else if (met_source == "ERA") {
	  # open first file in the sequence
	  m=1
	  input_file_1=paste(path_to_met_source,varid[1],"_",year_to_do,"0",m,".nc",sep="")
	  # open netcdf files
	  data1=nc_open(input_file_1) 

	  # read the met drivers
	  var1=ncvar_get(data1, infile_varid[1]) ; var1=var1[,,1:(dim(var1)[3])]
	  # assign correct error value
	  var1[which(is.na(var1) == TRUE)] = -9999

	  # close files after use
	  nc_close(data1)

	  # keep count of time steps
	  t_grid=dim(var1)[3]

	  # filter spatial extent
	  var1=var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
	  # move through time removing the "chaff"
	  for (i in seq(1, t_grid)) {if (i == 1) {tmp=as.vector(var1[,,i])[wheat]} else {tmp=append(tmp,as.vector(var1[,,i])[wheat])} }
	  # now append to the output variable
	  var1_out=tmp

	  # loop through months in the year
	  for (m in seq(2, 12)) {
	      # select correct file names
	      if (m > 9) {
		  input_file_1=paste(path_to_met_source,varid[1],"_",year_to_do,m,".nc",sep="")
	      } else {
		  input_file_1=paste(path_to_met_source,varid[1],"_",year_to_do,"0",m,".nc",sep="")
	      }
	      # open netcdf files
	      data1=nc_open(input_file_1)
	      # read the met drivers
	      var1=ncvar_get(data1, infile_varid) ; var1=var1[,,1:(dim(var1)[3])]
              # assign correct error value
              var1[which(is.na(var1))] = -9999
	      t_grid=t_grid+dim(var1)[3] ; tmp_t=dim(var1)[3]
	      # close files after use
	      nc_close(data1)

	      # update to the correct arrays
	      # filter spatial extent
	      var1=var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
	      # move through time removing the "chaff"
	      for (i in seq(1, tmp_t)) {if (i == 1) {tmp=as.vector(var1[,,i])[wheat]} else {tmp=append(tmp,as.vector(var1[,,i])[wheat])} }
	      # now append to the output variable
	      var1_out=append(var1_out,tmp)

	  } # looping through dataset

	  # clean up
	  rm(var1,tmp,i,m,tmp_t) ; gc(reset=TRUE,verbose=FALSE)

    } # end data source selection

    # return back to the user
    return(list(var_out=var1_out,t_grid=t_grid))

} # end function
## Use byte compile
load_met_function<-cmpfun(load_met_function)
