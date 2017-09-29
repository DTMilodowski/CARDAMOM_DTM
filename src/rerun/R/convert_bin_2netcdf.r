#!/usr/bin/Rscript

# Script which reads in R binary file and writes the data to netcdf files

# can be DALEC_GSI_BUCKET or DALEC_GSI_DFOL_CWD_FR 
modelname <- "DALEC_GSI_DFOL_CWD_FR"

project <- "BALI_GEMplots_daily"
run <- "001"
site <- ""
lat <- ""
lon <- ""
startyear <- 2011
endyear <- 2016

path2root <- "/home/dmilodow/DataStore_DTM/BALI/CARDAMOM_BALI/"
path2files <- paste(path2root,"projects/",project,"/rerun/",run,"/", sep="")
integer_count_of_sites <- 6

#############################################################################################
# create vector of site ids, 00001, 00002, etc.
vector_of_site_names_or_ids <- 0
for(i in 1:integer_count_of_sites) {
	i_tmp = formatC(i, width = 5, format = "d", flag = "0")
	vector_of_site_names_or_ids[i] = i_tmp
}

# read in binary files
bfile=paste(path2files,vector_of_site_names_or_ids[1],".RData",sep="")

print(paste("Reading in",bfile))
load(bfile)

lai=states_all$lai
gpp=states_all$gpp
evap=states_all$evap

ntsteps=length(lai[1,])
paramsets=length(lai[,1])
print(paste("Number of param sets is ",paramsets,sep=""))

nstats = 5 # mean, median, standard deviation, 25th percentile and 75th percentile

# compute statistics, for each timestep
lai_out=array(0., c(nstats,ntsteps)) 
gpp_out=array(0., c(nstats,ntsteps))
evap_out=array(0., c(nstats,ntsteps)) 

print("Calculating the mean, median, standard deviation, 25th percentile and 75th percentile values for each timestep")
for(i in 1:ntsteps){
	lai_out[1,i] <- mean(lai[,i]) 	
	lai_out[2,i] <- median(lai[,i])
	lai_out[3,i] <- sd(lai[,i]) # standard deviation
	lai_out[4,i] <- quantile(lai[,i], 0.25) # 25th percentile
	lai_out[5,i] <- quantile(lai[,i], 0.75) # 75th percentile

	gpp_out[1,i] <- mean(gpp[,i]) 	
	gpp_out[2,i] <- median(gpp[,i])
	gpp_out[3,i] <- sd(gpp[,i])
	gpp_out[4,i] <- quantile(gpp[,i], 0.25)
	gpp_out[5,i] <- quantile(gpp[,i], 0.75)

	evap_out[1,i] <- mean(evap[,i]) 	
	evap_out[2,i] <- median(evap[,i])
	evap_out[3,i] <- sd(evap[,i])
	evap_out[4,i] <- quantile(evap[,i], 0.25)
	evap_out[5,i] <- quantile(evap[,i], 0.75)
}

#############################################################################################
# create and write to netcdf
library(ncdf4)

f_out=paste(path2files,project, "_", startyear, "_", endyear, "_",modelname,".nc", sep="")

## write out all runs for lai
#f_out_lai=paste(path2files,site, "_weekly_Crop_", startyear, "_", endyear, "_bucket_lai",".nc", sep="")

xvals <- lon
yvals <- lat 
nx <- length(xvals)
ny <- length(yvals)

# define dimensions
londim <- ncdim_def("lon","degrees_east",xvals) 
latdim <- ncdim_def("lat","degrees_north",yvals) 
timedim <- ncdim_def("time","weeks",1:ntsteps,unlim=TRUE)
stats <- ncdim_def("stats","dimenisonless (1-5)",1:nstats) 

var_lai <- ncvar_def("lai", "m2m-2", list(latdim,londim,stats,timedim), longname="Leaf Area Index (LAI) calculated using DALECC-BUCKET")
var_gpp <- ncvar_def("gpp", "gC m-2day-1", list(latdim,londim,stats,timedim), longname="Gross Primary Productivity (GPP) calculated using DALECC-BUCKET")

ncnew <- nc_create(f_out, list(var_lai,var_gpp))

## write lai from all model dalecc runs to file
#parsdim <- ncdim_def("params","nparams",1:paramsets)
#var_lai_all <- ncvar_def("lai", "m2m-2", list(latdim,londim,parsdim,timedim), longname="Leaf Area Index (LAI) calculated using DALECC-BUCKET")
#ncnew_lai <- nc_create(f_out_lai, list(var_lai_all))

print("Writing data to file")
ncvar_put(ncnew,var_lai,lai_out)
ncvar_put(ncnew,var_gpp,gpp_out)

# write lai from all model runs
ncvar_put(ncnew_lai,var_lai_all,lai)
nc_close(ncnew_lai)

print(paste("The file has", ncnew$nvars,"variable(s): lai, gpp"))
print(paste("The file has", ncnew$ndim,"dimension(s): time, stats, lon, lat"))

# Don't forget to close the file
nc_close(ncnew)
print("File created!")






















