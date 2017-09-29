#!/usr/bin/Rscript

# Script which reads in R binary file and writes the data to netcdf files

# can be DALEC_GSI_BUCKET or DALEC_GSI_DFOL_CWD_FR 
modelname <- "DALEC_GSI_BUCKET"

site <- "Grignon"

# Aurade - 43.5496, 1.1061  
# Grignon - 48.8442, 1.9519 
# Gebesee - 51.1001, 10.9143 
# Oensingen - 47.2863, 7.7343
# Lonzee - 50.5516, 4.7461   
# Avignon - 43.9161, 4.8781 
# Klingenberg - 50.8929, 13.5225 
# Lamasquere - 43.4965, 1.2379
# Risbyholm - 55.5303, 12.0972

if (site=="Aurade"){
	startyear <- 2005
	endyear <- 2010
	lat <- 43.5496
	lon <- 1.1061  
} else if (site == "Avignon"){
	startyear <- 2004
	endyear <- 2007
	lat <- 43.9161
	lon <- 4.8781  
} else if (site == "Gebesee"){
	startyear <- 2002
	endyear <- 2010
	lat <- 51.1001
	lon <- 10.9143
} else if (site == "Grignon"){
	startyear <- 2004
	endyear <- 2007
	lat <- 48.8442
	lon <- 1.9519
} else if (site == "Klingenberg"){
	startyear <- 2004
	endyear <- 2012
	lat <- 50.8929
	lon <- 13.5225
} else if (site == "Lamasquere"){
	startyear <- 2005
	endyear <- 2010
	lat <- 43.4965
	lon <- 1.2379
} else if (site == "Lonzee"){
	startyear <- 2004
	endyear <- 2011
	lat <- 50.5516
	lon <- 4.7461
} else if (site == "Oensingen"){
	startyear <- 2004
	endyear <- 2011
	lat <- 47.2863
	lon <- 7.7343
} else if (site == "Risbyholm"){
	startyear <- 2004
	endyear <- 2008
	lat <- 55.5303
	lon <- 12.0972
}

path2files <- paste(site, "/DALECc_output/", sep="")

integer_count_of_sites <- 1

#############################################################################################
# create vector of site ids, 00001, 00002, etc.
vector_of_site_names_or_ids <- 0
for(i in 1:integer_count_of_sites) {
	i_tmp = formatC(i, width = 5, format = "d", flag = "0")
	vector_of_site_names_or_ids[i] = i_tmp
}

# read in binary files
bfile=paste(path2files,site, "_weekly_Crop_", startyear, "_", endyear, "_bucket_",vector_of_site_names_or_ids[1],".RData",sep="")

print(paste("Reading in ",bfile,sep=""))
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

f_out=paste(path2files,site, "_weekly_Crop_", startyear, "_", endyear, "_bucket",".nc", sep="")

# write out all runs for lai
f_out_lai=paste(path2files,site, "_weekly_Crop_", startyear, "_", endyear, "_bucket_lai",".nc", sep="")

print("Removing old files")
if (file.exists(paste(site,"/DALECc_output/",f_out,sep=""))) file.remove(paste(site,"/DALECc_output/",f_out,sep=""))

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
var_evap <- ncvar_def("et", "kgH2O m-2day-1", list(latdim,londim,stats,timedim), longname="Evapotranspiration (ET) calculated using DALECC-BUCKET")

ncnew <- nc_create(f_out, list(var_lai,var_gpp,var_evap))

# write lai from all model dalecc runs to file
parsdim <- ncdim_def("params","nparams",1:paramsets)
var_lai_all <- ncvar_def("lai", "m2m-2", list(latdim,londim,parsdim,timedim), longname="Leaf Area Index (LAI) calculated using DALECC-BUCKET")
ncnew_lai <- nc_create(f_out_lai, list(var_lai_all))

print("Writing data to file")
ncvar_put(ncnew,var_lai,lai_out)
ncvar_put(ncnew,var_gpp,gpp_out)
ncvar_put(ncnew,var_evap,evap_out)

# write lai from all model runs
ncvar_put(ncnew_lai,var_lai_all,lai)
nc_close(ncnew_lai)

print(paste("The file has", ncnew$nvars,"variable(s): lai, gpp, evap"))
print(paste("The file has", ncnew$ndim,"dimension(s): time, stats, lon, lat"))

# Don't forget to close the file
nc_close(ncnew)
print("File created!")






















