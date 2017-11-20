#!/usr/bin/Rscript

# Script which reads in R binary file and writes the data to netcdf files

# can be DALEC_GSI_BUCKET or DALEC_GSI_DFOL_CWD_FR 
modelname <- "DALEC_GSI_DFOL_CWD_FR"

project <- "BALI_GEMplots_daily"
run <- "002"
site <- c("MLA01","MLA02","SAF04","SAF03","SAF02","SAF01")
lat <- c(4.747, 4.437, 4.765, 4.690, 4.744, 4.729)
lon <- c(116.951, 116.951, 117.702, 117.586, 117.618, 117.618)
startyear <- 2011
endyear <- 2016

path2root <- "/home/dmilodow/DataStore_DTM/BALI/CARDAMOM_BALI/"
path2files <- paste(path2root,"projects/",project,"/rerun/",run,"/", sep="")
nsites <- 6

#############################################################################################
# create vector of site ids, 00001, 00002, etc.
vector_of_site_ids <- 0
for(i in 1:nsites) {
	i_tmp = formatC(i, width = 5, format = "d", flag = "0")
	vector_of_site_ids[i] = i_tmp
}

# Setup output arrays
nstats = 5 # mean, median, standard deviation, 25th percentile and 75th percentile

bfile=paste(path2files,vector_of_site_ids[1],".RData",sep="")
print(paste("Reading in",bfile))
load(bfile)
lai=states_all$lai
ntsteps = length(lai[1,])

# compute statistics, for each timestep
lai_out=array(0., c(nsites,nstats,ntsteps)) 
gpp_out=array(0., c(nsites,nstats,ntsteps))
nee_out=array(0., c(nsites,nstats,ntsteps))
Reco_out=array(0., c(nsites,nstats,ntsteps))
Rauto_out=array(0., c(nsites,nstats,ntsteps))
Rhet_out=array(0., c(nsites,nstats,ntsteps))
Cwoo_out=array(0., c(nsites,nstats,ntsteps))
Csom_out=array(0., c(nsites,nstats,ntsteps))
Cbio_out=array(0., c(nsites,nstats,ntsteps))
Croo_out=array(0., c(nsites,nstats,ntsteps))
Clit_out=array(0., c(nsites,nstats,ntsteps))
Clab_out=array(0., c(nsites,nstats,ntsteps))
Cfol_out=array(0., c(nsites,nstats,ntsteps))
gsi_out=array(0., c(nsites,nstats,ntsteps))
gsi_itemp_out=array(0., c(nsites,nstats,ntsteps))
gsi_iphoto_out=array(0., c(nsites,nstats,ntsteps))
gsi_ivpd_out=array(0., c(nsites,nstats,ntsteps))
Ccwd_out=array(0., c(nsites,nstats,ntsteps))
flux_fol_lit_out=array(0., c(nsites,nstats,ntsteps))
flux_wood_cwd_out=array(0., c(nsites,nstats,ntsteps))
flux_root_lit_out=array(0., c(nsites,nstats,ntsteps))
flux_cwd_lit_out=array(0., c(nsites,nstats,ntsteps))
Rhet_lit_out=array(0., c(nsites,nstats,ntsteps))
decomp_lit_out=array(0., c(nsites,nstats,ntsteps))


for(ss in 1:nsites) {
    print(paste("SITE: ", site[ss]))
    # read in binary files
    bfile=paste(path2files,vector_of_site_ids[ss],".RData",sep="")

    print(paste("Reading in",bfile))
    load(bfile)

    paramsets=length(states_all$lai[,1])
    print(paste("Number of param sets is ",paramsets,sep=""))

    print("Calculating the mean, median, standard deviation, 25th percentile and 75th percentile values for each timestep")
    for(i in 1:ntsteps){
	lai_out[ss,1,i] <- mean(states_all$lai[,i]) 	
	lai_out[ss,2,i] <- median(states_all$lai[,i])
	lai_out[ss,3,i] <- sd(states_all$lai[,i]) # standard deviation
	lai_out[ss,4,i] <- quantile(states_all$lai[,i], 0.25) # 25th percentile
	lai_out[ss,5,i] <- quantile(states_all$lai[,i], 0.75) # 75th percentile

	gpp_out[ss,1,i] <- mean(states_all$gpp[,i]) 	
	gpp_out[ss,2,i] <- median(states_all$gpp[,i])
	gpp_out[ss,3,i] <- sd(states_all$gpp[,i])
	gpp_out[ss,4,i] <- quantile(states_all$gpp[,i], 0.25)
	gpp_out[ss,5,i] <- quantile(states_all$gpp[,i], 0.75)
        
	nee_out[ss,1,i] <- mean(states_all$nee[,i]) 	
	nee_out[ss,2,i] <- median(states_all$nee[,i])
	nee_out[ss,3,i] <- sd(states_all$nee[,i])
	nee_out[ss,4,i] <- quantile(states_all$nee[,i], 0.25)
	nee_out[ss,5,i] <- quantile(states_all$nee[,i], 0.75)

	Reco_out[ss,1,i] <- mean(states_all$reco[,i]) 	
	Reco_out[ss,2,i] <- median(states_all$reco[,i])
	Reco_out[ss,3,i] <- sd(states_all$reco[,i])
	Reco_out[ss,4,i] <- quantile(states_all$reco[,i], 0.25)
	Reco_out[ss,5,i] <- quantile(states_all$reco[,i], 0.75)
        
	Rauto_out[ss,1,i] <- mean(states_all$rauto[,i]) 	
	Rauto_out[ss,2,i] <- median(states_all$rauto[,i])
	Rauto_out[ss,3,i] <- sd(states_all$rauto[,i])
	Rauto_out[ss,4,i] <- quantile(states_all$rauto[,i], 0.25)
	Rauto_out[ss,5,i] <- quantile(states_all$rauto[,i], 0.75)
        
	Rhet_out[ss,1,i] <- mean(states_all$rhet[,i]) 	
	Rhet_out[ss,2,i] <- median(states_all$rhet[,i])
	Rhet_out[ss,3,i] <- sd(states_all$rhet[,i])
	Rhet_out[ss,4,i] <- quantile(states_all$rhet[,i], 0.25)
	Rhet_out[ss,5,i] <- quantile(states_all$rhet[,i], 0.75)
        
	Cwoo_out[ss,1,i] <- mean(states_all$wood[,i]) 	
	Cwoo_out[ss,2,i] <- median(states_all$wood[,i])
	Cwoo_out[ss,3,i] <- sd(states_all$wood[,i]) 
	Cwoo_out[ss,4,i] <- quantile(states_all$wood[,i], 0.25) 
	Cwoo_out[ss,5,i] <- quantile(states_all$wood[,i], 0.75) 

	Csom_out[ss,1,i] <- mean(states_all$som[,i]) 	
	Csom_out[ss,2,i] <- median(states_all$som[,i])
	Csom_out[ss,3,i] <- sd(states_all$som[,i])
	Csom_out[ss,4,i] <- quantile(states_all$som[,i], 0.25)
	Csom_out[ss,5,i] <- quantile(states_all$som[,i], 0.75)

	Cbio_out[ss,1,i] <- mean(states_all$bio[,i]) 	
	Cbio_out[ss,2,i] <- median(states_all$bio[,i])
	Cbio_out[ss,3,i] <- sd(states_all$bio[,i])
	Cbio_out[ss,4,i] <- quantile(states_all$bio[,i], 0.25)
	Cbio_out[ss,5,i] <- quantile(states_all$bio[,i], 0.75)

	Croo_out[ss,1,i] <- mean(states_all$root[,i]) 	
	Croo_out[ss,2,i] <- median(states_all$root[,i])
	Croo_out[ss,3,i] <- sd(states_all$root[,i])
	Croo_out[ss,4,i] <- quantile(states_all$root[,i], 0.25)
	Croo_out[ss,5,i] <- quantile(states_all$root[,i], 0.75)

	Clit_out[ss,1,i] <- mean(states_all$lit[,i]) 	
	Clit_out[ss,2,i] <- median(states_all$lit[,i])
	Clit_out[ss,3,i] <- sd(states_all$lit[,i])
	Clit_out[ss,4,i] <- quantile(states_all$lit[,i], 0.25)
	Clit_out[ss,5,i] <- quantile(states_all$lit[,i], 0.75)
        
	Clab_out[ss,1,i] <- mean(states_all$lab[,i]) 	
	Clab_out[ss,2,i] <- median(states_all$lab[,i])
	Clab_out[ss,3,i] <- sd(states_all$lab[,i])
	Clab_out[ss,4,i] <- quantile(states_all$lab[,i], 0.25)
	Clab_out[ss,5,i] <- quantile(states_all$lab[,i], 0.75)

	Cfol_out[ss,1,i] <- mean(states_all$fol[,i]) 	
	Cfol_out[ss,2,i] <- median(states_all$fol[,i])
	Cfol_out[ss,3,i] <- sd(states_all$fol[,i])
	Cfol_out[ss,4,i] <- quantile(states_all$fol[,i], 0.25)
	Cfol_out[ss,5,i] <- quantile(states_all$fol[,i], 0.75)
        
	gsi_out[ss,1,i] <- mean(states_all$gsi[,i]) 	
	gsi_out[ss,2,i] <- median(states_all$gsi[,i])
	gsi_out[ss,3,i] <- sd(states_all$gsi[,i])
	gsi_out[ss,4,i] <- quantile(states_all$gsi[,i], 0.25)
	gsi_out[ss,5,i] <- quantile(states_all$gsi[,i], 0.75)

	gsi_itemp_out[ss,1,i] <- mean(states_all$gsi_itemp[,i]) 	
	gsi_itemp_out[ss,2,i] <- median(states_all$gsi_itemp[,i])
	gsi_itemp_out[ss,3,i] <- sd(states_all$gsi_itemp[,i])
	gsi_itemp_out[ss,4,i] <- quantile(states_all$gsi_itemp[,i], 0.25)
	gsi_itemp_out[ss,5,i] <- quantile(states_all$gsi_itemp[,i], 0.75)
        
	gsi_iphoto_out[ss,1,i] <- mean(states_all$gsi_iphoto[,i]) 	
	gsi_iphoto_out[ss,2,i] <- median(states_all$gsi_iphoto[,i])
	gsi_iphoto_out[ss,3,i] <- sd(states_all$gsi_iphoto[,i])
	gsi_iphoto_out[ss,4,i] <- quantile(states_all$gsi_iphoto[,i], 0.25)
	gsi_iphoto_out[ss,5,i] <- quantile(states_all$gsi_iphoto[,i], 0.75)
        
	gsi_ivpd_out[ss,1,i] <- mean(states_all$gsi_ivpd[,i]) 	
	gsi_ivpd_out[ss,2,i] <- median(states_all$gsi_ivpd[,i])
	gsi_ivpd_out[ss,3,i] <- sd(states_all$gsi_ivpd[,i])
	gsi_ivpd_out[ss,4,i] <- quantile(states_all$gsi_ivpd[,i], 0.25)
	gsi_ivpd_out[ss,5,i] <- quantile(states_all$gsi_ivpd[,i], 0.75)
        
	Ccwd_out[ss,1,i] <- mean(states_all$litwood[,i]) 	
	Ccwd_out[ss,2,i] <- median(states_all$litwood[,i])
	Ccwd_out[ss,3,i] <- sd(states_all$litwood[,i])
	Ccwd_out[ss,4,i] <- quantile(states_all$litwood[,i], 0.25)
	Ccwd_out[ss,5,i] <- quantile(states_all$litwood[,i], 0.75)

	flux_fol_lit_out[ss,1,i] <- mean(states_all$litfol_flx[,i]) 	
	flux_fol_lit_out[ss,2,i] <- median(states_all$litfol_flx[,i])
	flux_fol_lit_out[ss,3,i] <- sd(states_all$litfol_flx[,i])
	flux_fol_lit_out[ss,4,i] <- quantile(states_all$litfol_flx[,i], 0.25)
	flux_fol_lit_out[ss,5,i] <- quantile(states_all$litfol_flx[,i], 0.75)
        
	flux_wood_cwd_out[ss,1,i] <- mean(states_all$litwood_flx[,i]) 	
	flux_wood_cwd_out[ss,2,i] <- median(states_all$litwood_flx[,i])
	flux_wood_cwd_out[ss,3,i] <- sd(states_all$litwood_flx[,i])
	flux_wood_cwd_out[ss,4,i] <- quantile(states_all$litwood_flx[,i], 0.25)
	flux_wood_cwd_out[ss,5,i] <- quantile(states_all$litwood_flx[,i], 0.75)

	flux_root_lit_out[ss,1,i] <- mean(states_all$litroot_flx[,i]) 	
	flux_root_lit_out[ss,2,i] <- median(states_all$litroot_flx[,i])
	flux_root_lit_out[ss,3,i] <- sd(states_all$litroot_flx[,i])
	flux_root_lit_out[ss,4,i] <- quantile(states_all$litroot_flx[,i], 0.25)
	flux_root_lit_out[ss,5,i] <- quantile(states_all$litroot_flx[,i], 0.75)
        
	flux_cwd_lit_out[ss,1,i] <- mean(states_all$litcwd_flux[,i]) 	
	flux_cwd_lit_out[ss,2,i] <- median(states_all$litcwd_flux[,i])
	flux_cwd_lit_out[ss,3,i] <- sd(states_all$litcwd_flux[,i])
	flux_cwd_lit_out[ss,4,i] <- quantile(states_all$litcwd_flux[,i], 0.25)
	flux_cwd_lit_out[ss,5,i] <- quantile(states_all$litcwd_flux[,i], 0.75)

	Rhet_lit_out[ss,1,i] <- mean(states_all$Rhet_lit[,i]) 	
	Rhet_lit_out[ss,2,i] <- median(states_all$Rhet_lit[,i])
	Rhet_lit_out[ss,3,i] <- sd(states_all$Rhet_lit[,i])
	Rhet_lit_out[ss,4,i] <- quantile(states_all$Rhet_lit[,i], 0.25)
	Rhet_lit_out[ss,5,i] <- quantile(states_all$Rhet_lit[,i], 0.75)

	decomp_lit_out[ss,1,i] <- mean(states_all$decomp_lit[,i]) 	
	decomp_lit_out[ss,2,i] <- median(states_all$decomp_lit[,i])
	decomp_lit_out[ss,3,i] <- sd(states_all$decomp_lit[,i])
	decomp_lit_out[ss,4,i] <- quantile(states_all$decomp_lit[,i], 0.25)
	decomp_lit_out[ss,5,i] <- quantile(states_all$decomp_lit[,i], 0.75)
    }
}
#############################################################################################
# create and write to netcdf
library(ncdf4)

f_out=paste(path2files,project, "_", startyear, "_", endyear, "_",modelname,".nc", sep="")

# define dimensions
sitedim <- ncdim_def("sites","dimensionless",1:nsites) 
timedim <- ncdim_def("time","weeks",1:ntsteps,unlim=TRUE)
stats <- ncdim_def("stats","dimenisonless (1-5)",1:nstats) 

var_lai <- ncvar_def("lai", "m2m-2", list(sitedim,stats,timedim), longname="Leaf Area Index (LAI)")
var_gpp <- ncvar_def("gpp", "gC m-2day-1", list(sitedim,stats,timedim), longname="Gross Primary Productivity (GPP)")
var_nee <- ncvar_def("nee", "gC m-2day-1", list(sitedim,stats,timedim), longname="Net Ecosystem Exchange (NEE) ")
var_Reco <- ncvar_def("Reco", "gC m-2day-1", list(sitedim,stats,timedim), longname="Ecosystem respiration (Reco)")
var_Rauto <- ncvar_def("Rauto", "gC m-2day-1", list(sitedim,stats,timedim), longname="Autotrophic respiration (Rauto) ")
var_Rhet <- ncvar_def("Rhet", "gC m-2day-1", list(sitedim,stats,timedim), longname="Heterotrophic respiration (Rhet) ")
var_Cwoo <- ncvar_def("Cwoo", "gC m-2", list(sitedim,stats,timedim), longname="Wood carbon stock (Cwoo)")
var_Clit <- ncvar_def("Clit", "gC m-2", list(sitedim,stats,timedim), longname="Litter carbon stock (Clit)")
var_Clab <- ncvar_def("Clab", "gC m-2", list(sitedim,stats,timedim), longname="Labile carbon stock (Clab)")
var_Csom <- ncvar_def("Csom", "gC m-2", list(sitedim,stats,timedim), longname="Soil organic matter carbon stock (Csom)")
var_Croo <- ncvar_def("Croo", "gC m-2", list(sitedim,stats,timedim), longname="Root carbon stock (Croo)")
var_Cbio <- ncvar_def("Cbio", "gC m-2", list(sitedim,stats,timedim), longname="Aggregated biotic carbon stock (Cbio)")
var_Cfol <- ncvar_def("Cfol", "gC m-2", list(sitedim,stats,timedim), longname="Foliage carbon stock (Cfol)")
var_gsi <- ncvar_def("gsi", "unspecified units", list(sitedim,stats,timedim), longname="Growing Season Index (GSI)")
var_gsi_itemp <- ncvar_def("gsi_itemp", "unspecified units", list(sitedim,stats,timedim), longname="Growing Season Index (GSI) - temperature component")
var_gsi_iphoto <- ncvar_def("gsi_iphoto", "unspecified units", list(sitedim,stats,timedim), longname="Growing Season Index (GSI) - photoperiod component")
var_gsi_ivpd <- ncvar_def("gsi_ivpd", "unspecified units", list(sitedim,stats,timedim), longname="Growing Season Index (GSI) - vapour pressure deficit component")
var_Ccwd <- ncvar_def("Ccwd", "gC m-2", list(sitedim,stats,timedim), longname="Coarse woody debris carbon stock (Ccwd)")
var_flux_fol_lit <- ncvar_def("flux_fol_lit", "gC m-2day-1", list(sitedim,stats,timedim), longname="Carbon flux from foliage to litter")
var_flux_wood_cwd <- ncvar_def("flux_wood_cwd", "gC m-2day-1", list(sitedim,stats,timedim), longname="Carbon flux from wood to cwd")
var_flux_root_lit <- ncvar_def("flux_root_lit", "gC m-2day-1", list(sitedim,stats,timedim), longname="Carbon flux from roots to litter")
var_flux_cwd_lit <- ncvar_def("flux_cwd_lit", "gC m-2day-1", list(sitedim,stats,timedim), longname="Carbon flux from cwd to litter")
var_Rh_lit <- ncvar_def("Rh_lit", "gC m-2day-1", list(sitedim,stats,timedim), longname="Heterotrophic respiration flux from litter")
var_decomp_lit <- ncvar_def("decomp_lit", "gC m-2day-1", list(sitedim,stats,timedim), longname="decomposition flux from litter")
var_lat <- ncvar_def("latitude", "degrees_north", list(sitedim), longname="latitude of site")
var_lon <- ncvar_def("longitude", "decimal_east", list(sitedim), longname="longitude of site")
var_site <- ncvar_def("GEM_code", "text", list(sitedim), longname="GEM code for plot")

ncnew <- nc_create(f_out, list(var_lai,var_gpp,var_nee,var_Reco,var_Rauto,var_Rhet,var_Cwoo,var_Clab,var_Cfol,var_Croo,var_Clit,var_Ccwd,var_Csom,var_Cbio,var_gsi,var_gsi_itemp,var_gsi_iphoto,var_gsi_ivpd,var_flux_fol_lit,var_flux_wood_cwd,var_flux_root_lit,var_flux_cwd_lit,var_Rh_lit,var_decomp_lit,var_lat,var_lon,var_site))

print("Writing data to file")
ncvar_put(ncnew,var_lai,lai_out)
ncvar_put(ncnew,var_gpp,gpp_out)
ncvar_put(ncnew,var_nee,nee_out)
ncvar_put(ncnew,var_Reco,Reco_out)
ncvar_put(ncnew,var_Rauto,Rauto_out)
ncvar_put(ncnew,var_Rhet,Rhet_out)
ncvar_put(ncnew,var_Cwoo,Cwoo_out)
ncvar_put(ncnew,var_Clab,Clab_out)
ncvar_put(ncnew,var_Cfol,Cfol_out)
ncvar_put(ncnew,var_Croo,Croo_out)
ncvar_put(ncnew,var_Clit,Clit_out)
ncvar_put(ncnew,var_Ccwd,Ccwd_out)
ncvar_put(ncnew,var_Csom,Csom_out)
ncvar_put(ncnew,var_Cbio,Cbio_out)
ncvar_put(ncnew,var_gsi,gsi_out)
ncvar_put(ncnew,var_gsi_itemp,gsi_itemp_out)
ncvar_put(ncnew,var_gsi_iphoto,gsi_iphoto_out)
ncvar_put(ncnew,var_gsi_ivpd,gsi_ivpd_out)
ncvar_put(ncnew,var_flux_fol_lit,flux_fol_lit_out)
ncvar_put(ncnew,var_flux_wood_cwd,flux_wood_cwd_out)
ncvar_put(ncnew,var_flux_root_lit,flux_root_lit_out)
ncvar_put(ncnew,var_flux_cwd_lit,flux_cwd_lit_out)
ncvar_put(ncnew,var_Rh_lit,Rhet_lit_out)
ncvar_put(ncnew,var_decomp_lit,decomp_lit_out)
ncvar_put(ncnew,var_lat,lat)
ncvar_put(ncnew,var_lon,lon)
ncvar_put(ncnew,var_site,site)

print(paste("The file has", ncnew$ndim,"dimension(s): sites, time, stats"))

# Don't forget to close the file
nc_close(ncnew)
print("File created!")






















