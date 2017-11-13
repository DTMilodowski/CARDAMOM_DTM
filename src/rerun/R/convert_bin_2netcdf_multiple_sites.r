#!/usr/bin/Rscript

# Script which reads in R binary file and writes the data to netcdf files

# can be DALEC_GSI_BUCKET or DALEC_GSI_DFOL_CWD_FR 
modelname <- "DALEC_GSI_DFOL_CWD_FR"

project <- "BALI_GEMplots_daily"
run <- "001"
site <- c("MLA01","MLA02","SAF03","SAF04","SAF02","SAF01")
lat <- c(4.747, 4.437, 4.690, 4,765, 4,744, 4,729)
lon <- c(116.951, 116.951, 117.586, 117.702, 117.618, 117.618)
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
    # read in binary files
    bfile=paste(path2files,vector_of_site_ids[ss],".RData",sep="")

    print(paste("Reading in",bfile))
    load(bfile)

    lai=states_all$lai
    gpp=states_all$gpp
    nee=states_all$nee
    reco=states_all$reco
    rauto=states_all$rauto
    rhet=states_all$rhet
    wood=states_all$wood
    som=states_all$som
    bio=states_all$bio
    root=states_all$root
    lab=states_all$lab
    lit=states_all$lit
    fol=states_all$fol
    gsi=states_all$gsi
    gsi_itemp=states_all$gsi_itemp
    gsi_iphoto=states_all$gsi_iphoto
    litwood=states_all$litwood
    litfol_flx=states_all$litfol_flx
    litwood_flx=states_all$litwood_flx
    litroot_flx=states_all$litroot_flx
    litcwd_flx=states_all$litcwd_flux
    decomp_lit=states_all$decomp_lit
    Rhet_lit=states_all$Rhet_lit


    ntsteps=length(lai[1,])
    paramsets=length(lai[,1])
    print(paste("Number of param sets is ",paramsets,sep=""))

    print("Calculating the mean, median, standard deviation, 25th percentile and 75th percentile values for each timestep")
    for(i in 1:ntsteps){
	lai_out[ss,1,i] <- mean(lai[,i]) 	
	lai_out[ss,2,i] <- median(lai[,i])
	lai_out[ss,3,i] <- sd(lai[,i]) # standard deviation
	lai_out[ss,4,i] <- quantile(lai[,i], 0.25) # 25th percentile
	lai_out[ss,5,i] <- quantile(lai[,i], 0.75) # 75th percentile

	gpp_out[ss,1,i] <- mean(gpp[,i]) 	
	gpp_out[ss,2,i] <- median(gpp[,i])
	gpp_out[ss,3,i] <- sd(gpp[,i])
	gpp_out[ss,4,i] <- quantile(gpp[,i], 0.25)
	gpp_out[ss,5,i] <- quantile(gpp[,i], 0.75)
        
	nee_out[ss,1,i] <- mean(nee[,i]) 	
	nee_out[ss,2,i] <- median(nee[,i])
	nee_out[ss,3,i] <- sd(nee[,i])
	nee_out[ss,4,i] <- quantile(nee[,i], 0.25)
	nee_out[ss,5,i] <- quantile(nee[,i], 0.75)

	Reco_out[ss,1,i] <- mean(reco[,i]) 	
	Reco_out[ss,2,i] <- median(reco[,i])
	Reco_out[ss,3,i] <- sd(reco[,i])
	Reco_out[ss,4,i] <- quantile(reco[,i], 0.25)
	Reco_out[ss,5,i] <- quantile(reco[,i], 0.75)
        
	Rauto_out[ss,1,i] <- mean(rauto[,i]) 	
	Rauto_out[ss,2,i] <- median(rauto[,i])
	Rauto_out[ss,3,i] <- sd(rauto[,i])
	Rauto_out[ss,4,i] <- quantile(rauto[,i], 0.25)
	Rauto_out[ss,5,i] <- quantile(rauto[,i], 0.75)
        
	Rhet_out[ss,1,i] <- mean(rhet[,i]) 	
	Rhet_out[ss,2,i] <- median(rhet[,i])
	Rhet_out[ss,3,i] <- sd(rhet[,i])
	Rhet_out[ss,4,i] <- quantile(rhet[,i], 0.25)
	Rhet_out[ss,5,i] <- quantile(rhet[,i], 0.75)
        
	Cwoo_out[ss,1,i] <- mean(wood[,i]) 	
	Cwoo_out[ss,2,i] <- median(wood[,i])
	Cwoo_out[ss,3,i] <- sd(wood[,i]) 
	Cwoo_out[ss,4,i] <- quantile(wood[,i], 0.25) 
	Cwoo_out[ss,5,i] <- quantile(wood[,i], 0.75) 

	Csom_out[ss,1,i] <- mean(som[,i]) 	
	Csom_out[ss,2,i] <- median(som[,i])
	Csom_out[ss,3,i] <- sd(som[,i])
	Csom_out[ss,4,i] <- quantile(som[,i], 0.25)
	Csom_out[ss,5,i] <- quantile(som[,i], 0.75)

	Cbio_out[ss,1,i] <- mean(bio[,i]) 	
	Cbio_out[ss,2,i] <- median(bio[,i])
	Cbio_out[ss,3,i] <- sd(bio[,i])
	Cbio_out[ss,4,i] <- quantile(bio[,i], 0.25)
	Cbio_out[ss,5,i] <- quantile(bio[,i], 0.75)

	Croo_out[ss,1,i] <- mean(root[,i]) 	
	Croo_out[ss,2,i] <- median(root[,i])
	Croo_out[ss,3,i] <- sd(root[,i])
	Croo_out[ss,4,i] <- quantile(root[,i], 0.25)
	Croo_out[ss,5,i] <- quantile(root[,i], 0.75)

	Clit_out[ss,1,i] <- mean(lit[,i]) 	
	Clit_out[ss,2,i] <- median(lit[,i])
	Clit_out[ss,3,i] <- sd(lit[,i])
	Clit_out[ss,4,i] <- quantile(lit[,i], 0.25)
	Clit_out[ss,5,i] <- quantile(lit[,i], 0.75)
        
	Clab_out[ss,1,i] <- mean(lab[,i]) 	
	Clab_out[ss,2,i] <- median(lab[,i])
	Clab_out[ss,3,i] <- sd(lab[,i])
	Clab_out[ss,4,i] <- quantile(lab[,i], 0.25)
	Clab_out[ss,5,i] <- quantile(lab[,i], 0.75)

	Cfol_out[ss,1,i] <- mean(fol[,i]) 	
	Cfol_out[ss,2,i] <- median(fol[,i])
	Cfol_out[ss,3,i] <- sd(fol[,i])
	Cfol_out[ss,4,i] <- quantile(fol[,i], 0.25)
	Cfol_out[ss,5,i] <- quantile(fol[,i], 0.75)
        
	gsi_out[ss,1,i] <- mean(gsi[,i]) 	
	gsi_out[ss,2,i] <- median(gsi[,i])
	gsi_out[ss,3,i] <- sd(gsi[,i])
	gsi_out[ss,4,i] <- quantile(gsi[,i], 0.25)
	gsi_out[ss,5,i] <- quantile(gsi[,i], 0.75)

	gsi_itemp_out[ss,1,i] <- mean(gsi_itemp[,i]) 	
	gsi_itemp_out[ss,2,i] <- median(gsi_itemp[,i])
	gsi_itemp_out[ss,3,i] <- sd(gsi_itemp[,i])
	gsi_itemp_out[ss,4,i] <- quantile(gsi_itemp[,i], 0.25)
	gsi_itemp_out[ss,5,i] <- quantile(gsi_itemp[,i], 0.75)
        
	gsi_iphoto_out[ss,1,i] <- mean(gsi_iphoto[,i]) 	
	gsi_iphoto_out[ss,2,i] <- median(gsi_iphoto[,i])
	gsi_iphoto_out[ss,3,i] <- sd(gsi_iphoto[,i])
	gsi_iphoto_out[ss,4,i] <- quantile(gsi_iphoto[,i], 0.25)
	gsi_iphoto_out[ss,5,i] <- quantile(gsi_iphoto[,i], 0.75)
        
	Ccwd_out[ss,1,i] <- mean(litwood[,i]) 	
	Ccwd_out[ss,2,i] <- median(litwood[,i])
	Ccwd_out[ss,3,i] <- sd(litwood[,i])
	Ccwd_out[ss,4,i] <- quantile(litwood[,i], 0.25)
	Ccwd_out[ss,5,i] <- quantile(litwood[,i], 0.75)

	flux_fol_lit_out[ss,1,i] <- mean(litfol_flx[,i]) 	
	flux_fol_lit_out[ss,2,i] <- median(litfol_flx[,i])
	flux_fol_lit_out[ss,3,i] <- sd(litfol_flx[,i])
	flux_fol_lit_out[ss,4,i] <- quantile(litfol_flx[,i], 0.25)
	flux_fol_lit_out[ss,5,i] <- quantile(litfol_flx[,i], 0.75)
        
	flux_wood_cwd_out[ss,1,i] <- mean(litwood_flx[,i]) 	
	flux_wood_cwd_out[ss,2,i] <- median(litwood_flx[,i])
	flux_wood_cwd_out[ss,3,i] <- sd(litwood_flx[,i])
	flux_wood_cwd_out[ss,4,i] <- quantile(litwood_flx[,i], 0.25)
	flux_wood_cwd_out[ss,5,i] <- quantile(litwood_flx[,i], 0.75)

	flux_root_lit_out[ss,1,i] <- mean(litroot_flx[,i]) 	
	flux_root_lit_out[ss,2,i] <- median(litroot_flx[,i])
	flux_root_lit_out[ss,3,i] <- sd(litroot_flx[,i])
	flux_root_lit_out[ss,4,i] <- quantile(litroot_flx[,i], 0.25)
	flux_root_lit_out[ss,5,i] <- quantile(litroot_flx[,i], 0.75)
        
	flux_cwd_lit_out[ss,1,i] <- mean(litcwd_flx[,i]) 	
	flux_cwd_lit_out[ss,2,i] <- median(litcwd_flx[,i])
	flux_cwd_lit_out[ss,3,i] <- sd(litcwd_flx[,i])
	flux_cwd_lit_out[ss,4,i] <- quantile(litcwd_flx[,i], 0.25)
	flux_cwd_lit_out[ss,5,i] <- quantile(litcwd_flx[,i], 0.75)

	Rhet_lit_out[ss,1,i] <- mean(Rhet_lit[,i]) 	
	Rhet_lit_out[ss,2,i] <- median(Rhet_lit[,i])
	Rhet_lit_out[ss,3,i] <- sd(Rhet_lit[,i])
	Rhet_lit_out[ss,4,i] <- quantile(Rhet_lit[,i], 0.25)
	Rhet_lit_out[ss,5,i] <- quantile(Rhet_lit[,i], 0.75)

	decomp_lit_out[ss,1,i] <- mean(decomp_lit[,i]) 	
	decomp_lit_out[ss,2,i] <- median(decomp_lit[,i])
	decomp_lit_out[ss,3,i] <- sd(decomp_lit[,i])
	decomp_lit_out[ss,4,i] <- quantile(decomp_lit[,i], 0.25)
	decomp_lit_out[ss,5,i] <- quantile(decomp_lit[,i], 0.75)
    }
}
#############################################################################################
# create and write to netcdf
library(ncdf4)

f_out=paste(path2files,project, "_", startyear, "_", endyear, "_",modelname,".nc", sep="")

## write out all runs for lai
#f_out_lai=paste(path2files,site, "_weekly_Crop_", startyear, "_", endyear, "_bucket_lai",".nc", sep="")

xvals <- lon
yvals <- lat 
#nx <- length(xvals)
#ny <- length(yvals)

# define dimensions
#londim <- ncdim_def("lon","degrees_east",xvals) 
#latdim <- ncdim_def("lat","degrees_north",yvals) 
sitedim <- ncdim_def("sites","dimensionless",1:nsites) 
timedim <- ncdim_def("time","weeks",1:ntsteps,unlim=TRUE)
stats <- ncdim_def("stats","dimenisonless (1-5)",1:nstats) 

var_lai <- ncvar_def("lai", "m2m-2", list(stats,timedim,sitedim), longname="Leaf Area Index (LAI)")
var_gpp <- ncvar_def("gpp", "gC m-2day-1", list(stats,timedim,sitedim), longname="Gross Primary Productivity (GPP)")
var_nee <- ncvar_def("nee", "gC m-2day-1", list(stats,timedim,sitedim), longname="Net Ecosystem Exchange (NEE) ")
var_Reco <- ncvar_def("Reco", "gC m-2day-1", list(stats,timedim,sitedim), longname="Ecosystem respiration (Reco)")
var_Rauto <- ncvar_def("Rauto", "gC m-2day-1", list(stats,timedim,sitedim), longname="Autotrophic respiration (Rauto) ")
var_Rhet <- ncvar_def("Rhet", "gC m-2day-1", list(stats,timedim,sitedim), longname="Heterotrophic respiration (Rhet) ")
var_Cwoo <- ncvar_def("Cwoo", "gC m-2", list(stats,timedim,sitedim), longname="Wood carbon stock (Cwoo)")
var_Clit <- ncvar_def("Clit", "gC m-2", list(stats,timedim,sitedim), longname="Litter carbon stock (Clit)")
var_Clab <- ncvar_def("Clab", "gC m-2", list(stats,timedim,sitedim), longname="Labile carbon stock (Clab)")
var_Csom <- ncvar_def("Csom", "gC m-2", list(stats,timedim,sitedim), longname="Soil organic matter carbon stock (Csom)")
var_Croo <- ncvar_def("Croo", "gC m-2", list(stats,timedim,sitedim), longname="Root carbon stock (Croo)")
var_Cbio <- ncvar_def("Cbio", "gC m-2", list(stats,timedim,sitedim), longname="Aggregated biotic carbon stock (Cbio)")
var_Cfol <- ncvar_def("Cfol", "gC m-2", list(stats,timedim,sitedim), longname="Foliage carbon stock (Cfol)")
var_gsi <- ncvar_def("gsi", "unspecified units", list(stats,timedim,sitedim), longname="Growing Season Index (GSI)")
var_gsi_itemp <- ncvar_def("gsi_itemp", "unspecified units", list(stats,timedim,sitedim), longname="Growing Season Index (GSI) - temperature component")
var_gsi_iphoto <- ncvar_def("gsi_iphoto", "unspecified units", list(stats,timedim,sitedim), longname="Growing Season Index (GSI) - photoperiod component")
var_gsi_ivpd <- ncvar_def("gsi_ivpd", "unspecified units", list(stats,timedim,sitedim), longname="Growing Season Index (GSI) - vapour pressure deficit component")
var_Ccwd <- ncvar_def("Ccwd", "gC m-2", list(stats,timedim,sitedim), longname="Coarse woody debris carbon stock (Ccwd)")
var_flux_fol_lit <- ncvar_def("flux_fol_lit", "gC m-2day-1", list(stats,timedim,sitedim), longname="Carbon flux from foliage to litter")
var_flux_wood_cwd <- ncvar_def("flux_wood_cwd", "gC m-2day-1", list(,stats,timedim,sitedim), longname="Carbon flux from wood to cwd")
var_flux_root_lit <- ncvar_def("flux_root_lit", "gC m-2day-1", list(stats,timedim,sitedim), longname="Carbon flux from roots to litter")
var_flux_cwd_lit <- ncvar_def("flux_cwd_lit", "gC m-2day-1", list(stats,timedim,sitedim), longname="Carbon flux from cwd to litter")
var_Rh_lit <- ncvar_def("Rh_lit", "gC m-2day-1", list(stats,timedim,sitedim), longname="Heterotrophic respiration flux from litter")
var_decomp_lit <- ncvar_def("decomp_lit", "gC m-2day-1", list(stats,timedim,sitedim), longname="decomposition flux from litter")
var_lat <- ncvar_def("latitude", "degrees_north", list(sitedim), longname="latitude of site")
var_lon <- ncvar_def("longitude", "decimal_east", list(sitedim), longname="longitude of site")
var_site <- ncvar_def("GEM_code", "text", list(sitedim), longname="GEM code for plot")

ncnew <- nc_create(f_out, list(var_lai,var_gpp,var_nee,var_Reco,var_Rauto,var_Rhet,var_Cwoo,var_Clab,var_Cfol,var_Croo,var_Clit,var_Ccwd,var_Csom,var_Cbio,var_gsi,var_gsi_itemp,var_gsi_iphoto,var_gsi_ivpd,var_flux_fol_lit,var_flux_wood_cwd,var_flux_root_lit,var_flux_cwd_lit,var_Rh_lit,var_decomp_lit,var_lat,var_lon,var_site))

## write lai from all model dalecc runs to file
#parsdim <- ncdim_def("params","nparams",1:paramsets)
#var_lai_all <- ncvar_def("lai", "m2m-2", list(latdim,londim,parsdim,timedim), longname="Leaf Area Index (LAI) calculated using DALECC-BUCKET")
#ncnew_lai <- nc_create(f_out_lai, list(var_lai_all))

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

## write lai from all model runs
#ncvar_put(ncnew_lai,var_lai_all,lai)
#nc_close(ncnew_lai)

#print(paste("The file has", ncnew$nvars,"variable(s): lai, gpp"))
print(paste("The file has", ncnew$ndim,"dimension(s): sites, time, stats"))#lon, lat"))

# Don't forget to close the file
nc_close(ncnew)
print("File created!")






















