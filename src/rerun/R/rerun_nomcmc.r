#!/usr/bin/Rscript

# can be DALEC_GSI_BUCKET or DALEC_GSI_DFOL_CWD_FR 
modelname <- "DALEC_GSI_DFOL_CWD_FR"

site <- "SAFE"
project <- "BALI_GEMplots_daily"
run <- "001"

mcmcoutputdirectory <- paste("/home/dmilodow/DataStore_DTM/BALI/CARDAMOM_BALI/projects/",project,"/rerun/",run)

path2files <- paste(site, "/PARS/", sep="")
path2files <- "/home/dmilodow/DataStore_DTM/BALI/CARDAMOM_BALI/projects/BALI_GEMplots_daily/cardamom_output/001"

integer_count_of_sites <- 6

startyear <- 2011
endyear <- 2016

#############################################################################################
# load needed libraries
library(parallel)
library(compiler)

# set working directory
setwd("/exports/csce/datastore/geos/users/dmilodow/BALI/CARDAMOM_BALI/CARDAMOM_DTM/src/rerun/R")

# read in R functions
source("psrf_function.r")
source("read_binary_file_format.r")
source("read_parameter_chains.r")
source("simulate_all.r")
source("run_mcmc_results.r")
source("have_chains_converged.r")
source("acm_tessel_crop_r.r")

print(paste("Running MCMC results (",modelname,") at ",site,sep=""))

# set global variables needed
repair = 0 # overwrite exisiting file outputs
grid_override = TRUE # if running a gridded analysis save stock and flux information (TRUE)
use_parallel = FALSE # use parallel functions
numWorkers = 1 # if use parallel == TRUE, how many cores to use?
specific_pft = "pft_specific"

if (modelname == "DALEC_GSI_BUCKET"){
	npools = 9
	npars = 37
	nfluxes = 21	
} else if (modelname == "DALEC_GSI_DFOL_CWD_FR"){
	npools = 8
	npars = 35
	nfluxes = 16	
}

# create some project related information
nopools=array(npools,integer_count_of_sites) ;
nopars=array(npars,integer_count_of_sites) ;
nofluxes=array(nfluxes,integer_count_of_sites)

ctessel_pft = 1

# in this case, all sites are crops
vector_of_site_pfts_1_is_crops = as.vector(rep(1, integer_count_of_sites))

#print(vector_of_site_pfts_1_is_crops)

# create vector of site ids, 00001, 00002, etc.
vector_of_site_names_or_ids <- 0
for(i in 1:integer_count_of_sites) {
	i_tmp = formatC(i, width = 5, format = "d", flag = "0")
	vector_of_site_names_or_ids[i] = i_tmp
}

print("Removing old files")
if (modelname == "DALEC_GSI_BUCKET"){
	f_out = paste(site, "_weekly_Crop_", startyear, "_", endyear, "_bucket", sep="")		
} else if (modelname == "DALEC_GSI_DFOL_CWD_FR"){
	f_out = paste(site, "_weekly_Crop_", startyear, "_", endyear, "_no_bucket", sep="")
}

for(i in 1:integer_count_of_sites) {
	if (file.exists(paste(site,"/DALECc_output/",f_out,"_",vector_of_site_names_or_ids[i],".RData",sep=""))) file.remove(paste(site,"/DALECc_output/",f_out,"_",vector_of_site_names_or_ids[i],".RData",sep=""))
	if (file.exists(paste(site,"/DALECc_output/",f_out,"_",vector_of_site_names_or_ids[i],"_parameters.RData",sep=""))) file.remove(paste(site,"/DALECc_output/",f_out,"_",vector_of_site_names_or_ids[i],"_parameters.RData",sep=""))	
}

#if (specific_pft == "pft_specific") {nopars[which(ctessel_pft ==
#1)]=35+2 ; nofluxes[which(ctessel_pft == 1)]=21 ;
#nopools[which(ctessel_pft == 1)]=9}

model = list(name=modelname,nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)

# nsubsamples is the target number of accepted parameters you wanted to end up with in the *PARS file, typically this is 1000.

PROJECT=list(name = f_out
                   ,exepath = "source/"
                   ,datapath = paste(site, "/bin/", sep="")
                   ,results_processedpath = paste(site, "/DALECc_output/", sep="")
                   ,resultspath = paste(site, "/PARS", sep="")
                   ,nosites=integer_count_of_sites
                   ,sites = vector_of_site_names_or_ids
                   ,ctessel_pft = vector_of_site_pfts_1_is_crops
                   ,spatial_type = "site"
                   ,nsubsamples = 1000
                   ,latter_sample_frac = 0.5
		   ,parameter_type = "pft_specific"
                   ,model=model)

mcmc_results = run_mcmc_results(PROJECT)

#############################################################################################















