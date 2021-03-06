#!/usr/bin/Rscript

# can be DALEC_GSI_BUCKET or DALEC_GSI_DFOL_CWD_FR 
modelname <- "DALEC_GSI_DFOL_CWD_FR"

site <- "BALI GEM plots"
path2root <- "/home/dmilodow/DataStore_DTM/BALI/CARDAMOM_BALI/"
project <- "BALI_GEMplots_daily"
run <- "025"

path2rerun <- paste(path2root,"projects/",project,"/rerun/",run,"/", sep="")
path2mcmcfiles <- paste(path2root,"projects/",project,"/cardamom_output/",run,"/", sep="")
path2data <- paste(path2root,"projects/",project,"/data/",run,"/", sep="")
path2exe <- "/home/dmilodow/DataStore_DTM/BALI/CARDAMOM_BALI/CARDAMOM_DTM/src/CARDAMOM/trunk/LIBRARY/CARDAMOM_F/executable/"

integer_count_of_sites <- 6

startyear <- 2011
endyear <- 2017

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
	npools = 7
	npars = 39
	nfluxes = 24
}

# create some project related information
nopools=array(npools,integer_count_of_sites) ;
nopars=array(npars,integer_count_of_sites) ;
nofluxes=array(nfluxes,integer_count_of_sites)

ctessel_pft = 0   # distinction between crop model 1=crop, other number = generic (forest)

# in this case, all sites are crops
vector_of_site_pfts_1_is_crops = as.vector(rep(0, integer_count_of_sites))

# create vector of site ids, 00001, 00002, etc.
vector_of_site_names_or_ids <- 0
for(i in 1:integer_count_of_sites) {
    i_tmp = formatC(i, width = 5, format = "d", flag = "0")
    print(paste(i_tmp))
    vector_of_site_names_or_ids[i] = i_tmp
}

print("Removing old files")
for(i in 1:integer_count_of_sites) {
	if (file.exists(paste(path2rerun,vector_of_site_names_or_ids[i],".RData",sep=""))){
            file.remove(paste(path2rerun,vector_of_site_names_or_ids[i],".RData",sep=""))
        }
	if (file.exists(paste(path2rerun,vector_of_site_names_or_ids[i],"_parameters.RData",sep=""))){
            file.remove(paste(path2rerun,vector_of_site_names_or_ids[i],"_parameters.RData",sep=""))
        }
}

model = list(name=modelname,nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)

# nsubsamples is the target number of accepted parameters you wanted to end up with in the *PARS file, typically this is 1000.

PROJECT=list(name = project
                   ,exepath = path2exe
                   ,datapath = path2data
                   ,results_processedpath = path2rerun
                   ,resultspath = path2mcmcfiles 
                   ,nosites=integer_count_of_sites
                   ,sites = vector_of_site_names_or_ids
                   ,ctessel_pft = vector_of_site_pfts_1_is_crops
                   ,spatial_type = "site"
                   ,nsubsamples = 2000
                   ,latter_sample_frac = 0.5
		   ,parameter_type = "pft_specific"
                   ,model=model)

mcmc_results = run_mcmc_results(PROJECT)

# Plot up the parameter chains
#load(paste(path2rerun,"00006.RData",sep=""))
#par(mfrow=c(6,7))
#for (i in seq(1,40)) {plot(as.vector(parameters[i,,]),main=i)}
load(paste(path2rerun,"00001.RData",sep=""))
write.csv(file=paste(path2rerun,"00001.csv",sep=""),x=parameters[,501:1000,])
load(paste(path2rerun,"00002.RData",sep=""))
write.csv(file=paste(path2rerun,"00002.csv",sep=""),x=parameters[,501:1000,])
load(paste(path2rerun,"00003.RData",sep=""))
write.csv(file=paste(path2rerun,"00003.csv",sep=""),x=parameters[,501:1000,])
load(paste(path2rerun,"00004.RData",sep=""))
write.csv(file=paste(path2rerun,"00004.csv",sep=""),x=parameters[,501:1000,])
load(paste(path2rerun,"00005.RData",sep=""))
write.csv(file=paste(path2rerun,"00005.csv",sep=""),x=parameters[,501:1000,])
load(paste(path2rerun,"00006.RData",sep=""))
write.csv(file=paste(path2rerun,"00006.csv",sep=""),x=parameters[,501:1000,])
#############################################################################################











