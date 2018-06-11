
###
## Function to run CARDAMOM parameters via the choosen model
###

run_each_site<-function(n,PROJECT,stage,repair,grid_override,stage5modifiers) {

	if (stage == 5) {
	    stage5name = as.vector(unlist(stage5modifiers))
	    stage5name = paste(stage5name[1],stage5name[2],stage5name[3],stage5name[4],stage5name[5],stage5name[6],stage5name[7],stage5name[8],sep="_")
	    outfile=paste(PROJECT$results_processedpath,PROJECT$sites[n],"_",stage5name,".RData",sep="")
	    outfile1=paste(PROJECT$results_processedpath,PROJECT$sites[n],"_",stage5name,"_parameters.RData",sep="")
	} else { 
	    outfile=paste(PROJECT$results_processedpath,PROJECT$sites[n],".RData",sep="")
	    outfile1=paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")
	}

	if (file.exists(outfile1) == FALSE | repair == 1) {
	    # load only the desired latter fraction of the parameter vectors
	    # output is order dimensions(npar+1,iter,chain)
	    parameters=read_parameter_chains(PROJECT,n,3)

	    # ok so if we actually have some parameters we will 
            if (parameters[1] != -9999) {
		# test for convergence and whether or not there is any single chain which can be removed in they do not converge
		if (dim(parameters)[3] > 2) {
		    converged=have_chains_converged(parameters)
		    # if log-likelihood has passed then we are not interested
		    if (converged[length(converged)] == "FAIL") {
			notconv = TRUE ; i=1 ; max_likelihood=rep(NA, length.out=dim(parameters)[3]) ; CI90=rep(NA,length.out=c(2))
			while (notconv){
			    max_likelihood[i] = max(parameters[dim(parameters)[1],,i])
			    converged=have_chains_converged(parameters[,,-i]) ; i = i+1
			    # if removing one of the chains get convergence then great
			    if (converged[length(converged)] == "PASS") {
				# but we need to check for the possibility that the chain we have removed is actually better than the others
				CI90[1] = quantile(parameters[dim(parameters)[1],,(i-1)], prob=c(0.10)) ; CI90[2] = quantile(parameters[dim(parameters)[1],,-(i-1)], prob=c(0.90))
				# if the rejected chain is significantly better (at 90 % CI) than the converged chains then we have a problem
				if (CI90[1] > CI90[2]) {
				    # we we do nothing and allow the analysis to continue moving on
				} else {
				    # if the non-converged chain is worse or just the same in likelihood terms as the others then we will ditch it
				    notconv=FALSE ; i=i-1 # converged now?
				}
			    }
			    # if we have tried removing each chain and we still have not converged then give up anyway
			    #if (i > dim(parameters)[3] & notconv) {notconv=FALSE ; i=-9999}
			    # or actually just remove the lowest average likelihood chain
			    if (i > dim(parameters)[3] & notconv) {notconv=FALSE ; i=which(max_likelihood == min(max_likelihood)) ; i=i[1]}
			} # for removing chains
			# if we successfully found only chain to remove then remove it from the rest of the analysis.
			if (i != -9999) {parameters=parameters[,,-i] ; print(paste("chain rejected = ",i,sep=""))}
		    } # if likelihood not converged
		} # if more than 2 chains

		# load the met data for each site
		drivers=read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))
                # if Stage 5 sensitivity analysis is enacted - apply factor scaled adjustment to the drivers
                if (stage == 5) {

                    # do weather drivers first
                    drivers$met[,2]  = drivers$met[,2]  * stage5modifiers$airt_factor # min temperature
                    drivers$met[,3]  = drivers$met[,3]  * stage5modifiers$airt_factor # max temperature
                    drivers$met[,4]  = drivers$met[,4]  * stage5modifiers$swrad_factor
                    drivers$met[,5]  = drivers$met[,5]  * stage5modifiers$co2_factor
                    drivers$met[,7]  = drivers$met[,7]  * stage5modifiers$rainfall_factor
                    drivers$met[,15] = drivers$met[,15] * stage5modifiers$wind_spd_factor
                    drivers$met[,16] = drivers$met[,16] * stage5modifiers$vpd_factor
                    # do disturbance drivers next
                    disturbed_locations = which(drivers$met[,8] > 0)
                    drivers$met[disturbed_locations,8] = drivers$met[disturbed_locations,8] * stage5modifiers$deforestation_factor
                    disturbed_locations = which(drivers$met[,9] > 0)
                    drivers$met[disturbed_locations,9] = drivers$met[disturbed_locations,9] * stage5modifiers$burnt_area_factor
                    # do GSI related parameters next
                    drivers$met[,10] = drivers$met[,10] * stage5modifiers$airt_factor # 21 day rolling average daily minimum temperature

                }
		# now run the parameters to generate state variables
		# only use 500 for full propogation of states and fluxes
		if (dim(parameters)[2] > 500) {
		    print("Note only 500 parameter sets PER CHAIN will be used in full propogation of parameter / flux / state uncertainties")
		    sub_parameter=parameters[,sample(dim(parameters)[2],500,replace=FALSE),]
		} else {
		    sub_parameter=parameters
		}
		# run subsample of parameters for full results / propogation
                soil_info=c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
		states_all=simulate_all(n,PROJECT,PROJECT$model$name,drivers$met,sub_parameter[1:PROJECT$model$nopars[n],,],drivers$lat,PROJECT$ctessel_pft[n],PROJECT$parameter_type,PROJECT$exepath,soil_info)
		# pass to local variable for saving
		site_ctessel_pft=PROJECT$ctessel_pft[n]
		aNPP=states_all$aNPP 
		# store the results now in binary file
		save(parameters,drivers,sub_parameter,site_ctessel_pft,aNPP,file=outfile1)#, compress="gzip")
		if (PROJECT$spatial_type == "site" | grid_override == TRUE) {
		    save(parameters,drivers,states_all,sub_parameter,site_ctessel_pft,aNPP,file=outfile)#, compress="gzip")
		}
		dummy=0
	    } else {
		dummy=-1
	    } # parameters[1] != -9999
	} else {
	    dummy=-1
	    print('Already extracted result vectors (set repair = 1 if re-run is needed)')
	} # *parameters.RData already exists

	return(dummy)

} # end of run_each_site
## Use byte compile
run_each_site<-cmpfun(run_each_site)

run_mcmc_results <- function (PROJECT,stage,repair) {

    print('Welcome to RUN_MCMC_RESULTS!!')

    # bundle needed functions down the chain
    functions_list=c("read_parameter_chains","read_binary_file_format","simulate_all",
                     "read_binary_response_surface","crop_development_parameters","have_chains_converged","psrf")
    # start marker
    stime=proc.time()["elapsed"]	

    # how many plots in total do we have
    nos_plots=1:PROJECT$nosites
    
    # now check which ones we need to calculate, but only if override not in play
    if (repair != 1) {
	print("...beginning filterings for sites we have already processed")
	keep_list=0
	for (i in seq(1, length(nos_plots))) {
	    outfile1=paste(PROJECT$results_processedpath,PROJECT$sites[i],"_parameters.RData",sep="")
	    if (file.exists(outfile1) == FALSE) {keep_list=append(keep_list,i)}
	}
	# filter out the sites we already have then
	keep_list=keep_list[-1] ; print(paste("......removing ",length(nos_plots)-length(keep_list)," sites out of ",length(nos_plots)," from the analysis",sep="")) 
	nos_plots=nos_plots[keep_list] 
    }

    # Combine driver modifiers into a list to be passed to the run_each_site()
    stage5modifiers=list(airt_factor = airt_factor, swrad_factor = swrad_factor, co2_factor = co2_factor,
                         rainfall_factor = rainfall_factor, wind_spd_factor = wind_spd_factor, vpd_factor = vpd_factor,
                         deforestation_factor = deforestation_factor, burnt_area_factor = burnt_area_factor)

    # now request the creation of the plots
    if (use_parallel & length(nos_plots) > 1) {
	print("...beginning parallel operations")
	cl <- makeCluster(min(length(nos_plots),numWorkers), type = "PSOCK")
	clusterExport(cl,functions_list)
	# load R libraries in cluster
	clusterExport(cl,"load_r_libraries") ; clusterEvalQ(cl, load_r_libraries())
	dummy=parLapply(cl,nos_plots,fun=run_each_site,PROJECT=PROJECT,stage=stage,repair=repair,grid_override=grid_override,stage5modifiers=stage5modifiers) 
	stopCluster(cl)	
    } else {
	print("...beginning serial operations")
	# or use serial
	dummy=lapply(nos_plots,FUN=run_each_site,PROJECT=PROJECT,stage=stage,repair=repair,grid_override=grid_override,stage5modifiers=stage5modifiers)
    } # parallel option
    # tell me whats happening
    print(paste("...time to process ",round((proc.time()["elapsed"]-stime)/60,1)," minutes",sep=""))

} # end function run_mcmc_results
## Use byte compile
run_mcmc_results<-cmpfun(run_mcmc_results)
