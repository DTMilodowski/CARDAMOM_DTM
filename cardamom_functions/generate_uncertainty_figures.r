
###
## Function to extract state variables and direct the production of uncertainty plots for the key states and fluxes
###

generate_uncertainty_figures<-function(PROJECT,n) {

	# load locally needed library
	require(gplots)

	# generate file name of the output file created in stage 3
	loadfile=paste(PROJECT$results_processedpath,PROJECT$sites[n],".RData",sep="")

	if (file.exists(loadfile) == TRUE) {
	    #stime=proc.time()["elapsed"]
	    load(loadfile) ; print(paste("DALEC simulations will be loaded from ",loadfile,sep=""))
	    #print(paste("load dalec in ",proc.time()["elapsed"]-stime," seconds",sep=""))
	} else { 
	    # do we run the parameters yet for analysis
	    run_all=readline("Raw results have not been processed therefore we will do it now. Do you want to run all parameter vectors to generate confidence intervals? (y/n)")
	    if (run_all == "y") {
		PROJECT$latter_sample_frac=0.5 # readline("What (latter) fraction of accepted parameters to use (e.g. 0.5)?")
		run_mcmc_results(PROJECT)
	    } # if condition
	} # file exists statement

	# how many plots in total do we have
	nos_plots=0:12
        if (PROJECT$model$name == "DALEC_GSI_BUCKET") {nos_plots=c(-5,-4,-3,-2,nos_plots,15,22)}
	if (PROJECT$model$name == "DALECN_GSI_BUCKET") {nos_plots=c(-5,-4,-3,-2,nos_plots,15,22)}
#	if (PROJECT$model$name == "DALEC_GSI_FR_LABILE") {nos_plots=-1:12}
        if (PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR") {nos_plots=c(nos_plots,15,22)}
        if (PROJECT$model$name == "DALEC_GSI_DBio_FR") {nos_plots=0:16}
        if (PROJECT$model$name == "DALECN_GSI_FR" | PROJECT$model$name == "DALECN_GSI_DFOL_LABILE_FR" | PROJECT$model$name == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {nos_plots=c(0:12,15,17:21)}
	if (PROJECT$model$name == "ACM") {nos_plots=c(2,-2)}
	# now request the creation of the plots
	if (use_parallel & length(nos_plots) > 1) {
	    cl <- makeCluster(min(length(nos_plots),numWorkers), type = "PSOCK")    
	    # load R libraries in cluster
	    clusterExport(cl,c("load_r_libraries","rmse","gsi_controlling"))
	    clusterEvalQ(cl, load_r_libraries())
	    dummy=parLapply(cl,nos_plots,fun=uncertainty_figures,PROJECT=PROJECT,states_all=states_all,drivers=drivers,parameters=parameters,sub_parameter=sub_parameter,n=n,plotconfidence=plotconfidence) 
	    stopCluster(cl)	
	} else {
	    # or use serial
	    dummy=lapply(nos_plots,FUN=uncertainty_figures,PROJECT=PROJECT,states_all=states_all,drivers=drivers,parameters=parameters,sub_parameter=sub_parameter,n=n,plotconfidence=plotconfidence)
	} # parallel option

	# tidy before leaving
	gc(reset=TRUE, verbose=FALSE)

} 
## Use byte compile
generate_uncertainty_figures<-cmpfun(generate_uncertainty_figures)
