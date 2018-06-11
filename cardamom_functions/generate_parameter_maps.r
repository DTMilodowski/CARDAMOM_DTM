
###
## Function to extract state variables and direct the production of uncertainty plots for the key states and fluxes
###

generate_parameter_maps<-function(PROJECT) {


	# how many years of the analysis
	nos_years=length(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
	# load timesteps to local variable
	timestep_days=PROJECT$model$timestep_days ; seconds_per_day = 86400
	if (length(timestep_days) == 1) {
	    eh = TRUE ; n = 0
	    while(eh) {
		n = n + 1
		# generate file name of the output file created in stage 3
		loadfile=paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")
		if (file.exists(loadfile)) {load(loadfile) ; eh = FALSE}
	    }
	    timestep_days=rep(timestep_days, length.out=drivers$nodays)
	}

	# determine new output file for aggregated values
	outfile=paste(PROJECT$results_processedpath,PROJECT$name,"_parameter_maps.RData",sep="")

	if (file.exists(outfile) == FALSE | repair == 1) {
	    # loop through all the files to eventually build up a complete vector
	    for (n in seq(1, PROJECT$nosites)) {
		# generate file name of the output file created in stage 3
		loadfile=paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")

		if (file.exists(loadfile) == TRUE) {
		    load(loadfile) #; print(paste("DALEC simulations will be loaded from ",loadfile,sep=""))
                    if (length(which(is.na(as.vector(aNPP)) == TRUE)) > 0) {print(paste("site with NA is ",loadfile,sep=""))}
		    if (n < 100 | n%%100 == 0) {print(paste("...have loaded ",round((n/PROJECT$nosites)*100, digits=0),"% of pixels",sep=""))}
		    # create array we will be filling here if we have not already
		    if (exists("par_array_median") == FALSE) {
			max_pars=PROJECT$nochains*PROJECT$nsubsamples*0.5
			par_array_median=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,(max(PROJECT$model$nopars)+1)))
			par_array_unc=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,(max(PROJECT$model$nopars)+1)))
			par_array_converged=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,(max(PROJECT$model$nopars)+1)))
			pixel_initial_age_array=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
			pixel_yield_class_array=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
			pft_array=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
			# some allocation and residence time variables that are dependent on model states
			aNPP_array_median=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,3))
			aNPP_array_unc=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,3))
			resid_time_array_median=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,5))
			resid_time_array_unc=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,5))
			mean_temperature_array=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
			mean_vpd_array=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
			mean_radiation_array=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
			mean_precipitation_array=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
                        # extract non-model values
                        obs_lai_max=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
                        obs_lai_mean=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
                        obs_lai_sigma=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
                        obs_wood_estimate=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
                        obs_wood_uncertainty=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
		    }
		    # now load parameter information into the array
		    slot_j=as.numeric(PROJECT$sites[n])/PROJECT$long_dim
		    slot_i=as.numeric(PROJECT$sites[n])-(floor(slot_j)*PROJECT$long_dim)
		    if(slot_i == 0){slot_i = PROJECT$long_dim} ; slot_j=ceiling(slot_j)
		    mean_temperature_array[slot_i,slot_j]=mean((drivers$met[,2]+drivers$met[,3])*0.5)
		    mean_radiation_array[slot_i,slot_j]=mean(drivers$met[,4])
		    mean_vpd_array[slot_i,slot_j]=mean(drivers$met[,12])
		    mean_precipitation_array[slot_i,slot_j]=sum(drivers$met[,7]*timestep_days*seconds_per_day)/nos_years
                    obs_lai_max[slot_i,slot_j]=max(drivers$obs[,2])
                    obs_lai_mean[slot_i,slot_j]=mean(drivers$obs[which(drivers$obs[,2] != -9999),2])
                    obs_lai_sigma[slot_i,slot_j]=sd(drivers$obs[which(drivers$obs[,2] != -9999),2])
                    if (length(which(drivers$obs[,7] > -9999)) > 0) {
                        obs_wood_estimate[slot_i,slot_j]=max(drivers$obs[which(drivers$obs[,7] != -9999),7],na.rm=TRUE)
                        obs_wood_uncertainty[slot_i,slot_j]=max(drivers$obs[which(drivers$obs[,7] != -9999),18],na.rm=TRUE)
                    }
		    par_array_converged[slot_i,slot_j,]=0
		    converged=have_chains_converged(parameters)
		    par_array_converged[slot_i,slot_j,which(converged=="PASS")]=1
		    pixel_initial_age_array[slot_i,slot_j]=drivers$age
		    pixel_yield_class_array[slot_i,slot_j]=drivers$yield
		    pft_array[slot_i,slot_j]=drivers$ctessel_pft
		    # calculate some allocation (fol,wood,root) 
		    bob = quantile(aNPP[,1], prob=c(0.975,0.50,0.025),na.rm=TRUE)
		    aNPP_array_median[slot_i,slot_j,1]=bob[2] ; aNPP_array_unc[slot_i,slot_j,1]=bob[1]-bob[3]
		    bob = quantile(aNPP[,2], prob=c(0.975,0.50,0.025),na.rm=TRUE)
		    aNPP_array_median[slot_i,slot_j,2]=bob[2] ; aNPP_array_unc[slot_i,slot_j,2]=bob[1]-bob[3]
		    bob = quantile(aNPP[,3], prob=c(0.975,0.50,0.025),na.rm=TRUE)
		    aNPP_array_median[slot_i,slot_j,3]=bob[2] ; aNPP_array_unc[slot_i,slot_j,3]=bob[1]-bob[3]
		    # residence time variables (fol,wood,root,som,DeadOrg)
		    bob = quantile(aNPP[,4], prob=c(0.975,0.50,0.025),na.rm=TRUE)
		    resid_time_array_median[slot_i,slot_j,1]=bob[2] ; resid_time_array_unc[slot_i,slot_j,1]=bob[1]-bob[3]
		    bob = quantile(aNPP[,5], prob=c(0.975,0.50,0.025),na.rm=TRUE)
		    resid_time_array_median[slot_i,slot_j,2]=bob[2] ; resid_time_array_unc[slot_i,slot_j,2]=bob[1]-bob[3]
		    bob = quantile(aNPP[,6], prob=c(0.975,0.50,0.025),na.rm=TRUE)
		    resid_time_array_median[slot_i,slot_j,3]=bob[2] ; resid_time_array_unc[slot_i,slot_j,3]=bob[1]-bob[3]
		    bob = quantile(aNPP[,7], prob=c(0.975,0.50,0.025),na.rm=TRUE)
		    resid_time_array_median[slot_i,slot_j,4]=bob[2] ; resid_time_array_unc[slot_i,slot_j,4]=bob[1]-bob[3]
		    if (PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR" | PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALECN_GSI_BUCKET") {
			bob = quantile(aNPP[,8], prob=c(0.975,0.50,0.025),na.rm=TRUE)
			resid_time_array_median[slot_i,slot_j,5]=bob[2] ; resid_time_array_unc[slot_i,slot_j,5]=bob[1]-bob[3]
		    }
		    for (p in seq(1, dim(parameters)[1])) {
			bob = quantile(as.vector(parameters[p,,]), prob=c(0.975,0.50,0.025))
			par_array_median[slot_i,slot_j,p] = bob[2]
			par_array_unc[slot_i,slot_j,p] = abs(bob[1]-bob[3])
		    }

		} # 
	    } #  site loop

	    # tmp hack for avgN log10-normal
	    if (PROJECT$model$name == "DALECN_GSI_FR" | PROJECT$model$name == "DALECN_GSI_DFOL_FR" 
	      | PROJECT$model$name == "DALEC_GSI_BUCKET"| PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR"
	      | PROJECT$model$name == "DALECN_GSI_BUCKET") {
		par_array_median[,,11]=10**par_array_median[,,11]
		par_array_unc[,,11]=10**par_array_unc[,,11]
	    }

	} else {
	    # if we have loaded these before just load the final product
	    load(outfile)
	} # have these previously been loaded

        # inform the user
        print("......have finished loading - now beginning cluster analysis")

	if (file.exists(outfile) == FALSE | repair == 1) {
	    nos_uk_clusters = 1 ; uk_cluster = 1 ; uk_cluster_pft = 1
	    if (PROJECT$model$name == "DALEC_GSI_DFOL_FR" | PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR" | PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALECN_GSI_BUCKET") {

		# remove non-constrained parameters (i.e. those not actually used in this analysis)
		initial_conditions=c(18:23)
		par_array_median_normalised=par_array_median[,,-c(initial_conditions,28,29,30,31,32,33,35,36)]
		if (PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR") {
		    par_array_median_normalised=par_array_median[,,-c(initial_conditions,28,29,30,31,32,33,35,36,37)]
		} else if (PROJECT$model$name == "DALEC_GSI_BUCKET") {
		    par_array_median_normalised=par_array_median[,,-c(initial_conditions,28,29,30,31,32,33,35,36,37)]
		} else if (PROJECT$model$name == "DALECN_GSI_BUCKET") {
		    par_array_median_normalised=par_array_median[,,-c(initial_conditions,28,29,30,31,32,33,35,36,37)]
		}
		# now normalise the dataset
		for (i in seq(1,dim(par_array_median_normalised)[3])) {
		    min_par_val=min(par_array_median[,,i],na.rm=TRUE)
		    max_par_val=max(par_array_median[,,i],na.rm=TRUE)
		    par_array_median_normalised[,,i]=((par_array_median[,,i]-min_par_val)/(max_par_val-min_par_val))
		}
		#summary(as.vector(par_array_median_normalised)) ; hist(as.vector(par_array_median_normalised))
		par_array_tmp=array(NA,dim=c(prod(dim(par_array_median)[1:2]),dim(par_array_median_normalised)[3]))
		par_array_tmp[1:prod(dim(par_array_median_normalised)[1:2]),1:dim(par_array_median_normalised)[3]]=par_array_median_normalised
		actual_forests=which(is.na(par_array_tmp[,1]) == FALSE)
		par_array_tmp=par_array_tmp[actual_forests,]
		par_array_tmp=array(par_array_tmp,dim=c((length(par_array_tmp)/dim(par_array_median_normalised)[3]),dim(par_array_median_normalised)[3]))

		tmp = 0
                # Looping to find preference_input, the preferenceRange() returns 2 values, the first of which minises the number of clusters, 
                # while the seconds would return as many clusters as there are observations. 
                # It is the responsibility of the user to ensure the most appropriate use of these information to result in an appropriate number of clusters for error propagation
#		for (i in seq(1,100)) {
#		    tmp=append(tmp,preferenceRange(negDistMat(par_array_tmp[sample(1:dim(par_array_tmp)[1],0.05*dim(par_array_tmp)[1], replace=FALSE),],r=2))[1])
#		} ; preference_input=max(mean(tmp[-1]),median(tmp[-1]))
#		uk_clusters=apclusterL(negDistMat(r=2),par_array_tmp,frac=0.1,sweeps=10, p=preference_input,maxits=1000, convits=100)
                uk_clusters=apclusterL(negDistMat(r=2),par_array_tmp,frac=0.1,sweeps=10, q=0.5,maxits=1000, convits=100)
		nos_uk_clusters=length(uk_clusters@clusters) ; uk_clusters_exemplars=uk_clusters@exemplars

		uk_cluster_pft=array(NA,dim=c(dim(par_array_median_normalised)[1:2]))
		for (i in seq(1,length(uk_clusters@clusters))) {
		    uk_cluster_pft[actual_forests[uk_clusters@clusters[[i]]]] = i
		}
		uk_cluster_pft=array(uk_cluster_pft,dim=c(dim(par_array_median_normalised)[1:2]))
	    }
	}

        # some tidying which I need to check on later - i.e. THIS IS A HACK
        resid_time_array_unc[which(abs(as.vector(resid_time_array_unc)) == abs(log(0)) )]=NA
        resid_time_array_median[which(abs(as.vector(resid_time_array_median)) == abs(log(0)) )]=NA
	# inform the user
	print("......now generating parameter maps")

	# describe parameter numbers with log liklihood added on the end
	character_bit=rep("p",times=(PROJECT$model$nopars[1]))
	number_bit=1:(PROJECT$model$nopars[1])
	# merge the letter and numbers together
	par_names=c(paste(character_bit,number_bit,sep=""),"log-likelihood")
	# determine correct height and widths
	hist_height=4000 ; hist_width=7200
	fig_height=4000 ; fig_width=7200
	#fig_height=4000 ; fig_width=7200 ; proj="+init=epsg:4326"
	if (PROJECT$grid_type == "UK") { fig_height=8000 ; fig_width=7200 }
	# load colour palette
	colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral"))) 

	# calculate land mask
	#long1=seq(PROJECT$longitude[1],PROJECT$longitude[2], length.out=PROJECT$long_dim)
	#lat1=seq(PROJECT$latitude[1],PROJECT$latitude[2], length.out=PROJECT$lat_dim)
	landmask=array(PROJECT$landsea, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

        # assuming we have generated one lets create the Cluster analysis map
        if (PROJECT$model$name == "DALEC_GSI_DFOL_FR" | PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR" | PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALECN_GSI_BUCKET") {
            jpeg(file=paste(PROJECT$figpath,"Cluster_map_of_median_parameters_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
            par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
            image.plot(uk_cluster_pft, main=paste("Cluster analysis potential PFT map",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
            contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
            dev.off()
        }

	# create tempory array of parameters for the covariance analysis
	jpeg(file=paste(PROJECT$figpath,"parameter_correlations_median_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	tmp=array(par_array_median[,,1:(dim(par_array_median)[3])],dim=c(prod(dim(par_array_median)[1:2]),(dim(par_array_median)[3])))
	image.plot(cor(tmp,use="na.or.complete", method="spearman"), main="Median Parameter Correlation (Spearmans)")
	dev.off()
	# create tempory array of parameters for the covariance analysis
	jpeg(file=paste(PROJECT$figpath,"parameter_correlations_uncertainty_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	tmp=array(par_array_unc[,,1:(dim(par_array_unc)[3])],dim=c(prod(dim(par_array_unc)[1:2]),(dim(par_array_unc)[3])))
	image.plot(cor(tmp,use="na.or.complete", method="spearman"), main="Parameter Uncertainty Correlation (Spearmans)")
	dev.off()

	###
	## Some things that are really parameters but are dependent on the evolution of state variables so need more information to be calculated
	# assign this value once as aall arrays have the same number of values present
	colour_choices=colour_choices_upper(length(aNPP_array_median[,,1]))

	# create map of NPP allocations
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Cfoliar_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(aNPP_array_median[,,1], col=colour_choices, main=paste("Cfoliar NPP allocation median estimate",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Cwood_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(aNPP_array_median[,,2],col=colour_choices, main=paste("Cwood NPP allocation median estimate",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Croot_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(aNPP_array_median[,,3], col=colour_choices, main=paste("Croot NPP allocation median estimate",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_uncertainty_","Cfoliar_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(aNPP_array_unc[,,1], col=colour_choices, main=paste("Cfoliar NPP allocation uncertainty",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_uncertainty_","Cwood_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(aNPP_array_unc[,,2], col=colour_choices, main=paste("Cwood NPP allocation uncertainty",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_uncertainty_","Croot_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(aNPP_array_unc[,,3], col=colour_choices, main=paste("Croot NPP allocation uncertainty",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()

	# create map of residence times
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Cfoliar_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(resid_time_array_median[,,1], col=colour_choices, main=paste("Cfoliar residence time median estimate",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Cwood_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(resid_time_array_median[,,2], col=colour_choices, main=paste("Cwood residence time median estimate",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Croot_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(resid_time_array_median[,,3], col=colour_choices, main=paste("Croot residence time median estimate",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Csom_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(resid_time_array_median[,,4], col=colour_choices, main=paste("Csom residence time median estimate",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()	
	if (PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR" | PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALECN_GSI_BUCKET") {
	    jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","CDeadOrg_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    image.plot(resid_time_array_median[,,5], col=colour_choices, main=paste("CDeadOrg residence time median estimate",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	    contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	    dev.off()
	}
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_uncertainty_","Cfoliar_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(resid_time_array_unc[,,1], col=colour_choices, main=paste("Cfoliar residence time uncertainty",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,quantile(as.vector(resid_time_array_unc[,,1]),prob=c(0.975),na.rm=TRUE)))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_uncertainty_","Cwood_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(resid_time_array_unc[,,2], col=colour_choices, main=paste("Cwood residence time uncertainty",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,quantile(as.vector(resid_time_array_unc[,,2]),prob=c(0.975),na.rm=TRUE)))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_uncertainty_","Croot_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(resid_time_array_unc[,,3], col=colour_choices, main=paste("Croot residence time uncertainty",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,quantile(as.vector(resid_time_array_unc[,,3]),prob=c(0.975),na.rm=TRUE)))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	if (PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR" | PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALECN_GSI_BUCKET") {
	    jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_uncertainty_","CDeadOrg_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    image.plot(resid_time_array_unc[,,5], col=colour_choices, main=paste("CDeadOrg residence time uncertainty",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,quantile(as.vector(resid_time_array_unc[,,5]),prob=c(0.975),na.rm=TRUE)))
	    contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	    dev.off()
	}
	# create histograms
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_","Cfoliar_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(aNPP_array_median[,,1]), col=colour_choices, main=paste("Cfoliar NPP allocation median estimate",sep=""),   cex.main=2.4      )
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_","Cwood_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(aNPP_array_median[,,2]), col=colour_choices, main=paste("Cwood NPP allocation median estimate",sep=""),   cex.main=2.4      )
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_","Croot_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(aNPP_array_median[,,3]), col=colour_choices, main=paste("Croot NPP allocation median estimate",sep=""),   cex.main=2.4      )
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_uncertainty_","Cfoliar_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(aNPP_array_unc[,,1]), col=colour_choices, main=paste("Cfoliar NPP allocation uncertainty",sep=""),   cex.main=2.4      )
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_uncertainty_","Cwood_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(aNPP_array_unc[,,2]), col=colour_choices, main=paste("Cwood NPP allocation uncertainty",sep=""),   cex.main=2.4      )
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_uncertainty_","Croot_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(aNPP_array_unc[,,3]), col=colour_choices, main=paste("Croot NPP allocation uncertainty",sep=""),   cex.main=2.4      )
	dev.off()
	# create map of residence times
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_","Cfoliar_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(resid_time_array_median[,,1]), col=colour_choices, main=paste("Cfoliar residence time median estimate",sep=""),   cex.main=2.4      )
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_","Cwood_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(resid_time_array_median[,,2]), col=colour_choices, main=paste("Cwood residence time median estimate",sep=""),   cex.main=2.4      )
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_","Croot_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(resid_time_array_median[,,3]), main=paste("Croot residence time median estimate",sep=""),   cex.main=2.4      )
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_","Csom_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(resid_time_array_median[,,4]), col=colour_choices, main=paste("Csom residence time median estimate",sep=""),   cex.main=2.4      )
	dev.off()
	if (PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR" | PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALECN_GSI_BUCKET") {
	    jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_","CDeadOrg_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    hist(as.vector(resid_time_array_median[,,5]), col=colour_choices, main=paste("CDeadOrg residence time median estimate",sep=""),   cex.main=2.4      )
	    dev.off()
	}
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_uncertainty_","Cfoliar_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(resid_time_array_unc[,,1]), col=colour_choices, main=paste("Cfoliar residence time uncertainty",sep=""),   cex.main=2.4      )
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_uncertainty_","Cwood_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(resid_time_array_unc[,,2]), col=colour_choices, main=paste("Cwood residence time uncertainty",sep=""),   cex.main=2.4      )
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_uncertainty_","Croot_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(resid_time_array_unc[,,3]), col=colour_choices, main=paste("Croot residence time uncertainty",sep=""),   cex.main=2.4      )
	dev.off()
	jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_uncertainty_","Csom_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	hist(as.vector(resid_time_array_unc[,,4]), col=colour_choices, main=paste("Csom residence time uncertainty",sep=""),   cex.main=2.4      )
	dev.off()
	if (PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR" | PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALECN_GSI_BUCKET") {
	    jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_uncertainty_","CDeadOrg_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    hist(as.vector(resid_time_array_unc[,,5]), col=colour_choices, main=paste("CDeadOrg residence time uncertainty",sep=""),   cex.main=2.4      )
	    dev.off()
	}

	# now everything has been loaded into a nice array we begin plotting
	for (p in seq(1,dim(par_array_median)[3])) {

	    ###
	    ## spatial maps
	    jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_",par_names[p],"_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    image.plot(par_array_median[,,p], col=colour_choices, main=paste(par_names[p]," median estimate",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	    contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	    dev.off()
	    jpeg(file=paste(PROJECT$figpath,"parameter_maps_percent_uncertainty_",par_names[p],"_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    image.plot(abs((par_array_unc[,,p]/par_array_median[,,p])*100), col=colour_choices, main=paste(par_names[p]," uncertainty range % of median",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	    contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	    dev.off()
	    jpeg(file=paste(PROJECT$figpath,"parameter_maps_absolute_uncertainty_",par_names[p],"_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    image.plot(abs(par_array_unc[,,p]), col=colour_choices, main=paste(par_names[p]," uncertainty range",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	    contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	    dev.off()
	    jpeg(file=paste(PROJECT$figpath,"parameter_maps_converged_",par_names[p],"_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    image.plot(par_array_converged[,,p], col=colour_choices, main=paste(par_names[p]," Gelmen-Rubens convergence (1 = PASS / 0 = FALSE)",sep=""),axes=FALSE, zlim=c(0,1), cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
	    contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	    dev.off()
	    ###
	    ## histrograms
	    jpeg(file=paste(PROJECT$figpath,"parameter_hist_median_",par_names[p],"_",PROJECT$name,".jpg",sep=""), width=hist_width, height=hist_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    hist(as.vector(par_array_median[,,p]), main=paste(par_names[p]," median estimate",sep=""))
	    dev.off()
	    jpeg(file=paste(PROJECT$figpath,"parameter_hist_percent_uncertainty_",par_names[p],"_",PROJECT$name,".jpg",sep=""), width=hist_width, height=hist_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    hist(as.vector(abs((par_array_unc[,,p]/par_array_median[,,p])*100)), main=paste(par_names[p]," uncertainty range % of median",sep=""))
	    dev.off()
	    jpeg(file=paste(PROJECT$figpath,"parameter_hist_absolute_uncertainty_",par_names[p],"_",PROJECT$name,".jpg",sep=""), width=hist_width, height=hist_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    hist(as.vector(abs(par_array_unc[,,p])), main=paste(par_names[p]," uncertainty range",sep=""))
	    dev.off()
	    jpeg(file=paste(PROJECT$figpath,"parameter_hist_converged_",par_names[p],"_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    hist(as.vector(par_array_converged[,,p]), main=paste(par_names[p]," Gelmen-Rubens convergence (1 = PASS / 0 = FALSE)",sep=""))
	    dev.off()

	}

	# now finally for the age map
	if (length(which(as.vector(pixel_initial_age_array) > 0))) {
	    jpeg(file=paste(PROJECT$figpath,"age_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    image.plot(pixel_initial_age_array, col=colour_choices, main="Initial age (Years)",axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,max(as.vector(pixel_initial_age_array),na.rm=TRUE)))
	    contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	    dev.off()
	}
	# now finally for the yield class
	if (length(which(as.vector(pixel_yield_class_array) > 0))) {
	    jpeg(file=paste(PROJECT$figpath,"yield_class_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    image.plot(pixel_yield_class_array, col=colour_choices, main="Yield Class (m3.yr-1)",axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,max(as.vector(pixel_yield_class_array),na.rm=TRUE)))
	    contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	    dev.off()
	}
	# now any pft information we might have
	if (length(which(is.na((as.vector(pft_array))) == FALSE & as.vector(pft_array) > 0)) > 0) {
	    jpeg(file=paste(PROJECT$figpath,"pft_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	    image.plot(pft_array, main="CTESSEL PFT",axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,max(as.vector(pft_array),na.rm=TRUE)))
	    contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	    dev.off()
	}

	# mean temperature
	jpeg(file=paste(PROJECT$figpath,"mean_temperature_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(mean_temperature_array, col=colour_choices, main="Mean air temperature (oC)",axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(min(as.vector(mean_temperature_array),na.rm=TRUE),max(as.vector(mean_temperature_array),na.rm=TRUE)))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	# mean radiation
	jpeg(file=paste(PROJECT$figpath,"mean_radiation_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(mean_radiation_array, col=colour_choices, main="Mean radiation (MJ.m-2.day-1)",axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,max(as.vector(mean_radiation_array),na.rm=TRUE)))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	# mean vpd
	jpeg(file=paste(PROJECT$figpath,"mean_vpd_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(mean_vpd_array, col=colour_choices, main="Mean VPD (Pa)",axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,max(as.vector(mean_vpd_array),na.rm=TRUE)))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()
	# mean precipitation
	jpeg(file=paste(PROJECT$figpath,"mean_precipitation_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
	par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
	image.plot(mean_precipitation_array, col=colour_choices, main="Mean precipitation (mm.yr-1)",axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,max(as.vector(mean_precipitation_array),na.rm=TRUE)))
	contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
	dev.off()

	# output some aggragated values
	# probably best to add some aggregated met drivers to this concoction here
	save(landmask,nos_uk_clusters,uk_cluster,uk_cluster_pft,uk_clusters_exemplars
	    ,mean_precipitation_array,mean_temperature_array,mean_radiation_array,mean_vpd_array
	    ,par_array_median,par_array_unc,par_array_converged
            ,obs_lai_max,obs_lai_mean,obs_lai_sigma,obs_wood_estimate,obs_wood_uncertainty
	    ,pixel_yield_class_array,pixel_initial_age_array,pft_array,aNPP_array_median,aNPP_array_unc
	    ,resid_time_array_median,resid_time_array_unc,file=outfile)#, compress="gzip")

	# tidy before leaving
	gc(reset=TRUE, verbose=FALSE)

	print("...done with parameters...")

} # end function generate_parameter_maps
## Use byte compile
generate_parameter_maps<-cmpfun(generate_parameter_maps)
