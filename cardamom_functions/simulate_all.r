
###
## Function to generate mean state variable information by running the parameters and model choice
###

simulate_all<- function (site,PROJECT,model_name,met,pars,lat,pft,parameter_type,exepath,soil_info) {

      noedc=100
      output_dim=17 ; aNPP_dim=8

      # declare output variables
      # order is npar,chain,step
      laiall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
      gppall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
      rtotall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
      soilevapall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
      if (model_name != "ACM") {
	  neeall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  somall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  bioall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  rooall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  litall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  laball=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  folall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  rautoall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  rhetall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  harall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  gsiall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  gsitempall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  gsiphotoall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  gsivpdall=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  Clabslow=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  somfast=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  litroot=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  litwood=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  microact=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  microbial=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  litN=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  labN=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
	  DIN=array(0,dim=c(dim(pars)[2]*dim(pars)[3],dim(met)[1]))
      }

      # restructure pars
      pars_in=array(0,dim=c(dim(pars)[1],dim(pars)[2]*dim(pars)[3]))
      nos_iter=dim(pars)[2]*dim(pars)[3]
      for (n in seq(1,dim(pars)[1])){
	  pars_in[n,]=pars[n,,]
      }
      # loop through combinations
      if (model_name == "ACM") {
	  dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
	  output_dim=5
	  if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
	  tmp=.Fortran( "racm",output_dim=as.integer(output_dim),met=as.double(t(met)),pars=as.double(pars_in),out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
			              ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
	        	              ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
			              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
				      ,soil_frac_clay=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
				      ,soil_frac_sand=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4))) )
	  output=tmp$out_var
          output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
	  dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALEC_GSI_BUCKET") {
	  output_dim=23
          dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
	  crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
          if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
          tmp=.Fortran( "rdalecgsibucket",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                      ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                      ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                      ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
				      ,soil_frac_clay=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
				      ,soil_frac_sand=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4)))
				      ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
          output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
          aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
          dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALECN_GSI_BUCKET") {
	  output_dim=23
	  if (is.loaded("rdalecngsibucket") == FALSE) { dyn.load(paste(PROJECT$exepath,"/dalec.so", sep="")) }
	  crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
          if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
          tmp=.Fortran("rdalecngsibucket",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                      ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                      ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                      ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
				      ,soil_frac_clay=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
				      ,soil_frac_sand=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4)))
				      ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
          output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
          aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
          if (site == PROJECT$sites[length(PROJECT$sites)]) {dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))}
          rm(tmp) ; gc()
      } else if (model_name == "DALEC_CDEA") {
	  dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
	  pft_specific = 0
	  tmp=.Fortran( "rdaleccdea",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
				      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
			              ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
	        	              ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
			              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
	  output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
	  aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
	  dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALEC_CDEA_FR") {
	  dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
	  if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
	  tmp=.Fortran( "rdaleccdeafr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
				      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
			              ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
	        	              ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
			              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
	  output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
	  aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
	  dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALEC_GSI_FR_LABILE") {
	  output_dim=18
	  dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
	  if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
	  tmp=.Fortran( "rdalecgsifr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
				      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
			              ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
	        	              ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
			              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
	  output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
	  aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
	  dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALECN_GSI_FR") {
	  output_dim=22
	  dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
	  if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
	  tmp=.Fortran( "rdalecngsifr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
				      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
			              ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
	        	              ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
			              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
	  output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
	  aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
	  dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALEC_GSI_FR") {
	  dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
	  if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
	  tmp=.Fortran( "rdalecgsifr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
				      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
			              ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
	        	              ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
			              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
	  output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
	  aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
	  dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALEC_GSI_DFOL_FR") {
          dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
          if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
          tmp=.Fortran( "rdalecgsidfolfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
                                      ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                      ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                      ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
          output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
          aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
          dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALEC_GSI_DFOL_CWD_FR") {
          output_dim=19
          dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
	  crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
          if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
          tmp=.Fortran( "rdalecgsidfolcwdfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
                                      ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                      ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                      ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
				      ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
          output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
          aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
          dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALEC_GSI_DFOL_FROOT_FR") {
	  output_dim=19
          dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
          if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
          tmp=.Fortran( "rdalecgsidfolfrootfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
                                      ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                      ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                      ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
          output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
          aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
          dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALECN_GSI_DFOL_LABILE_FR") {
	  output_dim=24
          dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
          if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
          tmp=.Fortran( "rdalecngsidfollabfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
                                      ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                      ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                      ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
          output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
          aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
          dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {
	  output_dim=24
          dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
          if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
          tmp=.Fortran( "rdalecngsidfollabfrootfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
                                      ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                      ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                      ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
          output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
          aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
          dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALEC_GSI_DFOL_LABILE_FR") {
	  output_dim=20
          dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
          if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
          tmp=.Fortran( "rdalecgsidfollabfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
                                      ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                      ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                      ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
          output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
          aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
          dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALEC_GSI_MFOL_FR") {
	  dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
	  if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
	  tmp=.Fortran( "rdalecgsimfolfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
					 ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
			                 ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                         ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
	        	                 ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
			                 ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
	  output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
	  aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
	  dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "DALEC_GSI_DBio_FR") {
          output_dim=22
          dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
          if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
          tmp=.Fortran( "rdalecgsibiofr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
				      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
                                      ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                      ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                      ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                      ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
	  output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
	  aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
          dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else if (model_name == "AT_DALEC") {
	  randomForest=read_binary_response_surface(paste(exepath,"/gpp_emulator_parameters_",pft,".bin",sep=""))
	  crop_info=crop_development_parameters(paste(exepath,"/winter_wheat_development.csv",sep=""))
	  dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
	  if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
	  tmp=.Fortran( "ratdalec",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
				  ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(dim(pars)[2]*dim(pars)[3],aNPP_dim)))
			          ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                  ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
	        	          ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
			          ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),dim1=as.integer(randomForest$dim_1),dim2=as.integer(randomForest$dim_2)
                                  ,cdim=as.integer(crop_info$crop_dims),dummy_nos_trees=as.integer(randomForest$ntree)
			          ,dummy_nos_inputs=as.integer(randomForest$nos_inputs),nos_iter=as.integer(nos_iter)
                                  ,stock_seed_labile=as.double(crop_info$stock_seed_labile)
	     		          ,tmp1=as.double(randomForest$leftDaughter),tmp2=as.double(randomForest$rightDaughter)
                                  ,tmp3=as.double(randomForest$nodestatus),tmp4=as.double(randomForest$xbestsplit)
                                  ,tmp5=as.double(randomForest$nodepred),tmp6=as.double(randomForest$bestvar)
			          ,DS_shoot=as.double(crop_info$DS_shoot),fol_frac=as.double(crop_info$fol_frac),stem_frac=as.double(crop_info$stem_frac)
                                  ,DS_root=as.double(crop_info$DS_root),root_frac=as.double(crop_info$root_frac)
                                  ,DS_LRLV=as.double(crop_info$DS_LRLV),LRLV=as.double(crop_info$LRLV)
                                  ,DS_LRRT=as.double(crop_info$DS_LRRT),LRRT=as.double(crop_info$LRRT) )
	  output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
	  aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
	  dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
          rm(tmp) ; gc()
      } else {
	  stop(paste("Model choice (",model_name,") does not have corresponding R interface",sep=""))
      }
      # pass output (iter,nodays)
      laiall=output[,,1] #
      gppall=output[,,2] # gCm-2day-1
      evapall=output[,,3]# kgH2Om-2.day-1
      if (model_name == "ACM" & output_dim == 5) {soilevapall=output[,,4] ; rtotall=output[,,5]}
      if (model_name != "ACM") {
	  rautoall=output[,,3]#
	  rhetall=output[,,4]#
	  neeall=output[,,5]#
	  woodall=output[,,6]#
	  somall=output[,,7]#
	  bioall=output[,,8]#
	  rooall=output[,,9]#
	  litall=output[,,10]#
	  laball=output[,,11]#
	  folall=output[,,12]#
	  harall=output[,,13]#
	  if (grepl("GSI",model_name)) {
	      gsiall=output[,,14] # GSI combined
	      gsitempall=output[,,15]# GSI air temperature component
	      gsiphotoall=output[,,16] # GSI photoperiod
	      gsivpdall=output[,,17] # GSI VPD component
	  }
	  if (model_name == "DALEC_GSI_BUCKET" | model_name == "DALECN_GSI_BUCKET") {
	      evapall=output[,,18]   # kgH2Om-2.day-1}
	      rootwater=output[,,19] # kgH2Om-2 root zone
              wSWP=output[,,20] # Estimates total hydraulic resistance (MPa)
	      litwood=output[,,21]
              fire=output[,,23]
	  }
          if (model_name == "DALEC_GSI_DFOL_CWD_FR") {litwood=output[,,18] ; fire=output[,,19]}
	  if (model_name == "DALEC_GSI_FR_LABILE") {labile_slow=output[,,18]}
	  if (model_name == "DALEC_GSI_DFOL_FROOT_FR") {labroot=output[,,18] ; labwood=output[,,19]}
	  if (model_name == "DALEC_GSI_DFOL_LABILE_FR") {labroot=output[,,18] ; labwood=output[,,19] ; litwood=output[,,20]}
	  if (model_name == "DALECN_GSI_DFOL_LABILE_FR" | model_name == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {
	      labroot=output[,,18] ; labwood=output[,,19] ; litwood=output[,,20]
     	      litN=output[,,21] ; labN=output[,,22] ; DIN=output[,,23]
	      N_mineralisation=output[,,24]
	  }
	  if (model_name == "DALEC_GSI_DBio_FR") {
	    somfast=output[,,18]
	    litroot=output[,,19]
	    litwood=output[,,20]
	    microact=output[,,21]
	    microbial=output[,,22]
	  }
	  if (model_name == "DALECN_GSI_FR") {
	      litwood=output[,,18]
	      litN=output[,,19]
	      labN=output[,,20]
	      DIN=output[,,21]
	      N_mineralisation=output[,,22]
	  }
      } # not acm

	if (model_name != "ACM") {
	    # combine outputs
	    states_all=list(lai=laiall
			    ,gpp=gppall
			    ,nee=neeall
			    ,reco=neeall+gppall
			    ,rauto=rautoall
			    ,rhet=rhetall
			    ,wood=woodall
			    ,som=somall
			    ,bio=bioall
			    ,root=rooall
			    ,lit=litall
			    ,lab=laball
			    ,fol=folall
			    ,harvest_C=harall
			    ,gsi=gsiall
			    ,gsi_itemp=gsitempall
			    ,gsi_iphoto=gsiphotoall
			    ,gsi_ivpd=gsivpdall
			    ,aNPP=aNPP)
      #		      ,states_mean=states_mean)

	    if (model_name == "DALEC_GSI_FR_LABILE") {states_all$Clabslow=labile_slow}
	    if (model_name == "DALEC_GSI_DFOL_FROOT_FR") {states_all$labroot=labroot ; states_all$labwood=labwood}
            if (model_name == "DALEC_GSI_DFOL_CWD_FR") {states_all$litwood=litwood ; states_all$fire=fire}
	    if (model_name == "DALEC_GSI_DFOL_LABILE_FR") {states_all$labroot=labroot ; states_all$labwood=labwood ; states_all$litwood=litwood}
	    if (model_name == "DALECN_GSI_DFOL_LABILE_FR" | model_name == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {
		states_all$labroot=labroot ; states_all$labwood=labwood
		states_all$litwood=litwood ; states_all$litN=litN
		states_all$labN=labN ; states_all$DIN=DIN
		states_all$N_mineralisation=N_mineralisation
	    }

	    # something species for DBio model at the moment
	    if (model_name == "DALEC_GSI_DBio_FR") {
		states_all$somfast=somfast
		states_all$litroot=litroot
		states_all$litwood=litwood
		states_all$microact=microact
		states_all$microbial=microbial
	    }
	    if (model_name == "DALECN_GSI_FR") {
	      states_all$litwood=litwood
	      states_all$litN=litN
	      states_all$labN=labN
	      states_all$DIN=DIN
	      states_all$N_mineralisation=N_mineralisation
	    }
          if (model_name == "DALEC_GSI_BUCKET" | model_name == "DALECN_GSI_BUCKET") {
	      states_all$litwood=litwood
              states_all$evap=evapall   # kgH2Om-2.day-1}
              states_all$rootwater=rootwater # kgH2Om-2 root zone
              states_all$wSWP=wSWP
              states_all$fire=fire
          }
      } else {
	    # combine outputs
	    states_all=list(lai=laiall,gpp=gppall,evap=evapall,soilevap=soilevapall,Rtot = rtotall)
      }

      #print("NOTE: some variables are missing prior to completion of matlab transcription")
      # return state variable means
      return(states_all) ; gc(verbose=FALSE)

} # end of function
## Use byte compile
simulate_all<-cmpfun(simulate_all)
