
###
## Function to describe the PROJECT requirements
###

cardamom_project_setup <- function (paths,PROJECT) {

      # local paths to be set
      typepath=paste(paths$cardamom_output,"/",PROJECT$type,"/",sep="")
      localpath=paste(paths$cardamom_output,"/",PROJECT$type,"/",PROJECT$name,"/",sep="")
      datapath=paste(localpath,"DATA/",sep="")
      resultspath=paste(localpath,"RESULTS/",sep="")
      results_processedpath=paste(localpath,"RESULTS_PROCESSED/",sep="")
      figpath=paste(localpath,"FIGURES/",sep="")
      exepath=paste(localpath,"EXECUTABLE/",sep="")

      # some useful variables
      modelname=PROJECT$model$name
      parameter_type=PROJECT$parameter_type
      project_src=PROJECT$source
      project_type=PROJECT$type

      # create local paths if they do no exist already
      failed=TRUE
      while(failed) {
	  yesno=readline("Do you want to automatically generate filepath directories (y/n)")
	  if (yesno != "y" & yesno != "n") {failed=TRUE} else {failed=FALSE}
      }
      if (yesno == "y") {
	  if (file.exists(typepath) == FALSE ){system(paste("mkdir ",typepath,sep=""))}
	  if (file.exists(localpath) == FALSE ){system(paste("mkdir ",localpath,sep=""))}
	  if (file.exists(datapath) == FALSE ){system(paste("mkdir ",datapath,sep=""))}
	  if (file.exists(resultspath) == FALSE ){system(paste("mkdir ",resultspath,sep=""))}
	  if (file.exists(results_processedpath) == FALSE ){system(paste("mkdir ",results_processedpath,sep=""))}
	  if (file.exists(figpath) == FALSE ){system(paste("mkdir ",figpath,sep=""))}
	  if (file.exists(exepath) == FALSE ){system(paste("mkdir ",exepath,sep=""))}
      }

      # number of chains desired?
      failed=TRUE
      while (failed) {
	  nochains=as.integer(readline("How many chains?"))
	  if (nochains > 1 & nochains < 8) {failed=FALSE} else {failed=TRUE} 
      }
      # number of samples per chain / accepted parameters
      failed=TRUE
      while(failed) {
	  nsamples=as.integer(readline("How many parameters to accept?"))
	  if (nsamples >= 1e4 & nsamples <= 30e6) {failed=FALSE} else {failed=TRUE} 
      }
      failed=TRUE
      while(failed) {
	  nsubsamples=as.integer(readline("How many accepted parameters to keep per chain (recommend = 1000)?"))
	  if (nsubsamples > 1e2 & nsubsamples < 1e6 & nsubsamples < nsamples) {failed=FALSE} else {failed=TRUE} 
      }

      # how much of sample chain to be kept
      latter_sample_frac=0.5
      # approximate chain runtime estimate in hours (approx 1 our for 1 million iterations)
      cre=nsamples/1.0e6

      # ask the user how long they want to set the simulation to run for
      chain_runtime=readline(paste("What is the maximum expected runtime for each chain (whole hours)? (NOTE: Given ",nsamples," required parameter vectors per chain, a 1000 timestep chain will take approximately ",cre," hours to run. However, allow for at least twice this time period if possible.)"))
      if (as.numeric(chain_runtime) > 48 | as.numeric(chain_runtime) < 0.5) {chain_runtime=readline(paste("Maximum number of hours to be submitted to Eddie is 48 hours, please re-select the number of hours"))}

      # PROJECT discription
      description=readline("Any other comments?")

      # calculate the sample rate
      samplerate=nsamples/nsubsamples

      # creation date
      date=Sys.time()
      # remove the spaces and other characters
      date=gsub("-", "",date)
      date=gsub(" ", "_",date)
      date=gsub(":", "_",date)

      # define executable name
      exe=paste(PROJECT$name,".exe",sep="")

      use_eddie=readline("Will you run this PROJECT on Eddie (y/n)")
      if (use_eddie == "y") {
	  use_eddie = TRUE 
	  # do we want an email to notify you of eddie works
	  email=readline("Enter your email address for ECDF notification (if you want)")
      } else { 
	  use_eddie = FALSE
      }

      # do I compile on eddie?
      if (use_eddie) {
	  # are we using eddie or not
	  eddiepath=paste(paths$cardamom_ecdf,project_type,PROJECT$name,sep="/")
	  ecdf_source=paste(paths$cardamom_ecdf,"LIBRARY/",sep="/")
	  # declare the eddie specific paths
	  edatapath=paste(eddiepath,"DATA/",sep="/")
	  eresultspath=paste(eddiepath,"RESULTS/",sep="/")
	  eoestreampath=paste(eddiepath,"OUTPUT_ERROR_STREAM/",sep="/")
	  eexepath=paste(eddiepath,"EXECUTABLES/",sep="/")

	  # check current host address
	  #home_computer=Sys.info()["nodename"]

	  # generate cardamom submit scripts
#	  generate_eddie_submit_script(paths)

	  # create directories on eddie and copy some important shell scripts
	  commands=c(paste("mkdir ",paths$cardamom_ecdf,sep="")
		    ,paste("mkdir ",ecdf_source,sep="")
		    ,paste("mkdir ",paths$cardamom_ecdf,"/",project_type,sep="")
		    ,paste("mkdir ",paths$cardamom_ecdf,"/",project_type,"/",PROJECT$name,sep="")
		    ,paste("mkdir ",edatapath,sep="")
		    ,paste("mkdir ",eresultspath,sep="")
		    ,paste("mkdir ",eoestreampath,sep="")
		    ,paste("mkdir ",eexepath,sep="")
		    ,paste("scp ",username,"@",home_computer,":",paths$cardamom,"/cardamom_functions/CARDAMOM_ECDF_SUBMIT_BUNDLES.sh ",eexepath,sep="")
		    ,paste("chmod +x ",eexepath,"/CARDAMOM_ECDF_SUBMIT_BUNDLES.sh",sep=""))

	  print("WARNING: Do not update the source code while jobs are currently running or they will fail")
	  comline=readline("Copy and compile any source code updates to Eddie (y/n)?")
	  if (comline == "y") {
	      print("Backup source code currently on eddie first")
	      print("Then copying source code to eddie")
	      print("Finally compile source code on eddie")
	      if (project_src == "C") {
		  commands=append(commands,c(paste("mv ",ecdf_source,"CARDAMOM_C ",ecdf_source,"CARDAMOM_C_BKP",sep="")
			    ,paste("scp -r ",username,"@",home_computer,":",paths$cardamom,"LIBRARY/CARDAMOM_C ",ecdf_source,sep="")
			    ,paste("gcc ",ecdf_source,"CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/DALEC_CDEA_TEMPLATE.c -o ",ecdf_source,"CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/a.out -lm",sep="")
			    ,paste("cp ",ecdf_source,"CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/a.out ",eexepath,"/",exe,sep="")))
	      } else if (project_src == "Fortran") {
		  # compiler options 
		  compiler_options=""#"-xhost -ipo -no-ftz"
		  if (timing) {compiler_options=paste(compiler_options," -pg",sep="")}
		  if (debug) {compiler_options=paste(compiler_options," -debug",sep="")}
		  commands=append(commands,c(paste("rm -r ",ecdf_source,"CARDAMOM_F_BKP",sep="")
			    ,paste("mv ",ecdf_source,"CARDAMOM_F ",ecdf_source,"CARDAMOM_F_BKP",sep="")
			    ,paste("scp -r ",username,"@",home_computer,":",paths$cardamom,"LIBRARY/CARDAMOM_F ",ecdf_source,sep="")
			    ,paste("cd ",ecdf_source,"CARDAMOM_F/executable",sep="")
			    ,paste("rm cardamom.exe")
			    ,paste("rm ",eexepath,"/",exe,sep="")
			    ,paste(compiler," -O2 ",compiler_options," ../misc/math_functions.f90 ../misc/oksofar.f90 ../model/",modelname,"/src/",modelname,".f90",
			    " ../model/",modelname,"/src/",modelname,"_CROP.f90",
			    " ../general/cardamom_structures.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90",
			    " ../model/",modelname,"/src/",modelname,"_PARS.f90 ../general/cardamom_io.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC.f90",
			    " ../model/",modelname,"/likelihood/MODEL_LIKELIHOOD.f90 ../general/cardamom_main.f90 -o cardamom.exe",sep="")
			    ,paste("cp ",ecdf_source,"CARDAMOM_F/executable/cardamom.exe ",eexepath,"/",exe,sep="")))
		  if (modelname == "AT_DALEC") {
		      commands=append(commands,paste("cp ",ecdf_source,"CARDAMOM_F/model/AT_DALEC/src/gpp_emulator_parameters* ",eexepath,"/",sep=""))
		      system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/model/AT_DALEC/src/gpp_emulator_parameters* ",exepath,"/",sep=""))
		  } #  only for AT_DALEC use the GAM emulator
		  if (modelname == "AT_DALEC" & parameter_type == "pft_specific") {
		      commands=append(commands,paste("cp ",ecdf_source,"CARDAMOM_F/model/AT_DALEC/src/winter_wheat_development.csv ",eexepath,"/",sep=""))
		      system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/model/AT_DALEC/src/winter_wheat_development.csv ",exepath,"/",sep=""))
		  } #  only for AT_DALEC use the GAM emulator
		  if (modelname == "DALEC_GSI_DFOL_CWD_FR" & parameter_type == "pft_specific") {
		      commands=append(commands,paste("cp ",ecdf_source,"CARDAMOM_F/model/DALEC_GSI_DFOL_CWD_FR/src/winter_wheat_development.csv ",eexepath,"/",sep=""))
		      system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/model/DALEC_GSI_DFOL_CWD_FR/src/winter_wheat_development.csv ",exepath,"/",sep=""))
		  } #  only for AT_DALEC use the GAM emulator
		  if (modelname == "DALEC_GSI_BUCKET" & parameter_type == "pft_specific") {
		      commands=append(commands,paste("cp ",ecdf_source,"CARDAMOM_F/model/DALEC_GSI_BUCKET/src/winter_wheat_development.csv ",eexepath,"/",sep=""))
		      system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/model/DALEC_GSI_BUCKET/src/winter_wheat_development.csv ",exepath,"/",sep=""))
		  } #  only for AT_DALEC use the GAM emulator
		  if (modelname == "DALECN_GSI_BUCKET" & parameter_type == "pft_specific") {
		      commands=append(commands,paste("cp ",ecdf_source,"CARDAMOM_F/model/DALEC_GSI_BUCKET/src/winter_wheat_development.csv ",eexepath,"/",sep=""))
		      system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/model/DALEC_GSI_BUCKET/src/winter_wheat_development.csv ",exepath,"/",sep=""))
		  } #  only for AT_DALEC use the GAM emulator
	      } else {
		  stop('Source code language has not been specified')
	      }
	  }
	  # issue commands to eddie
	  ecdf_execute(commands)

      } 
      # then compile locally
      comline="y" 
      if (comline == "y") {
	  print("Finally compile source code locally")
	  if (project_src == "C") {
	      system(paste("gcc ",paths$cardamom,"LIBRARY/CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/DALEC_CDEA_TEMPLATE.c -o ",paths$cardamom,"LIBRARY/CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/a.out -lm",sep=""))
	      system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/a.out ",exepath,"/",exe,sep=""))
	  } else if (project_src == "Fortran") {
	      # compiler options
              compiler_options=""#"-xhost -ipo -no-ftz"
	      if (timing) {compiler_options=paste(compiler_options," -pg",sep="")}
	      if (debug) {compiler_options=paste(compiler_options," -debug -traceback",sep="")}
		# if executables are present either in source library or in project folder remove them
		#if (file.exists(paste(paths$cardamom,"LIBRARY/CARDAMOM_F/executable/cardamom.exe",sep=""))) { system(paste("rm cardamom.exe")) }
		if (file.exists(paste(exepath,"/",exe,sep=""))) {system(paste("rm ",exepath,"/",exe,sep=""))}
		# store current working directory so that we can leave it briefly but return later
		cwd=getwd()
		setwd(paste(paths$cardamom,"LIBRARY/CARDAMOM_F/executable/",sep=""))
		# issue compile commands
		system(paste(compiler," -O2 ",compiler_options," ../misc/math_functions.f90 ../misc/oksofar.f90 ../model/",modelname,"/src/",modelname,".f90",
		" ../model/",modelname,"/src/",modelname,"_CROP.f90",
		" ../general/cardamom_structures.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90",
		" ../model/",modelname,"/src/",modelname,"_PARS.f90 ../general/cardamom_io.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC.f90",
		" ../model/",modelname,"/likelihood/MODEL_LIKELIHOOD.f90 ../general/cardamom_main.f90 -o cardamom.exe",sep=""))
		system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/executable/cardamom.exe ",exepath,"/",exe,sep=""))
		# now also generate the shared library needed later by R
		system(paste("gfortran -O2 -shared ../model/",modelname,"/src/",modelname,".f90 ",
		"../model/",modelname,"/src/",modelname,"_CROP.f90 ",
		"../model/",modelname,"/src/",modelname,"_R_interface.f90 ","-o dalec.so -fPIC",sep=""))
		system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/executable/dalec.so ",exepath,"/dalec.so",sep=""))
		if (modelname == "AT_DALEC") {
		    system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/model/AT_DALEC/src/gpp_emulator_parameters* ",exepath,"/",sep=""))
		} #  only for AT_DALEC use the GAM emulator
		if (modelname == "AT_DALEC" & parameter_type == "pft_specific") {
		    system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/model/AT_DALEC/src/winter_wheat_development.csv ",exepath,"/",sep=""))
		} #  only for AT_DALEC use the GAM emulator
		if ((modelname == "DALEC_GSI_DFOL_CWD_FR" | modelname == "DALEC_GSI_BUCKET" | modelname == "DALECN_GSI_BUCKET") & parameter_type == "pft_specific") {
		    system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/model/DALEC_GSI_DFOL_CWD_FR/src/winter_wheat_development.csv ",exepath,"/",sep=""))
		} #  only for AT_DALEC use the GAM emulator
		# return to original working directory
		setwd(cwd)
	  } else {
	      stop('Source code language has not been specified')
	  }
      } # copy and compile to eddie

      # prepare output
      PROJECT$ecdf=use_eddie
      PROJECT$localpath=localpath
      PROJECT$datapath=datapath
      PROJECT$resultspath=resultspath
      PROJECT$results_processedpath=results_processedpath
      PROJECT$figpath=figpath
      PROJECT$exepath=exepath
      PROJECT$nochains=nochains
      PROJECT$nsamples=nsamples
      PROJECT$latter_sample_frac=latter_sample_frac
      PROJECT$nsubsamples=nsubsamples
      PROJECT$chain_runtime=chain_runtime
      PROJECT$description=description
      PROJECT$samplerate=samplerate
      PROJECT$date=date
      PROJECT$exe=exe
	# eddie specific information
	if (use_eddie) {
	    PROJECT$eddiepath=eddiepath
	    PROJECT$edatapath=edatapath
	    PROJECT$eresultspath=eresultspath
	    PROJECT$eoestreampath=eoestreampath
	    PROJECT$eexepath=eexepath
	    PROJECT$email=email
	}
	return(PROJECT)
} # function end

