
###
## Function to submit processes to eddie
###

submit_processes_to_local_machine<-function (PROJECT_in) {

    print('PREPARING TO SUBMIT MCMC TO LOCAL MACHINE')
    print('This function should be valid for all CARDAMOM compatible DALEC MCMC functions')

    # do we want to remove any previous output files?
    delete_old=readline("Delete any previous output files for this project name?(y/n)")    
    if (delete_old == "y") {system(paste("rm ",PROJECT_in$resultspath,"/*",sep=""))}
    background=readline("Run jobs in background?(y/n)")
    # begin writing out the file contents
    # construct the file now
    cwd=getwd()
    setwd(PROJECT_in$exepath)
    for (n in seq(1, PROJECT_in$nosites)) {
	for (c in seq(1, PROJECT_in$nochains)) {
	    infile=paste(PROJECT_in$datapath,PROJECT_in$name,"_",PROJECT_in$sites[n],".bin",sep="")
	    output=paste(PROJECT_in$resultspath,PROJECT_in$name,"_",PROJECT_in$sites[n],"_",c,"_",sep="")
	    if (background == "y") {
		system(paste(PROJECT_in$exepath,PROJECT_in$exe," ",infile," ",output," ",as.integer(PROJECT_in$nsamples)," 0 ",as.integer(PROJECT_in$samplerate)," & ",sep=""))
		#print(paste(PROJECT_in$exepath,PROJECT_in$exe," ",infile," ",output," ",as.integer(PROJECT_in$nsamples)," 0 ",as.integer(PROJECT_in$samplerate)," & ",sep=""))
	    } else {
		system(paste(PROJECT_in$exepath,PROJECT_in$exe," ",infile," ",output," ",as.integer(PROJECT_in$nsamples)," 0 ",as.integer(PROJECT_in$samplerate),sep=""))
	    }
	}
    }
    setwd(cwd) ; rm(cwd)
    print("Command issued to local machine")

} # end of function
