int WRITE_RESULTS(double *PARS, double PROB,PARAMETER_INFO PI,MCMC_OPTIONS MCO){

int n;

FILE *fileout=fopen(MCO.outfile,"ab");
FILE *filestep=fopen(MCO.stepfile,"wb");
for (n=0;n<PI.npars;n++){
        fwrite(&PARS[n],1,sizeof(double),fileout);
       	fwrite(&PI.stepsize[n],1,sizeof(double),filestep);}
    
/*writing likelyhood*/
fwrite(&PROB,1,sizeof(double),fileout);
fclose(fileout);
fclose(filestep);


return 0;


}
