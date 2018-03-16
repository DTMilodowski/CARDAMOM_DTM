





#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../DALEC_CODE/DALEC_CDEA/PARS_INFO_CDEA.c"
#include "../../CARDAMOM_GENERAL/CARDAMOM_READ_BINARY_DATA.c"
#include "PROJECT_FUN/FIND_EDC_INITIAL_VALUES.c"
/*DALEC_SYNTHETIC SETUP*/







READ_PARI_DATA(PARAMETER_INFO *PI, DATA *DATA,MCMC_OUTPUT *MCOUT, char *CLA[]){
/*READING IN DALEC_SYNTHETIC DATA*/
/*opening file*/


/*CHANGE4EDDIE*/
/*this only applies to the native environment - change accordingly to add default file!*/
char filename[200];
if (atoi(CLA[0])<1){strcpy(filename,"/home/abloom/pdramatlab/dalec_fun/dalec_data_template.bin");}
else{strcpy(filename,CLA[1]);}
/*strcat(filename,CLA[1]);*/
/*strcat(filename,".bin");*/


printf("File %s\n",filename);
INITIALIZE_DATA_STRUCT(DATA);
CARDAMOM_READ_BINARY_DATA(filename,DATA);
oksofar("Site level data read ok");

/*Allocating memory for model output*/
/*SPECIFIC to DALEC_CDEA*/
int nopools=6;
DATA->nopools=nopools;
DATA->M_LAI=calloc(DATA->nodays,sizeof(double));
DATA->M_GPP=calloc(DATA->nodays,sizeof(double));
DATA->M_NEE=calloc(DATA->nodays,sizeof(double));
DATA->M_FLUXES=calloc(DATA->nodays*16,sizeof(double));
DATA->M_POOLS=calloc((DATA->nodays+1)*nopools,sizeof(double));

oksofar("Created fields for Model output");


/*PARAMETER INFO*/
/*READING DALEC INFO FROM FILE!!*/

/*defining parameter structure*/
/*Note - old structure quite similar to new, worth standardizing!*/
PI->npars=23;
PI->parmin=calloc(PI->npars,sizeof(double));
PI->parmax=calloc(PI->npars,sizeof(double));
PI->parini=calloc(PI->npars,sizeof(double));
PI->parfix=calloc(PI->npars,sizeof(double));
PI->stepsize=calloc(PI->npars,sizeof(double));


/*defining maximum and minimum parameter values*/
PARS_INFO_CDEA(PI);


/*defining step size*/
int n;
for (n=0;n<PI->npars;n++){PI->stepsize[n]=0.01;}

oksofar("Done with parameter definitions");

/*defining initial values*
 * Need to perform MCMC run to determine this*/
FIND_EDC_INITIAL_VALUES(*DATA,PI);



for (n=0;n<PI->npars;n++){PI->stepsize[n]=0.001;}
oksofar("Done with initial parameters");

INITIALIZE_MCMC_OUTPUT(*PI,MCOUT);




}


/*MCMC OPTIONS*/






READ_MCOPT(MCMC_OPTIONS *MCOPT, char *CLA[]){

/*number of command line imputs*/
int ncli=atoi(CLA[0]);

/*defining MCMC_OPTIONS structure*/
MCOPT->APPEND=1;
MCOPT->nADAPT=100;
MCOPT->fADAPT=0.5;
/*command line (or default) values*/
if (ncli<3){MCOPT->nOUT=1000;}else{MCOPT->nOUT=atoi(CLA[3]);};
if (ncli<4){MCOPT->nPRINT=1000;}else{MCOPT->nPRINT=atoi(CLA[4]);};
if (ncli<5){MCOPT->nWRITE=1000;}else{MCOPT->nWRITE=atoi(CLA[5]);};

MCOPT->randparini=0;
MCOPT->returnpars=0;
MCOPT->fixedpars=0;
char outfile[200], stepfile[200];
if (ncli<2){strcpy(outfile,"MOUT_");strcpy(stepfile,"MOUT_");}
else{strcpy(outfile,CLA[2]); strcpy(stepfile,CLA[2]);}
strcat(outfile,"PARS");
strcat(stepfile,"STEP");

/*directory*/
strcpy(MCOPT->outfile,outfile);
strcpy(MCOPT->stepfile,stepfile);
}


/*Enter all fields originally defined with MALLOC*/
int MEMORY_CLEANUP(DATA DATA, PARAMETER_INFO PI, MCMC_OPTIONS MCOPT, MCMC_OUTPUT MCOUT){

free(PI.parmin);
free(PI.parmax);
free(PI.parini);
free(PI.parfix);
free(PI.stepsize);
free(DATA.MET);
free(DATA.LAI);
free(DATA.NEE);
free(DATA.WOO);
free(DATA.GPP);


free(DATA.M_FLUXES);
free(DATA.M_LAI);
free(DATA.M_NEE);
free(DATA.M_POOLS);
free(DATA.M_GPP);

if (DATA.ngpp>0){free(DATA.gpppts);}
if (DATA.nlai>0){free(DATA.laipts);}
if (DATA.nnee>0){free(DATA.neepts);}
if (DATA.nwoo>0){free(DATA.woopts);}


free(MCOUT.best_pars);


return 0;

}








