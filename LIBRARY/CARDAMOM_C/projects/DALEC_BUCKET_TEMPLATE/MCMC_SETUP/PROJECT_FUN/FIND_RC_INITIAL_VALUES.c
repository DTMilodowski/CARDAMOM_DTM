#include "../../../../mcmc_fun/MHMCMC/MCMC_FUN/MHMCMC.c"
#include "RC_MODEL_LIKELYHOOD.c"

int FIND_RC_INITIAL_VALUES(DATA DATA,PARAMETER_INFO *PI){


/*This MCMC is designed to find the best fitting parameters ONLY*/

MCMC_OPTIONS MCOPT;
MCMC_OUTPUT MCOUT;
INITIALIZE_MCMC_OUTPUT(*PI,&MCOUT);


MCOPT.APPEND=0;
MCOPT.nADAPT=20;
MCOPT.fADAPT=0.5;
MCOPT.nOUT=2000;
MCOPT.nPRINT=0;
MCOPT.nWRITE=0;
/*randparini = 0*/
/*this means all PI.parini values must either be given values or entered as -9999*/
MCOPT.randparini=1;
MCOPT.returnpars=1;
/*setting fixedpars option to 1*/
MCOPT.fixedpars=1;


/*status update here*/
oksofar("starting MCMC for RC inipars");
int n;




for (n=0;n<PI->npars;n++){
PI->stepsize[n]=0.02;
/*PI->stepsize[n]=0.0005;*/
PI->parini[n]=DATA.parpriors[n];
PI->parfix[n]=0;
if (PI->parini[n]!=-9999 & DATA.rc_random_search<1) {PI->parfix[n]=1;}

}



/*done*/
/*PRC is the log likelihood*/
double PRC=log(0);
int count=0;
while (PRC<0){
oksofar("RC attempt");
for (n=0;n<PI->npars;n++){PI->stepsize[n]=0.0005;}
/*insert prior value option here!*/
MHMCMC(RC_MODEL_LIKELYHOOD,DATA,*PI,MCOPT,&MCOUT);
for (n=0;n<PI->npars;n++){PI->parini[n]=MCOUT.best_pars[n];}

PRC=RC_MODEL_LIKELYHOOD(DATA,*PI, PI->parini);count=count+1;
/*in case one RC missing*/
if (PRC<0 & count%3==0){PI->parini[n]=DATA.parpriors[n];}

}


/*for (n=0;n<PI->npars;n++){printf("%8.6f  ",PI->parini[n]);}printf("\n");
for (n=0;n<PI->npars;n++){PI->stepsize[n]=0.01;}*/
/*resetting fixed pars to zero for main r*/
for (n=0;n<PI->npars;n++){PI->parfix[n]=0;}

/*clearing MCOUT fields*/
free(MCOUT.best_pars);



}

