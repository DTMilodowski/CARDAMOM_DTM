
#include "../../auxi_fun/oksofar.c"

/*defines all the structures, i.e. DATA, MCOPT, PI*/

#include "../../mcmc_fun/MHMCMC/MCMC_FUN/MCMCOPT.c"
#include "MCMC_SETUP/MCMC_MODULES.c"


#include "../../mcmc_fun/MHMCMC/MCMC_FUN/MHMCMC.c"

/*MCMC subroutine*/
#include "MCMC_SETUP/MODEL_LIKELIHOOD.c"


int main(int argc,char *CLA[]){
/*To correctly set-up the MHMCMC*/

/*inputs*/
/*1. file in*/
/*2. file out*/
/*3. number of MCMC solutions requested*/
/*4. print-to-screen frequency*/
/*5. write-to-file frequency*/




/*SETTING number of command line imputs as char in CLA[0]*/
sprintf(CLA[0],"%d",argc-1);




/*defining data structure*/
DATA DATA;
/*defining parameter structure*/
PARAMETER_INFO PI;
/*defining MCMC_OPTIONS structure*/
MCMC_OPTIONS MCOPT;
/*defining the MODEL_LIKELIHOOD*/
/*MODEL_LIKELIHOOD MODEL;*/

MCMC_OUTPUT MCOUT;


/*These lines guarantee high frequency random generator seeding*/
struct timeval time;
gettimeofday(&time,NULL);
srand((time.tv_sec * 1e6) + (time.tv_usec));





/*Defining all MCMC components*/
/*USER DEFINED: SETUP MCMC - templates provides*/
READ_PARI_DATA(&PI, &DATA, &MCOUT, CLA);
READ_MCOPT(&MCOPT,CLA);

/*the function MODEL should be defined in the subroutines*/





MHMCMC(MODEL_LIKELIHOOD,DATA,PI,MCOPT,&MCOUT);

/*???????*/
/*User Defined function needed to clean up memory*/
MEMORY_CLEANUP(DATA,PI,MCOPT,MCOUT);

return 0;

}





