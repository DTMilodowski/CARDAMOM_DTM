#include <math.h>
#include "../../../DALEC_CODE/EDC/EDC1_CDEA.c"
#include "../../../DALEC_CODE/EDC/EDC2_CDEA.c"
#include "../../../DALEC_CODE/DALEC_CDEA/DALEC_CDEA.c"
#include "../MODEL_LIKELIHOOD.c"




double EDC_MODEL_LIKELIHOOD(DATA DATA_IN,PARAMETER_INFO PI, double *PARS){

struct EDCDIAGNOSTIC EDCD;
EDCD.DIAG=1;
int EDC, n;
double P=0;

DATA *DATA;DATA=&DATA_IN;

EDC=EDC1_CDEA(PARS,&EDCD, DATA->meantemp,DATA->meanrad);

/*running model*/
DALEC_CDEA(DATA->MET, PARS, DATA->deltat,DATA->nodays, DATA->LAT, DATA->M_LAI, DATA->M_NEE, DATA->M_FLUXES, DATA->M_POOLS);






/*EDC2 check*/
EDC=EDC*EDC2_CDEA(PARS, DATA->MET, DATA->M_LAI, DATA->M_NEE, DATA->M_POOLS,DATA->M_FLUXES,DATA->nodays,PI,&EDCD,DATA->meantemp);

/*LIKELIHOOD*/
int tot_exp=0;
for (n=0;n<EDCD.nedc;n++){
	tot_exp+=1-EDCD.PASSFAIL[n];
}




P=-0.5*((double)tot_exp*10)*DATA->EDC;


/*overriding if model likelyhood is zero*/
double ML=MODEL_LIKELIHOOD(*DATA,PI,PARS);
if (DATA->EDC==0 & (isinf(ML)==-1 || isnan(ML))){P=P-0.5*10;}


/*adding an exponential decay related term to find starting point quicker!*/
double PC=0,C;
double co=-log(2)/(DATA->nodays*DATA->deltat);
for (n=0;n<6;n++){
C=expdecay2(DATA->M_POOLS,n,DATA->nodays+1,DATA->deltat,DATA->nopools);
if (C<co & finite(C)){PC=PC-0.5*pow((C-co)/(co*10),2);}}

P=P+PC;


/*Log likelihood*/
return P;

}




