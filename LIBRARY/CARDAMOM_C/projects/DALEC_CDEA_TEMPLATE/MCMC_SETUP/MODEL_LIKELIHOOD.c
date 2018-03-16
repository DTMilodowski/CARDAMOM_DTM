#pragma once
#include <math.h>
#include "../../DALEC_CODE/EDC/EDC1_CDEA.c"
#include "../../DALEC_CODE/EDC/EDC2_CDEA.c"
#include "../../DALEC_CODE/DALEC_CDEA/DALEC_CDEA.c"
#include "PROJECT_FUN/LIKELIHOOD.c"
#include "PROJECT_FUN/LIKELIHOOD_P.c"





double MODEL_LIKELIHOOD(DATA DATA,PARAMETER_INFO PI,double *PARS){

struct EDCDIAGNOSTIC EDCD;
EDCD.DIAG=0;
int EDC,n;
double P=0,P_p;

EDC=pow(EDC1_CDEA(PARS,&EDCD,DATA.meantemp,DATA.meanrad),DATA.EDC);
P=P+log((double)EDC);
if (EDC==1){

/*PARAMETER LOG LIKELIHOOD*/
P=P+LIKELIHOOD_P(DATA,PARS);
P_p=P;

/*running model*/
DALEC_CDEA(DATA.MET, PARS, DATA.deltat,DATA.nodays, DATA.LAT, DATA.M_LAI, DATA.M_NEE, DATA.M_FLUXES, DATA.M_POOLS);
/*storing GPP*/
for (n=0;n<DATA.nodays;n++){DATA.M_GPP[n]=DATA.M_FLUXES[n*16];}


/*EDC2 check*/
EDC=EDC2_CDEA(PARS, DATA.MET, DATA.M_LAI, DATA.M_NEE, DATA.M_POOLS,DATA.M_FLUXES,DATA.nodays,PI,&EDCD,DATA.meantemp);
EDC=pow(EDC,DATA.EDC);



/*LIKELIHOOD*/
P=P+log((double)EDC);
if (EDC==1){P=P+LIKELIHOOD(DATA);}}


/*Returning the log likelihood P*/
return P;


}




