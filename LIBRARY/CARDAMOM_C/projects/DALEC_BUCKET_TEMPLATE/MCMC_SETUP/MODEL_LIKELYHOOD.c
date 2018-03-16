#pragma once
#include <math.h>
#include "../../DALEC_CODE/RC/RC1_BUCKET.c"
#include "../../DALEC_CODE/RC/RC2_BUCKET.c"
#include "../../DALEC_CODE/DALEC_BUCKET/DALEC_BUCKET.c"
#include "PROJECT_FUN/LIKELYHOOD.c"
#include "PROJECT_FUN/LIKELYHOOD_P.c"





double MODEL_LIKELYHOOD(DATA DATA,PARAMETER_INFO PI,double *PARS){

struct RCDIAGNOSTIC RCD;
RCD.DIAG=0;
int RC,n;
double P=0,P_p;

RC=pow(RC1_BUCKET(PARS,&RCD,DATA.meantemp,DATA.meanrad),DATA.RC);
P=P+log((double)RC);
if (RC==1){

/*PARAMETER LOG LIKELYHOOD*/
P=P+LIKELYHOOD_P(DATA,PARS);
P_p=P;

/*running model*/
DALEC_BUCKET(DATA.MET, PARS, DATA.deltat,DATA.nodays, DATA.LAT, DATA.M_LAI, DATA.M_NEE, DATA.M_FLUXES, DATA.M_POOLS);
/*storing GPP*/
for (n=0;n<DATA.nodays;n++){DATA.M_GPP[n]=DATA.M_FLUXES[n*16];}


/*RC2 check*/
RC=RC2_BUCKET(PARS, DATA.MET, DATA.M_LAI, DATA.M_NEE, DATA.M_POOLS,DATA.M_FLUXES,DATA.nodays,PI,&RCD,DATA.meantemp);
RC=pow(RC,DATA.RC);



/*LIKELYHOOD*/
P=P+log((double)RC);
if (RC==1){P=P+LIKELYHOOD(DATA);}}


/*Returning the log likelihood P*/
return P;


}




