#include <math.h>
#include "../../../DALEC_CODE/RC/RC1_FIRES.c"
#include "../../../DALEC_CODE/RC/RC2_FIRES.c"
#include "../../../DALEC_CODE/DALEC_FIRES/DALEC_FIRES.c"
#include "../MODEL_LIKELYHOOD.c"




double RC_MODEL_LIKELYHOOD(DATA DATA_IN,PARAMETER_INFO PI, double *PARS){

struct RCDIAGNOSTIC RCD;
RCD.DIAG=1;
int RC, n;
double P=0;

DATA *DATA;DATA=&DATA_IN;

RC=RC1_FIRES(PARS,&RCD, DATA->meantemp,DATA->meanrad);

/*running model*/
DALEC_FIRES(DATA->MET, PARS, DATA->deltat,DATA->nodays, DATA->LAT, DATA->M_LAI, DATA->M_NEE, DATA->M_FLUXES, DATA->M_POOLS);






/*RC2 check*/
RC=RC*RC2_FIRES(PARS, DATA->MET, DATA->M_LAI, DATA->M_NEE, DATA->M_POOLS,DATA->M_FLUXES,DATA->nodays,PI,&RCD,DATA->meantemp);

/*LIKELYHOOD*/
int tot_exp=0;
for (n=0;n<RCD.nrc;n++){
	tot_exp+=1-RCD.PASSFAIL[n];
}

P=-0.5*((double)tot_exp*10)*DATA->RC;


/*overriding if model likelyhood is zero*/
double ML=MODEL_LIKELYHOOD(*DATA,PI,PARS);
if (DATA->RC==0 & (isinf(ML)==-1 || isnan(ML))){P=P-0.5*10;}

/*for (n=0;n<RCD.nrc;n++){printf("%d ",RCD.PASSFAIL[n]);};printf("\n");*/




/*Log likelihood*/
return P;

}




