
typedef struct{
double parmin[23];
double parmax[23];
int npars;}PARAMETER_INFO;

#include "mex.h"
#include "math.h"
#include "DALEC_CDEA.c"
#include "PARS_INFO_CDEA.c"

#include "../EDC/EDC1_CDEA.c"
#include "../EDC/EDC2_CDEA.c"




/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *consts,*pars,*MET,*LAI,*FLUXES,*POOLS,*deltat,*NEE,*lat,*EDC,*EDCID;
  mwSize nr,nc;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=4) mexErrMsgTxt("Five inputs (MET,pars,deltat,lat) required.");
  
  /* check to make sure the first input argument is a scalar 
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      mxGetN(prhs[0])*mxGetM(prhs[0])!=1 ) {
    mexErrMsgTxt("Input x must be a scalar.");
  }*/
  
  
  /*function in MATLAB form 
 * [LAI,FLUXES,POOLS] = DALEC_CDE_C (MET,pars,deltat)*/
  MET = mxGetPr(prhs[0]);
 
/*pars = 17 of them*/
 pars = mxGetPr(prhs[1]);

  deltat = mxGetPr(prhs[2]);
  lat = mxGetPr(prhs[3]);

/*number of rows/cols in input matrix*/
  nr = mxGetN(prhs[0]);
  nc = mxGetM(prhs[0]);

/*defining constants here*/


int nopools=6;



/*defining EDC diagnostic vector*/
/*1s and 0s for each EDC*/
struct EDCDIAGNOSTIC EDCD;
EDCD.nedc=100;

  
 /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(nr,1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(nr,1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(16,nr, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(nopools,nr+1, mxREAL);
  plhs[4] = mxCreateDoubleMatrix(1,1, mxREAL);
  plhs[5] = mxCreateDoubleMatrix(1,EDCD.nedc, mxREAL);
 
  
  /*  create a C pointer to a copy of the output matrix */
  /*GPP = mxGetPr(plhs[0]);*/
  LAI= mxGetPr(plhs[0]);
  NEE= mxGetPr(plhs[1]);
  FLUXES = mxGetPr(plhs[2]);
  POOLS = mxGetPr(plhs[3]);
  EDC = mxGetPr(plhs[4]);
  EDCID = mxGetPr(plhs[5]);
  /*  call the C subroutine */
  /*gppresponse(pars,consts,GPP);*/
  
  /*assigning values to pools*/
  /*LR,L,F,R,W,Lit,SOM*/
 
int n;


 DALEC_CDEA(MET,pars,deltat[0],nr,lat[0],LAI,NEE,FLUXES,POOLS); 


int m;



  /*pointer to pointer to double values*/
 


PARAMETER_INFO parinfo;
PARS_INFO_CDEA(&parinfo);

for (n=0;n<EDCD.nedc;n++){EDCD.PASSFAIL[n]=1;}


/*mean temp and mean rad*/
double meantemp=0, meanrad=0;
for(n=0;n<nr;n++){
meantemp+=(MET[n*6+1]*0.5+MET[n*6+2]*0.5)/(double)nr;
meanrad+=MET[n*6+3]/(double)nr;}



 EDCD.DIAG=1;
EDC[0]=1;
 EDC[0]=EDC[0]*EDC1_CDEA(pars,&EDCD,meantemp,meanrad);


 EDC[0]=EDC[0]*(double)EDC2_CDEA(pars,MET,LAI,NEE,POOLS, FLUXES,nr,parinfo,&EDCD,meantemp);
 /*assigning all PASSFAIL values to EDCID*/
 if (EDCD.DIAG==1){for (n=0;n<EDCD.nedc;n++){EDCID[n]=(double)EDCD.PASSFAIL[n];}}




}
