#include "mex.h"
#include "math.h"
#include "DALEC_BUCKET.c"

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *consts,*pars,*MET,*LAI,*FLUXES,*POOLS,*deltat,*NEE,*lat,*RC;
  mwSize nr,nc;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=4) mexErrMsgTxt("Four inputs (MET,pars,deltat,lat) required.");
  
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






  
 /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(nr,1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(nr,1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(16,nr, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(7,nr+1, mxREAL);
  plhs[4] = mxCreateDoubleMatrix(1,1, mxREAL);
 
  
  /*  create a C pointer to a copy of the output matrix */
  /*GPP = mxGetPr(plhs[0]);*/
  LAI= mxGetPr(plhs[0]);
  NEE= mxGetPr(plhs[1]);
  FLUXES = mxGetPr(plhs[2]);
  POOLS = mxGetPr(plhs[3]);
  RC = mxGetPr(plhs[4]);
  /*  call the C subroutine */
  /*gppresponse(pars,consts,GPP);*/
  
  /*assigning values to pools*/
  /*LR,L,F,R,W,Lit,SOM*/
  
  /*CHECKING ALL INPUTS HERE*/

  
  DALEC_BUCKET(MET,pars,deltat[0],nr,lat[0],LAI,NEE,FLUXES,POOLS); 



}
