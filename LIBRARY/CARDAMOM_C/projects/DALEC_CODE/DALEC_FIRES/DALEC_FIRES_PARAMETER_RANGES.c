
typedef struct{
double parmin[23];
double parmax[23];
int npars;}PARAMETER_INFO;

#include "mex.h"
#include "math.h"
#include "PARS_INFO_FIRES.c"





/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
 /* if(nrhs!=0) mexErrMsgTxt("no inputs required.");*/
  


/*defining parameter info structure ptr here*/
  PARAMETER_INFO PI;
  PARS_INFO_FIRES(&PI);
 
 /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(1,PI.npars, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,PI.npars, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
 

  
  /*  create a C pointer to a copy of the output matrix */
  double *nopars, *parmin,*parmax;
 
  nopars= mxGetPr(plhs[2]);
  parmin= mxGetPr(plhs[0]);
  parmax= mxGetPr(plhs[1]);
 
  nopars[0]=(double)PI.npars;


 int n;
  for (n=0;n<PI.npars;n++){parmin[n]=PI.parmin[n];parmax[n]=PI.parmax[n];}



}
