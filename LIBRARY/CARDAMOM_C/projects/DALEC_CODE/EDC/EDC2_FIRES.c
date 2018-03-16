#pragma once
#include "../../../math_fun/std.c"




/*Component of expdecay (next function)*/
int expdecay2(double const *POOLS, int pool, int nodays,int deltat,int nopools)
/*only accepting pool and number of days, the rest done here*/
{
int noyears=((nodays*deltat)/365)+1, EDC,n,count=0,y=0;
double *PY;


PY=calloc(noyears,sizeof(double));
double P0=POOLS[pool];

/*aggregating to number of years*/
/*calculating mean for every year*/
for (n=0;n<nodays;n++){
PY[y]=PY[y]+POOLS[n*nopools+pool];
if (count>=364/deltat){PY[y]=PY[y]/(count+1); count=-1;y=y+1;}
count=count+1;
}

/*performing exponential test*/
/*only with full years*/


/*CONTINUE FROM HERE*/
EDC=numexpdecay(P0,PY,y);
/*EDC=expdecay(PY,y)*a simpler version with gradients and percentages*/


free(PY);
return EDC;
}



int numexpdecay(double const P0,double const *PY,int const y)
{

/*compare and sync with dalec_expdecay.m*/

double n,a,b,c,A,B,C;
int EDC=1;
/*denominators*/
double dty=floor((double)y/2);
double dt0=y-floor((double)y/2)-1;
/*diff numerators/denominators*/
double dpdty=(PY[y-1]-PY[y-y/2-1])/(dty);
double dpdt0=(PY[y-y/2-1]-PY[0])/(dt0);





c=log(dpdty/dpdt0)/(floor((double)y/2));
b=(dpdt0*dt0)/(exp(c*(floor((double)y/2)))-1);
a=PY[0]-b;


C=c/365;
A=a;
B=P0-A;

/*half life must be less than 10 years*/

if (fabs(-log(2)/C)<y*365.25 & C<0 & finite(C) & abs(B)>abs(A)*0.01){EDC=0;}

/*C = a +b e ^{ct}*/
/*turnover in years*/
/* dCdt = b*c*exp(c*t)*/
/*log(dcdt(2)/dcdt(1)) = c(t2-t1)*/
/*note that t is t=[0,1,2,etc.]*/
/*
printf("dpdty = %f\n",dpdty);
printf("dpdt0 = %f\n",dpdt0);
printf("a = %f\n",a);
printf("b = %f\n",b);
printf("c = %f\n",c);
printf("A = %f\n",A);
printf("B = %f\n",B);
printf("C = %f\n",C);
printf("P0 = %f\n",P0);
printf("PY[0] = %f\n",PY[0]);
printf("EDC = %i\n",EDC);
*/



return EDC;


}


/*mean matrix from double pointer routine*/
double mean_pool(double *PA,int p,int nc, int nopools)
{
int r,c;
double meanpool=0;

for (c=0;c<nc;c++){
meanpool=meanpool+PA[c*nopools+p]/(double)nc;
}


return meanpool;

}


int EDC2_FIRES(double *pars, double *MET, double *LAI, double *NEE, double *POOLS,double *FLUXES,int nodays,PARAMETER_INFO PI,struct EDCDIAGNOSTIC *EDCD,double meantemp)
{



/*THESE EDC2 checks are for DALEC_FIRES*/
int EDC=1,n=0,edc=0;
int faili=0;
int DIAG=EDCD->DIAG;/*1 or 0*/


/*FIRES*/
int nomet=6,nopools=6;
int nofluxes=16;
int done=0;
int deltat=(int)(MET[nomet]-MET[0]);
int k=0;



/*deriving mean pools here!*/
double *MPOOLS;
MPOOLS=malloc(nopools*sizeof(double));
for (n=0;n<nopools;n++){MPOOLS[n]=mean_pool(POOLS,n,nodays+1,nopools);};



/***********************EDCs start here****************************/



/*EDC no 5:*/
/*MAX(Cfol)/clma <LAImax*/
faili=7;
double LAImax=10;
if (EDC==1 || DIAG==1){
for (n=0;n<nodays;n++){
if (LAI[n]>LAImax){EDC=0;EDCD->PASSFAIL[5-1]=0;}}}
/*MLAI=mean_pointer_vector(LAI,nodays);*/



/*EDC no 8*/
/*0.2*Cf < Cr < 5*Cf*/
/*Cfoliar : Croot = 5:1 or 1:5*/
if ((EDC==1 || DIAG==1) & (MPOOLS[1]>MPOOLS[2]*5 | MPOOLS[1]*5<MPOOLS[2])){EDC=0;EDCD->PASSFAIL[8-1]=0; }


/*EDC 9*/
/*G10 = 2 - growth factor in 10 years*/
/*for any further constraint please incorporate in likelyhood*/
double noyears=nodays*deltat/365.25;
double G=1.1;
for (n=0;n<nopools;n++){
/*Rapid POOL Growth: cannot increase by more than factor of Gn over N yrs*/
/*note - this is GROWTH ONLY!!!*/
if ((EDC==1 || DIAG==1) & (POOLS[n + nopools*nodays]/POOLS[n]>G*(double)noyears)){EDC=0;EDCD->PASSFAIL[9-1]=0;}}



/*EDC no 10*/
/*POOL EXPONENTIAL DECAY*/
/*Performed on all seven pools*/
for (n=0;n<nopools;n++){if (EDC==1 || DIAG==1) {
EDC=expdecay2(POOLS,n,nodays+1,deltat,nopools);
if (EDC==0){EDCD->PASSFAIL[10-1]=0;}}};




double const fauto=pars[1];
double const ffol=(1-fauto)*pars[2];
double const flab=(1-fauto-ffol)*pars[12];
double const froot=(1-fauto-ffol-flab)*pars[3];
double const fwood=1-fauto-ffol-flab-froot;

/*fraction of GPP som under equilibrium conditions*/
double const fsom=fwood+(froot+flab+ffol)*pars[0]/(pars[0]+pars[7]);
double const flit=(froot+flab+ffol);



/*SOM attractor - must be within a factor of 2 from Csom0*/
/*half the TR is balanced by x2 Csom*/

/*equilibrium factor (in comparison to C_initial)*/
double EQF=10.0;
double meangpp=0;
for (n=0;n<nodays;n++){meangpp+=FLUXES[n*nofluxes]/(double)nodays;}




/*EDC 11 - SOM quasi-steady-state*/
if ((EDC==1 || DIAG==1) & meangpp*fsom/(pars[8]*exp(pars[9]*meantemp))>pars[22]*EQF){EDC=0;EDCD->PASSFAIL[11-1]=0;}
if ((EDC==1 || DIAG==1) & meangpp*fsom/(pars[8]*exp(pars[9]*meantemp))<pars[22]/EQF){EDC=0;EDCD->PASSFAIL[11-1]=0;}

/*EDC 12 - LIT quasi-steady-state*/
if ((EDC==1 || DIAG==1) & meangpp*flit/(pars[7]*exp(pars[9]*meantemp))>pars[21]*EQF){EDC=0;EDCD->PASSFAIL[12-1]=0;}
if ((EDC==1 || DIAG==1) & meangpp*flit/(pars[7]*exp(pars[9]*meantemp))<pars[21]/EQF){EDC=0;EDCD->PASSFAIL[12-1]=0;}

/*EDC 13 - WOO quasi-steady-state*/
if ((EDC==1 || DIAG==1) & meangpp*fwood/(pars[5])>pars[20]*EQF){EDC=0;EDCD->PASSFAIL[13-1]=0;}
if ((EDC==1 || DIAG==1) & meangpp*fwood/(pars[5])<pars[20]/EQF){EDC=0;EDCD->PASSFAIL[13-1]=0;}

/*EDC 14 - ROO quasi-steady-state*/
if ((EDC==1 || DIAG==1) & meangpp*froot/(pars[6])>pars[19]*EQF){EDC=0;EDCD->PASSFAIL[14-1]=0;}
if ((EDC==1 || DIAG==1) & meangpp*froot/(pars[6])<pars[19]/EQF){EDC=0;EDCD->PASSFAIL[14-1]=0;}





/***********************EDCs done here****************************/



/*Additional faults can be stored in positions 35-40*/

/*PRIOR RANGES - ALL POOLS MUST CONFORM*/

for (n=0;n<6;n++){if ((EDC==1 || DIAG==1) & ((MPOOLS[n])>PI.parmax[n+17])){EDC=0;EDCD->PASSFAIL[35-1]=0;}}

int PEDC;
/*ensuring minimum of each pool is zero & finite*/
if (EDC==1 || DIAG==1)
{double min; int nn;n=0;
while (n<nopools & (EDC==1 || DIAG==1))
{nn=0;PEDC=1;while (nn<nodays+1 & (PEDC==1))
{if (POOLS[n+nn*nopools]<0 || isnan(POOLS[36-1])==1)
{EDC=0;PEDC=0;EDCD->PASSFAIL[35+n]=0;}nn=nn+1;};
n=n+1;}}









/*FREE MEMORY*/
free(MPOOLS);





/*final check confirming EDC = 1 or 0*/
int Num_EDC=100;
if (DIAG==1){for (n=0;n<Num_EDC;n++){if (EDCD->PASSFAIL[n]==0){EDC=0;}}}

/*Returning EDC */
return EDC;












}
























