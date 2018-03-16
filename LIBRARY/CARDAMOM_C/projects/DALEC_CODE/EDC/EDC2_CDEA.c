#pragma once
#include "../../../math_fun/std.c"




/*Component of expdecay (next function)*/
double expdecay2(double const *POOLS, int pool, int nodays,int deltat,int nopools)
/*only accepting pool and number of days, the rest done here*/
{
/*using 365 day averaging window for each year!*/

/*explicitly calculating number of years*/
int noyears=floor(nodays*deltat/365);
int EDC=1,n,count,y;

/*yearly means and yearly means with offset Ds = 100 days*/
double P0=POOLS[pool];
/*Deriving exponential decay coefficients a,b and c in equation
 * Cexp = a + b*exp(c*t)*/
/*four mean pool values to be derived*/
/* MP0 = mean pool (year 1 to year end-2)
 * MP1 = mean pool (year 2 to year end-1)
 * MP0os = mean pool (year 1+os to year end-2+os)
 * MP1os = mean pool (year 1+os to year end-2+os)*/
double MP0=0,MP1=0,MP0os=0,MP1os=0;
/*averaging window*/
int aw=floor(365./(double)deltat);
/*offset in days*/
/*OFFSET = 1 is ideal to capture rapid exponential decay without compromising quality of fit*/
int os=1;
/*deriving coefficients to *
 * in Cexp = A + B exp(C*t);*/
double A,B,b,C;

/*mean pool within each averaging window*/
for (n=0;n<aw;n++){MP0=MP0+POOLS[n*nopools+pool];};MP0=MP0/(double)aw;
for (n=aw;n<aw*2;n++){MP1=MP1+POOLS[n*nopools+pool];};MP1=MP1/(double)aw;
for (n=os;n<aw+os;n++){MP0os=MP0os+POOLS[n*nopools+pool];};MP0os=MP0os/(double)aw;
for (n=aw+os;n<aw*2+os;n++){MP1os=MP1os+POOLS[n*nopools+pool];};MP1os=MP1os/(double)aw;

/*deriving mean gradient ratio dcdt1/dcdt0*/
/*dcdt1 is the numeric gradient between n+1 and n+365+1*/
/*dcdt0 is the numeric gradient between n and n+365*/
double dcdtr=0,dcdt1,dcdt0;

double dcdty1=(MP1os-MP0os);
double dcdty0=(MP1-MP0);

/*denominators*/

/*using multiyear mean to determine c*/
C=log(dcdty1/dcdty0)/((double)os*deltat);
/*deriving final exp. decay fit with startpoint = startpoint of pool*/

/*compared and validated against dalec_expdecay3.m*/
/*
printf("noyears = %d\n",noyears);
printf("MP0 = %f\n",MP0);
printf("MP1 = %f\n",MP1);
printf("MP0os = %f\n",MP0os);
printf("MP1os = %f\n",MP1os);
printf("A = %e\n",A);
printf("B = %e\n",B);
printf("C = %e\n",C);
printf("aw = %d\n",aw);
printf("dcdty1 = %f\n",dcdt1);
printf("dcdty0 = %f\n",dcdt0);
*/

/*half life must be more than noyears*/

/*if (fabs(-log(2)/C)<noyears*365.25 & C<0 & finite(C) & abs(B)>abs(A)*0.01){EDC=0;}*/


return C;


}



/*mean matrix from double pointer routine*/
double mean_annual_pool(double *POOLS, int year, int pool, int nopools,int deltat){
/*inputs
 * POOLS: Pools double pointer, as output from DALEC
 * year: year for which to average (first year = 0)
 * pool: the specific pool 
 * nc
/*declarations*/
int r,c;
double meanpool=0;
/*deriving mean of pool p*/
int stday=floor(365.25*year/deltat);
int enday=floor(365.25*(year+1)/deltat);
for (c=stday;c<enday;c++){
meanpool=meanpool+POOLS[c*nopools+pool]/(enday-stday);}
/*returing meanpool value*/
return meanpool;}


/*mean matrix from double pointer routine*/
double mean_pool(double *PA,int p,int nc, int nopools){
/*declarations*/
int r,c;
double meanpool=0;
/*deriving mean of pool p*/
for (c=0;c<nc;c++){
meanpool=meanpool+PA[c*nopools+p]/(double)nc;
}
/*returing meanpool value*/
return meanpool;}


int EDC2_CDEA(double *pars, double *MET, double *LAI, double *NEE, double *POOLS,double *FLUXES,int nodays,PARAMETER_INFO PI,struct EDCDIAGNOSTIC *EDCD,double meantemp)
{



/*THESE EDC2 checks are for DALEC_CDEA*/
int EDC=1,n=0,edc=0;
int DIAG=EDCD->DIAG;/*1 or 0*/


/*CDEA*/
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



/*EDC LAI related:*/
/*now obsolete*/
/*MAX(Cfol)/clma <LAImax*/
/*double LAImax=10;
if (EDC==1 || DIAG==1){
for (n=0;n<nodays;n++){
if (LAI[n]>LAImax){EDC=0;EDCD->PASSFAIL[5-1]=0;}}}*/
/*MLAI=mean_pointer_vector(LAI,nodays);*/



/*EDC no 6*/
/*0.2*Cf < Cr < 5*Cf*/
/*Cfoliar : Croot = 5:1 or 1:5*/
if ((EDC==1 || DIAG==1) & (MPOOLS[1]>MPOOLS[2]*5 | MPOOLS[1]*5<MPOOLS[2])){EDC=0;EDCD->PASSFAIL[6-1]=0; }


/*EDC 7*/
/*G10 = 2 - growth factor in 10 years*/
/*for any further constraint please incorporate in likelyhood*/
int noyears=floor(nodays*deltat/365);
double G=0.1;
int y;


/*mean annual pool array*/
double *MAPOOL=malloc(noyears*sizeof(double));
for (n=0;n<nopools;n++){
/*Rapid POOL Growth: cannot increase by more than factor of Gn over N yrs*/
/*note - this is GROWTH ONLY!!!*/

/*step 1 - deriving mean annual pools*/
for (y=0;y<noyears;y++){MAPOOL[y]=mean_annual_pool(POOLS,y,n,nopools,deltat);}

/*step 2 - checking pool growth*/
/*if ((EDC==1 || DIAG==1) & (POOLS[n + nopools*nodays]/POOLS[n]>G*(double)noyears)){EDC=0;EDCD->PASSFAIL[7-1]=0;}}*/
if ((EDC==1 || DIAG==1) & (MAPOOL[noyears-1]/MAPOOL[0]>(1+G*(double)noyears))){EDC=0;EDCD->PASSFAIL[7-1]=0;}

}




/*EDC no 8*/
/*POOL EXPONENTIAL DECAY*/
/*Performed on all seven pools*/

/*exporting the decay coefficient here*/
double C;
for (n=0;n<nopools;n++){if (EDC==1 || DIAG==1) {
C=expdecay2(POOLS,n,nodays+1,deltat,nopools);
if (fabs(-log(2)/C)<365.25*noyears & C<0 & finite(C)){EDC=0;}

if (EDC==0){EDCD->PASSFAIL[8-1]=0;}}};




double const fauto=pars[1];
double const ffol=(1-fauto)*pars[2];
double const flab=(1-fauto-ffol)*pars[12];
double const froot=(1-fauto-ffol-flab)*pars[3];
double const fwood=1-fauto-ffol-flab-froot;

/*fraction of GPP to som and lit under equilibrium conditions*/
double const fsom=fwood+(froot+flab+ffol)*pars[0]/(pars[0]+pars[7]);
double const flit=(froot+flab+ffol);



/*SOM attractor - must be within a factor of 2 from Csom0*/
/*half the TR is balanced by x2 Csom*/

/*equilibrium factor (in comparison to C_initial)*/
double EQF=10.0;
double meangpp=0;
for (n=0;n<nodays;n++){meangpp+=FLUXES[n*nofluxes]/(double)nodays;}




/*EDC 9 - SOM steady-state attractor proximity attractor proximity*/
if ((EDC==1 || DIAG==1) & meangpp*fsom/(pars[8]*exp(pars[9]*meantemp))>pars[22]*EQF){EDC=0;EDCD->PASSFAIL[9-1]=0;}
if ((EDC==1 || DIAG==1) & meangpp*fsom/(pars[8]*exp(pars[9]*meantemp))<pars[22]/EQF){EDC=0;EDCD->PASSFAIL[9-1]=0;}

/*EDC 10 - LIT  -steady-state attractor proximity*/
if ((EDC==1 || DIAG==1) & meangpp*flit/(pars[7]*exp(pars[9]*meantemp))>pars[21]*EQF){EDC=0;EDCD->PASSFAIL[10-1]=0;}
if ((EDC==1 || DIAG==1) & meangpp*flit/(pars[7]*exp(pars[9]*meantemp))<pars[21]/EQF){EDC=0;EDCD->PASSFAIL[10-1]=0;}

/*EDC 11 - WOO  -steady-state attractor proximity*/
if ((EDC==1 || DIAG==1) & meangpp*fwood/(pars[5])>pars[20]*EQF){EDC=0;EDCD->PASSFAIL[11-1]=0;}
if ((EDC==1 || DIAG==1) & meangpp*fwood/(pars[5])<pars[20]/EQF){EDC=0;EDCD->PASSFAIL[11-1]=0;}

/*EDC 12 - ROO  -steady-state attractor proximity*/
if ((EDC==1 || DIAG==1) & meangpp*froot/(pars[6])>pars[19]*EQF){EDC=0;EDCD->PASSFAIL[12-1]=0;}
if ((EDC==1 || DIAG==1) & meangpp*froot/(pars[6])<pars[19]/EQF){EDC=0;EDCD->PASSFAIL[12-1]=0;}





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
free(MAPOOL);





/*final check confirming EDC = 1 or 0*/
int Num_EDC=100;
if (DIAG==1){for (n=0;n<Num_EDC;n++){if (EDCD->PASSFAIL[n]==0){EDC=0;}}}

/*Returning EDC */
return EDC;












}
























