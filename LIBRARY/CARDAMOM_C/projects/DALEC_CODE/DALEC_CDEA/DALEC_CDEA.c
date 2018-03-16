#pragma once
#include "ACM.c"

/*

DALEC - All Biomes: Crops, Deciduous, Evergreen,Forests
This routine is C only 


*/




/*
/* GPP RESPONSE CURVE HERE
double gppresponse(double const *pars, double const *consts)
{
  /*double gc=0,pp=0,qq=0,ci=0,e0=0,mult=0,dayl=0,cps=0,dec=0;
  double gc,pp,qq,ci,e0,mult,dayl,cps,dec,GPP;
  /*pars= &pars;  
  consts= &consts;

  gc=(double)pow(fabs(pars[8]),consts[9])/(consts[5] * pars[9] + 0.5 * ( pars[1]- pars[2]));
  pp=(double)pars[0]*pars[3]/gc*consts[0]*exp(consts[7]*pars[1]);
  qq=(double)consts[2]-consts[3];
  ci=(double)0.5*(pars[4]+qq-pp+pow(pow(pars[4]+qq-pp,2)-4*(pars[4]*qq-pp*consts[2]),0.5));
  e0=(double)consts[6]*pow(pars[0],2)/(pow(pars[0],2)+consts[8]);
  dec=(double)-23.4*cos((360.*(pars[5]+10.)/365.)*pars[10]/180.)*pars[10]/180.;
  mult=(double)tan(pars[6]*pars[10]/180)*tan(dec);
  if (mult>=1){ 
   dayl=0.;}  
  else if(mult<=-1)
  dayl=24.;
  else{ 
  dayl=(double)24.*acos(-mult) / pars[10];}
  

  cps=(double)e0*pars[7]*gc*(pars[4]-ci)/(e0*pars[7]+gc*(pars[4]-ci));
  GPP=cps*(consts[1]*dayl+consts[4]);
  return GPP;
}

*/


double offset(double const L, double const w)
{

/*solution to t = f*exp(-t^2)*/
/*see dalecstepfunction2*/

double mxc[7]={0.000023599784710, 0.000332730053021,    0.000901865258885,  -0.005437736864888,  -0.020836027517787,   0.126972018064287,   -0.188459767342504};

double lf=log(L-1);
double os=mxc[0]*pow(lf,6) + mxc[1]*pow(lf,5) + mxc[2]*pow(lf,4) + mxc[3]*pow(lf,3) + mxc[4]*pow(lf,2) + mxc[5]*lf +mxc[6];

os=os*w;

return os;
}

void DALEC_CDEA(double const *MET,double const *pars, double const deltat,int const nr, double const lat,
double *LAI, double *NEE, double *FLUXES, double *POOLS)
{

double gpppars[11],pi;
/*C-pools, fluxes, meteorology indices*/
int p,f,m,nxp, i;
int n=0;
pi=3.1415927;


/*constant gpppars terms*/
gpppars[3]=1;
gpppars[6]=lat;
gpppars[8]=-2.0;
gpppars[9]=1.0;
gpppars[10]=pi;


 double constants[10]={pars[10],0.0156935,4.22273,208.868,0.0453194,0.37836,7.19298, 0.011136,2.1001,0.789798};

  /*assigning values to pools*/
  /*L,F,R,W,Lit,SOM*/
  POOLS[0]=pars[17];
  POOLS[1]=pars[18];
  POOLS[2]=pars[19];
  POOLS[3]=pars[20];
  POOLS[4]=pars[21];
  POOLS[5]=pars[22];



/* NOTES FOR POOLS AND FLUXES
MET[:,0]: projday
MET[:,1]: mintemp
MET[:,2]: maxtemp
MET[:,3]: rad
MET[:,4]: co2
MET[:,5]: yearday



  POOLS[0,0]=pars(8);L
  POOLS[0,1]=pars(5);F
  POOLS[0,2]=pars(6);R
  POOLS[0,3]=pars(3);W
  POOLS[0,4]=pars(2);Litter
  POOLS[0,5]=pars(2);Som


        %fluxes - other*********
        0.GPP
        1.temprate
        2.respiration_auto
        3.leaf_production
        4.labile_production
        5.root_production
        6.wood_production
        7.labile_release
        8.leaffall_factor
        9.leaflitter_production
        10.woodlitter_production  
        11.rootlitter_production         
     	12.respiration_het_litter
  	13.respiration_het_som
  	14.litter2som
  	15.labrelease_factor
*/



/*constants for exponents of leaffall and labrelease factors*/
/*width*/
double wf=pars[15]*sqrt(2)/2;
double wl=pars[13]*sqrt(2)/2;


/*factor*/
double ff=(log(pars[4])-log(pars[4]-1))/2;
double fl=(log(1.001)-log(0.001))/2;




/*additional offset*/
double osf=offset(pars[4],wf);
double osl=offset(1.001,wl);


/*scaling to biyearly sine curve*/
double sf=365.25/pi;







/*repeating loop for each timestep*/
for ( n=0; n < nr; n++){
/*LAI*/
p=6*n;nxp=6*(n+1);m=6*n;f=16*n;
LAI[n]=POOLS[p+1]/pars[16]; 
/*GPP*/
gpppars[0]=LAI[n];
gpppars[1]=MET[m+2];
gpppars[2]=MET[m+1];
gpppars[4]=MET[m+4];
gpppars[5]=MET[m+5];
gpppars[7]=MET[m+3];



/*GPP - daily timestep*/
FLUXES[f+0]=ACM(gpppars,constants);
/*temprate - now comparable to Q10 - factor at 0C is 1*/
FLUXES[f+1]=exp(pars[9]*0.5*(MET[m+2]+MET[m+1]));
/*respiration auto*/
FLUXES[f+2]=pars[1]*FLUXES[f+0];
/*leaf production rate*/
FLUXES[f+3]=(FLUXES[f+0]-FLUXES[f+2])*pars[2];
/*labile production*/
FLUXES[f+4] = (FLUXES[f+0]-FLUXES[f+2]-FLUXES[f+3])*pars[13-1];              

/*root production*/        
FLUXES[f+5] = (FLUXES[f+0]-FLUXES[f+2]-FLUXES[f+3]-FLUXES[f+4])*pars[4-1];            

/*wood production rate*/       
FLUXES[f+6] = FLUXES[f+0]-FLUXES[f+2]-FLUXES[f+3]-FLUXES[f+5]-FLUXES[f+4]; 


/*leaf fall factor*/
/*FLUXES[f+8]=pow((0.75+cos((MET[m+5]-pars[15-1])*(2*pi/365.25))/4),10)*pars[5-1];*/
/*FLUXES[f+8]=1- pow((1 +0.00001 - 1/pars[4]),nf*exp(nf*sin((MET[m+0]-pars[14])/sf)*sf)/pow(1+exp(nf*sin((MET[m+0]-pars[14])/sf)*sf),2));*/
FLUXES[f+8] = (2/sqrt(pi))*(ff/wf)*exp(-pow(sin((MET[m+0]-pars[14]+osf)/sf)*sf/wf,2));


/*Labrelease factor*/
FLUXES[f+15]=(2/sqrt(pi))*(fl/wl)*exp(-pow(sin((MET[m+0]-pars[11]+osl)/sf)*sf/wl,2));
 

/*total labile release - mean at daily timestep */
FLUXES[f+7]=POOLS[p+0]*(1-pow(1-FLUXES[f+15],deltat))/deltat;
     
/*total leaf litter production*/       
FLUXES[f+9] = POOLS[p+1]*(1-pow(1-FLUXES[f+8],deltat))/deltat;                     
 
/*total wood litter production*/       
FLUXES[f+10] = POOLS[p+3]*(1-pow(1-pars[6-1],deltat))/deltat;                    

/*root litter production*/
FLUXES[f+11] = POOLS[p+2]*(1-pow(1-pars[7-1],deltat))/deltat;                    
        

/*respiration heterotrophic litter*/
FLUXES[f+12] = POOLS[p+4]*(1-pow(1-FLUXES[f+1]*pars[8-1],deltat))/deltat;
        
/*respiration heterotrophic SOM*/
FLUXES[f+13] = POOLS[p+5]*(1-pow(1-FLUXES[f+1]*pars[9-1],deltat))/deltat;         

/*litter to SOM*/
FLUXES[f+14] = POOLS[p+4]*(1-pow(1-pars[1-1]*FLUXES[f+1],deltat))/deltat; 


/*all fluxes are at a daily timestep*/
NEE[n]=(-FLUXES[f+0]+FLUXES[f+2]+FLUXES[f+12]+FLUXES[f+13]);

        POOLS[nxp+0] = POOLS[p+0] + (FLUXES[f+4]-FLUXES[f+7])*deltat;
        POOLS[nxp+1] =  POOLS[p+1] + (FLUXES[f+3] - FLUXES[f+9] + FLUXES[f+7])*deltat;
        POOLS[nxp+3] = POOLS[p+3] +  (FLUXES[f+6] - FLUXES[f+10])*deltat;
        POOLS[nxp+2] = POOLS[p+2] + (FLUXES[f+5] - FLUXES[f+11])*deltat;
        POOLS[nxp+4] = POOLS[p+4] + (FLUXES[f+9] + FLUXES[f+11] - FLUXES[f+12] - FLUXES[f+14])*deltat;        
        POOLS[nxp+5]= POOLS[p+5]+ (FLUXES[f+14] - FLUXES[f+13]+FLUXES[f+10])*deltat;                                



}
}







