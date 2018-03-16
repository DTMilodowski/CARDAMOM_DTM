#pragma once
#include "ACM_BUCKET.c"
/*

DALEC - All Bucket: Crops, Deciduous, Evergreen,Forests
This routine is C only 
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


/*main function*/

void DALEC_BUCKET(double const *MET,double const *pars, double const deltat,int const nr, double const lat,
double *LAI, double *NEE, double *FLUXES, double *POOLS)
{

double gpppars[13],pi;
/*C-pools, fluxes, meteorology indices*/
int p,f,m,nxp, i;
int n=0;
pi=3.1415927;


/*constant gpppars terms*/
gpppars[3]=1;
gpppars[6]=lat;
/*defining pi*/
gpppars[10]=pi;



 double constants[11]={pars[10],0.0156935,4.22273,208.868,0.0453194,0.37836,7.19298, 0.011136,2.1001,0.789798,pars[23]};



/*
% p(1) Litter to SOM conversion rate  - m_r
% p(2) Fraction of GPP respired - f_a
% p(3) Fraction of NPP alloco2ted to foliage - f_f 
% p(4) Fraction of NPP2 alloco2ted to roots - f_r
% p(5) Leaf lifespan - L_f
% p(6) Turnover rate of wood - t_w
% p(7) Turnover rate of roots - t_r
% p(8) Litter turnover rate - t_l
% p(9) SOM turnover rate  - t_S
% p(10) Parameter in exponential term of temperature - \theta
% p(11) Canopy efficiency parameter - C_eff
% p(12) = date of Clab release - B_day  
% p(13) = Fraction allocated to Clab - f_l
% p(16) = leaf fall duration period - R_f
% p(15) = date of leaf fall - F_day
% p(14) = lab release duration period - R_l
% p(17) = LMA
% p(24) = Water limitation exponent
% p(25) = Water use efficiency
*/



  /*assigning values to pools*/
  /*L,F,R,W,Lit,SOM*/
  POOLS[0]=pars[17];
  POOLS[1]=pars[18];
  POOLS[2]=pars[19];
  POOLS[3]=pars[20];
  POOLS[4]=pars[21];
  POOLS[5]=pars[22];
  POOLS[6]=pars[25];



/* NOTES FOR POOLS AND FLUXES
MET[:,0]: projday
MET[:,1]: mintemp
MET[:,2]: maxtemp
MET[:,3]: rad
MET[:,4]: co2
MET[:,5]: yearday
MET[:,6]: prec



  POOLS[0,0]=pars(18);L
  POOLS[0,1]=pars(19);F
  POOLS[0,2]=pars(20);R
  POOLS[0,3]=pars(21);W
  POOLS[0,4]=pars(22);Litter
  POOLS[0,5]=pars(23);Som
  POOLS[0,6]=pars(26);Wat


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

int nn;





/*repeating loop for each timestep*/
for ( n=0; n < nr; n++){
/*LAI*/
p=7*n;nxp=7*(n+1);m=8*n;f=16*n;
LAI[n]=POOLS[p+1]/pars[16]; 
/*GPP*/
gpppars[0]=LAI[n];
gpppars[1]=MET[m+2];
gpppars[2]=MET[m+1];
gpppars[4]=MET[m+4];
gpppars[5]=MET[m+5];
gpppars[7]=MET[m+3];
gpppars[11]=POOLS[p+6];



/*psi rtot*/
gpppars[8]=-4.0;
gpppars[9]=1.0;




FLUXES[f+0]=ACM_BUCKET(gpppars,constants);
/*temprate*/
FLUXES[f+1]=0.5*exp(pars[9]*0.5*(MET[m+2]+MET[m+1]));
/*respiration auto*/
FLUXES[f+2]=pars[1]*FLUXES[f+0];
/*leaf production*/
FLUXES[f+3]=(FLUXES[f+0]-FLUXES[f+2])*pars[2];
/*labile production*/
FLUXES[f+4] = (FLUXES[f+0]-FLUXES[f+2]-FLUXES[f+3])*pars[13-1];              

/*root production*/        
FLUXES[f+5] = (FLUXES[f+0]-FLUXES[f+2]-FLUXES[f+3]-FLUXES[f+4])*pars[4-1];            

/*wood production*/       
FLUXES[f+6] = FLUXES[f+0]-FLUXES[f+2]-FLUXES[f+3]-FLUXES[f+5]-FLUXES[f+4]; 


/*leaf fall factor*/
/*FLUXES[f+8]=pow((0.75+cos((MET[m+5]-pars[15-1])*(2*pi/365.25))/4),10)*pars[5-1];*/
/*FLUXES[f+8]=1- pow((1 +0.00001 - 1/pars[4]),nf*exp(nf*sin((MET[m+0]-pars[14])/sf)*sf)/pow(1+exp(nf*sin((MET[m+0]-pars[14])/sf)*sf),2));*/
FLUXES[f+8] = (2/sqrt(pi))*(ff/wf)*exp(-pow(sin((MET[m+0]-pars[14]+osf)/sf)*sf/wf,2));


/*Labrelease factor*/
FLUXES[f+15]=(2/sqrt(pi))*(fl/wl)*exp(-pow(sin((MET[m+0]-pars[11]+osl)/sf)*sf/wl,2));
 

/*labile release - FIXFIXFIX*/
FLUXES[f+7]=POOLS[p+0]*FLUXES[f+15];
     
/*leaf litter production*/       
FLUXES[f+9] = POOLS[p+1]*FLUXES[f+8];                     
 
/*wood litter production*/       
FLUXES[f+10] = pars[6-1]*POOLS[p+3];                    

/*root litter production*/
FLUXES[f+11] = pars[7-1]*POOLS[p+2];                    
/*respiration heterotrophic litter*/
FLUXES[f+12] = pars[8-1]*POOLS[p+4]*FLUXES[f+1];
/*respiration heterotrophic SOM*/
FLUXES[f+13] = pars[9-1]*POOLS[p+5]*FLUXES[f+1];         
/*litter to SOM*/
FLUXES[f+14] = pars[1-1]*POOLS[p+4]*FLUXES[f+1]; 



NEE[n]=-FLUXES[f+0]+FLUXES[f+2]+FLUXES[f+12]+FLUXES[f+13];

        POOLS[nxp+0] = POOLS[p+0] + (FLUXES[f+4]-FLUXES[f+7])*deltat;
        POOLS[nxp+1] =  POOLS[p+1] + (FLUXES[f+3] - FLUXES[f+9] + FLUXES[f+7])*deltat;
        POOLS[nxp+3] = POOLS[p+3] +  (FLUXES[f+6] - FLUXES[f+10])*deltat;
        POOLS[nxp+2] = POOLS[p+2] + (FLUXES[f+5] - FLUXES[f+11])*deltat;
        POOLS[nxp+4] = POOLS[p+4] + (FLUXES[f+9] + FLUXES[f+11] - FLUXES[f+12] - FLUXES[f+14])*deltat;    
        POOLS[nxp+5]= POOLS[p+5] + (FLUXES[f+14] - FLUXES[f+13]+FLUXES[f+10])*deltat;                                
        /*POOLS[nxp+6]= POOLS[p+6]*1 - (FLUXES[f+0]*MET[m+7]/pars[24])*deltat + MET[m+6]*deltat;                                */

        POOLS[nxp+6]= POOLS[p+6]*(1-pars[24]*POOLS[p+6])  + MET[m+6]*deltat;                                




}

/*

       for (nn=0;nn<8;nn++){printf("%6.6f     ",MET[nn]);}printf("\n");
       for (nn=0;nn<26;nn++){printf("%6.6f ",pars[nn]);}printf("\n");
       for (nn=0;nn<1;nn++){printf("%6.6f ",deltat);}printf("\n");
       for (nn=0;nn<1;nn++){printf("%6.6f ",lat);}printf("\n");
       for (nn=0;nn<1;nn++){printf("%d ",nr);}printf("\n");
*/









}







