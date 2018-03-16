#pragma once
#include "ACM.c"

/*

DALEC - All Biomes: Crops, Deciduous, Evergreen,Forests
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

void DALEC_FIRES(double const *MET,double const *pars, double const deltat,int const nr, double const lat,
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

/*fully combusted fire flux (all pools except SOM)*/
double CFF[5];

/*fully combusted fire flux (all pools except SOM)*/
double NCFF[5];

/*combustion efficiencies*/
double combust_eff[5]={0.1,0.9,0.1,0.5,0.3};

/*resilience factor*/
double rfac=0.5;


/*number of MET drivers*/
int nomet=7;

/*number of DALEC pools*/
int nopools=6;

/*number of DALEC fluxes to store*/
int nofluxes=16;


/*repeating loop for each timestep*/
for ( n=0; n < nr; n++){
/*ppol index*/
p=nopools*n;
/*next pool index*/
nxp=nopools*(n+1);
/*met index*/
m=nomet*n;
/*flux array index*/
f=nofluxes*n;



/*LAI*/
LAI[n]=POOLS[p+1]/pars[16]; 
/*GPP*/
gpppars[0]=LAI[n];
gpppars[1]=MET[m+2];
gpppars[2]=MET[m+1];
gpppars[4]=MET[m+4];
gpppars[5]=MET[m+5];
gpppars[7]=MET[m+3];




FLUXES[f+0]=ACM(gpppars,constants);
/*temprate - now comparable to Q10 - factor at 0C is 1*/
FLUXES[f+1]=exp(pars[9]*0.5*(MET[m+2]+MET[m+1]));
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





/*total pool transfers (no fires yet)*/

        POOLS[nxp+0] = POOLS[p+0] + (FLUXES[f+4]-FLUXES[f+7])*deltat;
        POOLS[nxp+1] =  POOLS[p+1] + (FLUXES[f+3] - FLUXES[f+9] + FLUXES[f+7])*deltat;
        POOLS[nxp+3] = POOLS[p+3] +  (FLUXES[f+6] - FLUXES[f+10])*deltat;
        POOLS[nxp+2] = POOLS[p+2] + (FLUXES[f+5] - FLUXES[f+11])*deltat;
        POOLS[nxp+4] = POOLS[p+4] + (FLUXES[f+9] + FLUXES[f+11] - FLUXES[f+12] - FLUXES[f+14])*deltat;        
        POOLS[nxp+5]= POOLS[p+5]+ (FLUXES[f+14] - FLUXES[f+13]+FLUXES[f+10])*deltat;                                

/*total pool transfers - WITH FIRES*/
/*first fluxes*/

	/*LABILE*/
	CFF[0] = POOLS[nxp+0]*MET[m+6]*combust_eff[0];
	NCFF[0] = POOLS[nxp+0]*MET[m+6]*(1-combust_eff[0])*(1-rfac);

	/*foliar*/
	CFF[1] = POOLS[nxp+1]*MET[m+6]*combust_eff[1];
	NCFF[1] = POOLS[nxp+1]*MET[m+6]*(1-combust_eff[1])*(1-rfac);

	/*root*/
	CFF[2] = POOLS[nxp+2]*MET[m+6]*combust_eff[2];
	NCFF[2] = POOLS[nxp+2]*MET[m+6]*(1-combust_eff[2])*(1-rfac);

	/*wood*/
	CFF[3] = POOLS[nxp+3]*MET[m+6]*combust_eff[3];
	NCFF[3] = POOLS[nxp+3]*MET[m+6]*(1-combust_eff[3])*(1-rfac);

	/*litter*/
	CFF[4] = POOLS[nxp+4]*MET[m+6]*combust_eff[4];
	NCFF[4] = POOLS[nxp+4]*MET[m+6]*(1-combust_eff[4])*(1-rfac);


/*Adding all fire pool transfers here*/
	POOLS[nxp+0]=POOLS[nxp+0]-CFF[0]-NCFF[0];
	POOLS[nxp+1]=POOLS[nxp+1]-CFF[1]-NCFF[1];
	POOLS[nxp+2]=POOLS[nxp+2]-CFF[2]-NCFF[2];
	POOLS[nxp+3]=POOLS[nxp+3]-CFF[3]-NCFF[3];
	POOLS[nxp+4]=POOLS[nxp+4]-CFF[4]-NCFF[4]+NCFF[0]+NCFF[1]+NCFF[2];
	POOLS[nxp+5]=POOLS[nxp+5]+NCFF[3]+NCFF[4];


NEE[n]=-FLUXES[f+0]+FLUXES[f+2]+FLUXES[f+12]+FLUXES[f+13]+(CFF[0]+CFF[1]+CFF[2]+CFF[3]+CFF[4])/deltat;


}
}







