#include "../../../../auxi_fun/filediag.c"
#include "../../../../auxi_fun/oksofar.c"

/*USER DEFINED STRUCTURE*/
/*This is the ONLY structure the user should alter, adapt, etc.,
 *  * in order to optimally stored all required model drivers, and measurements
 *   * observations, etc. etc.*/
typedef struct{
/*DRIVERS*/
double *MET;
double meantemp;
double meanrad;
/*OBS*/
double *GPP;
double *NEE;
double *LAI;
double *WOO;
int *gpppts;
int *neepts;
int *woopts;
int *laipts;
int ngpp;
int nnee;
int nlai;
int nwoo;
/*saving computational speed by allocating memory to model output*/
double *M_GPP;
double *M_NEE;
double *M_LAI;
double *M_FLUXES;
double *M_POOLS;
double *C_POOLS;
/*static data*/
int nodays;
double deltat;
double LAT;
int ID;
int noobs;
int nomet;
int nopools;
int RC;
/*binary file mcmc options (need to add all options HERE except inout files)*/
int rc_random_search;

/*priors*/
double parpriors[50];
double parpriorunc[50];
double otherpriors[50];
double otherpriorunc[50];
}DATA;





int READ_FIRES_SITELEVEL_DATA(char *filename,DATA *DATA)
{

/*NOTE: this function reads data as written by DALEC_FLUXCOM_MCMC_PROJECT_SITELEVEL.m
 * For any adaptations to this function make sure to keep in sync with matlab function*/

/*TEMPLATE FOR ALL DALEC MCMC DATA files*/
/*Static Elements: 1-100 - use as many as needed*/
/*Parameter Priors: 101-150*/
/*Parameter prior uncertainty: 151-200*/
/*Other priors & uncertainties: 201-300*/
/*TEMPORAL DRIVERS & DATA: 301-end*/




FILE *fid=fopen(filename,"rb");
filediag(fid,filename);

/*READING STATIC DATA*/

double statdat[100];
/*reading 1-100*/
fread(statdat,sizeof(double),100,fid);
/*DALEC model run info*/
DATA->ID=(int)statdat[0];
DATA->LAT=statdat[1];
DATA->nodays=(int)statdat[2];
DATA->nomet=(int)statdat[3];
DATA->noobs=(int)statdat[4];
DATA->RC=(int)statdat[5];

/*DALEC MCMC run info*/
/*set to 1 for (a) few and (b) well constrained priors, otherwise 0*/
DATA->rc_random_search=(int)statdat[10];



/*UP TO USER to read data and allocate it to DATA structure*/
double parpriors[50],parpriorunc[50],otherpriors[50],otherpriorunc[50];
fread(parpriors,sizeof(double),50,fid);
fread(parpriorunc,sizeof(double),50,fid);
fread(otherpriors,sizeof(double),50,fid);
fread(otherpriorunc,sizeof(double),50,fid);


/*For universal data structure, DATA contains 50 parameter prior spaces
 * use as many as needed!*/
memcpy(DATA->parpriors,parpriors,50*sizeof(double));
memcpy(DATA->parpriorunc,parpriorunc,50*sizeof(double));
memcpy(DATA->otherpriors,otherpriors,50*sizeof(double));
memcpy(DATA->otherpriorunc,otherpriorunc,50*sizeof(double));



/*READING TEMPORAL DATA*/

if (DATA->MET==0){DATA->MET=calloc(DATA->nomet*DATA->nodays,sizeof(double));}
if (DATA->GPP==0){DATA->GPP=calloc(DATA->nodays,sizeof(double));}
if (DATA->NEE==0){DATA->NEE=calloc(DATA->nodays,sizeof(double));}
if (DATA->LAI==0){DATA->LAI=calloc(DATA->nodays,sizeof(double));}
if (DATA->WOO==0){DATA->WOO=calloc(DATA->nodays,sizeof(double));}

DATA->ngpp=0;
DATA->nlai=0;
DATA->nnee=0;
DATA->nwoo=0;



/*6 met fields for DALEC FIRES*/

int n,nn;
double *metline, *obsline;
metline=calloc(DATA->nomet,sizeof(double));
obsline=calloc(DATA->noobs,sizeof(double));


for (n=0;n<DATA->nodays;n++){
	fread(metline,sizeof(double),DATA->nomet,fid);
	fread(obsline,sizeof(double),DATA->noobs,fid);
	for (nn=0;nn<DATA->nomet;nn++){DATA->MET[n*DATA->nomet+nn]=metline[nn];}
	DATA->GPP[n]=obsline[0];
	DATA->LAI[n]=obsline[1];
	DATA->NEE[n]=obsline[2];
	if (obsline[0]>-9998){DATA->ngpp=DATA->ngpp+1;}
	if (obsline[1]>-9998){DATA->nlai=DATA->nlai+1;}
	if (obsline[2]>-9998){DATA->nnee=DATA->nnee+1;}
	if (DATA->noobs>3){DATA->WOO[n]=obsline[3];if (obsline[3]>-9998){DATA->nwoo=DATA->nwoo+1;}}
}
DATA->deltat=DATA->MET[DATA->nomet]-DATA->MET[0];



fclose(fid);

if (DATA->ngpp>0){DATA->gpppts=calloc(DATA->ngpp,sizeof(int));}
if (DATA->nlai>0){DATA->laipts=calloc(DATA->nlai,sizeof(int));}
if (DATA->nnee>0){DATA->neepts=calloc(DATA->nnee,sizeof(int));}
if (DATA->nwoo>0){DATA->woopts=calloc(DATA->nwoo,sizeof(int));}

/*Deriving laipts, gpppts, neepts*/
int c;
c=0;for (n=0;n<DATA->nodays;n++){if (DATA->GPP[n]>-9998){DATA->gpppts[c]=n;c=c+1;}}
c=0;for (n=0;n<DATA->nodays;n++){if (DATA->LAI[n]>-9998){DATA->laipts[c]=n;c=c+1;}}
c=0;for (n=0;n<DATA->nodays;n++){if (DATA->NEE[n]>-9998){DATA->neepts[c]=n;c=c+1;}}
if (DATA->noobs>3){c=0;for (n=0;n<DATA->nodays;n++){if (DATA->WOO[n]>-9998){DATA->woopts[c]=n;c=c+1;}}}


/*deriving mean temp and mean rad*/
DATA->meantemp=0;
DATA->meanrad=0;

printf("Mean Rad = %f\n",DATA->meanrad);
printf("Mean Temp = %f\n",DATA->meantemp);
printf("No days = %d\n",DATA->nodays);
for (n=0;n<DATA->nodays;n++){DATA->meantemp+=0.5*DATA->MET[DATA->nomet*n+1]/(double)DATA->nodays;}
for (n=0;n<DATA->nodays;n++){DATA->meantemp+=0.5*DATA->MET[DATA->nomet*n+2]/(double)DATA->nodays;}
for (n=0;n<DATA->nodays;n++){DATA->meanrad+=DATA->MET[DATA->nomet*n+3]/(double)DATA->nodays;}

printf("Mean Rad = %f\n",DATA->meanrad);
printf("Mean Temp = %f\n",DATA->meantemp);



free(metline);
free(obsline);

return 0;



}



int INITIALIZE_DATA_STRUCT(DATA *DATA){

/*initialising array pointers as zero*/
DATA->MET=0;
DATA->GPP=0;
DATA->LAI=0;
DATA->NEE=0;
DATA->WOO=0;
return 0;


}


int FREE_DATA_STRUCT(DATA DATA){


if (DATA.ngpp>0){free(DATA.gpppts);}
if (DATA.nlai>0){free(DATA.laipts);}
if (DATA.nnee>0){free(DATA.neepts);}
if (DATA.nwoo>0){free(DATA.woopts);}

free(DATA.MET);
free(DATA.LAI);
free(DATA.NEE);
free(DATA.WOO);
free(DATA.GPP);



return 0;

}



