#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../DALEC_CODE/DALEC_CDEA/DALEC_CDEA.c"
#include "../DALEC_CDEA_TEMPLATE/MCMC_SETUP/PROJECT_FUN/READ_CDEA_SITELEVEL_DATA.c"


typedef struct{
double *pars;
double *prob;
int nsol;
int npar;
}PARS;



int READ_PARS_FILE(char *path,int nf,PARS *PARS){


int nout=0;
char file[200];
char num[8];
sprintf(num,"%05d",nf);
/* format "file" as .../.../.../DG_01234PARS*/
strcpy(file,path);
strcat(file,"DG_");
strcat(file,num);
strcat(file,"PARS");


FILE *fd=fopen(file,"rb");
filediag(fd,file);
/*NEED TO READ PARAMETERS HERE*/

if (fd!=0){
/*identifying number of rows in file*/
fseek(fd,0,SEEK_END);
PARS->nsol=ftell(fd)/(sizeof(double)*(PARS->npar+1));
fseek(fd,0,SEEK_SET);
/*allocating memory for all parameters*/





/*reading all data*/
double b;
int n;
for (n=0;n<PARS->nsol;n++){
fread(PARS->pars+n*(PARS->npar),sizeof(double),PARS->npar,fd);
/*fread(pl,sizeof(double),PARS->npar,fd);*/
/*reading - discarding likelyhood*/
fread(PARS->prob+n,sizeof(double),1,fd);}

fclose(fd);nout=1;
/*otherwise there's no file*/


}




else{nout=0;}


return nout;

}




int main(int NCLI, char *CLI[]){

/*INPUTS
 *CLI[1] = pathdir
 *CLI[2] = number of pixels in global run
 *CLI[3] = maximum number of mcmc solutions
 *
 *OUTPUTS
 *files with all (1)parameter (2)fluxes, (3)pools and (4) any other dynamics, etc.
 *
 */

int npar=23;
int nmet=6;
int npoo=6;


char projectpath[200], resultspath[200],outputpath[200],metpath[200];

/*copying first argument with path*/
strcpy(projectpath, CLI[1]);
strcpy(resultspath, CLI[1]);
strcpy(outputpath, CLI[1]);
strcpy(metpath, CLI[1]);

/*results in RESULTS folder*/
/*outputs in RESULTS_PROCESSED*/

strcat(resultspath,"/RESULTS/");
strcat(outputpath,"/RESULTS_PROCESSED/");
strcat(metpath,"/DATA/");

int n, NPX,np,m,p,npa;
/*input files*/
char metfile[200],parsfile[200],filenum[8];
/*declaring output files*/
char neefile[200],gppfile[200],poofile[200],parfile[200],laifile[200];

/*defining output files*/
strcpy(neefile,outputpath);strcat(neefile,"NEE.bin");
strcpy(laifile,outputpath);strcat(laifile,"LAI.bin");
strcpy(gppfile,outputpath);strcat(gppfile,"GPP.bin");
strcpy(poofile,outputpath);strcat(poofile,"POO.bin");
strcpy(parfile,outputpath);strcat(parfile,"PAR.bin");


/*defining data structure and reading met file to fill it in!*/
DATA D;INITIALIZE_DATA_STRUCT(&D);
sprintf(filenum,"%05d",1);
strcpy(metfile,metpath);strcat(metfile,"DG_");strcat(metfile,filenum);strcat(metfile,".bin");
READ_CDEA_SITELEVEL_DATA(metfile,&D);


double *M_LAI=malloc(D.nodays*sizeof(double));
double *M_GPP=malloc(D.nodays*sizeof(double));
double *M_NEE=malloc(D.nodays*sizeof(double));
double *M_FLUXES=malloc(D.nodays*sizeof(double)*16);
double *M_POOLS=malloc((D.nodays+1)*sizeof(double)*npoo);


/*printing stuff to make sure everything's OK*/

double Pmax,Wtot;




/*number of pixels*/
NPX=atoi(CLI[2]);

double *NEE,*GPP,*POO,*PAR,*LAI;
/*output fields*/
NEE=calloc(D.nodays*NPX,sizeof(double));
GPP=calloc(D.nodays*NPX,sizeof(double));
LAI=calloc(D.nodays*NPX,sizeof(double));
POO=calloc(D.nodays*NPX*npoo,sizeof(double));
PAR=calloc(npar*NPX,sizeof(double));
/*************/

/*defining parameter structure*/
PARS PARS;
PARS.npar=npar;
PARS.pars=malloc(PARS.npar*atoi(CLI[3])*sizeof(double));
PARS.prob=malloc(atoi(CLI[3])*sizeof(double));
/*data structure*/


/*loading file*/
for (n=0;n<NPX;n++){

printf("Doing File %d\n",n+1);

sprintf(filenum,"%05d",n+1);
printf("filenum = %s\n",filenum);





/*loading parameters file*/

strcpy(parsfile,resultspath);strcat(parsfile,"DG_");strcat(parsfile,filenum);strcat(parsfile,"PARS");
int fileisthere=READ_PARS_FILE(resultspath,n+1,&PARS);



/*ONLY CONTINUING if file exists*/
if (fileisthere==1){
	/*loading met file*/
	strcpy(metfile,metpath);strcat(metfile,"DG_");strcat(metfile,filenum);strcat(metfile,".bin");
	READ_CDEA_SITELEVEL_DATA(metfile,&D);
	/*runnning all pars with met data*/
	/*looping through all parameters*/

	/*deriving max likelyhood - using this to determine if it is worth running model*/
	Pmax=PARS.prob[0];for (np=0;np<PARS.nsol;np++){if (PARS.prob[np]>Pmax){Pmax=PARS.prob[np];}}
	Wtot=0;for (np=0;np<PARS.nsol;np++){if (PARS.prob[np]>Pmax-log(10000)){Wtot+=1;}}
	printf("Wtot = %e\n",Wtot);

	for (np=0;np<PARS.nsol;np++){
		/*RUNNING DALEC*/
	
		if (PARS.prob[np]>Pmax-log(10000)){
		/*better to run likelyhood function to ensure everything is OK!*/
		DALEC_CDEA(D.MET, PARS.pars +np*PARS.npar, D.deltat,D.nodays, D.LAT, M_LAI, M_NEE, M_FLUXES, M_POOLS);



		/*************store/write all data here***********/
		/*fluxes*/
		for (m=0;m<D.nodays;m++){
		NEE[n*D.nodays+m]+=M_NEE[m]/Wtot;
		LAI[n*D.nodays+m]+=M_LAI[m]/Wtot;
		GPP[n*D.nodays+m]+=M_FLUXES[m*16]/Wtot;

		/*pools*/
		for (p=0;p<npoo;p++){POO[p*D.nodays*NPX+n*D.nodays+m]+=log(M_POOLS[(m+1)*npoo+p])/Wtot;}}
	

		/*storing mean pars logarithmically*/
		/*Weigh correctly with parbrob*/
		for (npa=0;npa<PARS.npar;npa++){PAR[n*npar+npa]+=log(PARS.pars[PARS.npar*np+npa])/Wtot;}
	}}
      /*************DONE STORING D****************/

	}
/*done with if*/
}


/*calculating probabilities sum*/

/*transforming log averaged quantities*/
for (n=0;n<npoo*D.nodays*NPX;n++){POO[n]=exp(POO[n]);}
for (n=0;n<npar*NPX;n++){PAR[n]=exp(PAR[n]);}
for (n=0;n<D.nodays*NPX;n++){GPP[n]=GPP[n];}
for (n=0;n<D.nodays*NPX;n++){NEE[n]=NEE[n];}
for (n=0;n<D.nodays*NPX;n++){LAI[n]=LAI[n];}



/*output files*/
FILE *neefd=fopen(neefile,"wb");
FILE *laifd=fopen(laifile,"wb");
FILE *gppfd=fopen(gppfile,"wb");
FILE *poofd=fopen(poofile,"wb");
FILE *parfd=fopen(parfile,"wb");
/**********/



/*writing all data*/
printf("D.nodays*NPX = %d\n", D.nodays*NPX);
printf("M_NEE[0] = %f\n",M_NEE[0]);
printf("M_GPP[0] = %f\n",M_GPP[0]);
printf("M_POOLS[0] = %f\n",M_POOLS[0]);
printf("NEE[0] = %f\n",NEE[0]);
printf("GPP[0] = %f\n",GPP[0]);
printf("POO[0] = %f\n",POO[0]);
fwrite(NEE,sizeof(double),D.nodays*NPX,neefd);fclose(neefd);
fwrite(LAI,sizeof(double),D.nodays*NPX,laifd);fclose(laifd);
fwrite(GPP,sizeof(double),D.nodays*NPX,gppfd);fclose(gppfd);
fwrite(POO,sizeof(double),D.nodays*NPX*npoo,poofd);fclose(poofd);
fwrite(PAR,sizeof(double),npar*NPX,parfd);fclose(parfd);


/*free all memory*/
free(PARS.pars);free(PARS.prob);free(NEE);free(GPP);free(POO);free(PAR);free(LAI);


oksofar("about to free data struct");

free(M_LAI);
free(M_GPP);
free(M_NEE);
free(M_FLUXES);
free(M_POOLS);

FREE_DATA_STRUCT(D);


return 0;}




