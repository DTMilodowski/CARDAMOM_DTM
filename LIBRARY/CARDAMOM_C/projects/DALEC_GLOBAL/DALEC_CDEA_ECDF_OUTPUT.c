#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../DALEC_CODE/DALEC_CDEA/DALEC_CDEA.c"
#include "../CARDAMOM_GENERAL/CARDAMOM_READ_BINARY_DATA.c"


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
int nflu=16;


char projectpath[200], resultspath[200],outputpath[200],metpath[200];

/*copying first argument with path*/
strcpy(projectpath, CLI[1]);
strcpy(resultspath, CLI[1]);
strcpy(outputpath, CLI[1]);
strcpy(metpath, CLI[1]);

int maxnp=atoi(CLI[3]);

/*outputs in RESULTS_PROCESSED*/

strcat(resultspath,"/RESULTS/");
strcat(outputpath,"/RESULTS_PROCESSED/");
strcat(metpath,"/DATA/");

int n, NPX,np,m,p,npa;
/*input files*/
char metfile[200],parsfile[200],filenum[8];
/*declaring output files*/
char neefile[200],gppfile[200],poofile[200],parfile[200],laifile[200],nppfile[200],rsofile[200];

/*declaring mean output files*/
char neemfile[200],gppmfile[200],poomfile[200],nppmfile[200];

/*pool change files*/
char dpoofile[200];


/*defining output files*/
strcpy(neefile,outputpath);strcat(neefile,"NEE.bin");
strcpy(rsofile,outputpath);strcat(rsofile,"RSO.bin");
strcpy(laifile,outputpath);strcat(laifile,"LAI.bin");
strcpy(gppfile,outputpath);strcat(gppfile,"GPP.bin");
strcpy(nppfile,outputpath);strcat(nppfile,"NPP.bin");
strcpy(poofile,outputpath);strcat(poofile,"POO.bin");
strcpy(parfile,outputpath);strcat(parfile,"PAR.bin");

/*defining mean output files*/
strcpy(neemfile,outputpath);strcat(neemfile,"NEEM.bin");
strcpy(gppmfile,outputpath);strcat(gppmfile,"GPPM.bin");
strcpy(nppmfile,outputpath);strcat(nppmfile,"NPPM.bin");
strcpy(poomfile,outputpath);strcat(poomfile,"POOM.bin");
/*pool change files*/
strcpy(dpoofile,outputpath);strcat(dpoofile,"DPOO.bin");


/*defining data structure and reading met file to fill it in!*/
DATA D;INITIALIZE_DATA_STRUCT(&D);
sprintf(filenum,"%05d",1);
strcpy(metfile,metpath);strcat(metfile,"DG_");strcat(metfile,filenum);strcat(metfile,".bin");
CARDAMOM_READ_BINARY_DATA(metfile,&D);

printf("maxnp = %d\n",maxnp);
double *M_LAI=malloc(D.nodays*sizeof(double)*maxnp);
double *M_GPP=malloc(D.nodays*sizeof(double)*maxnp);
double *M_NEE=malloc(D.nodays*sizeof(double)*maxnp);
double *M_FLUXES=malloc(D.nodays*sizeof(double)*nflu*maxnp);
double *M_POOLS=malloc((D.nodays+1)*sizeof(double)*npoo*maxnp);
int ND=D.nodays;

/*printing stuff to make sure everything's OK*/

double Pmax,Pthreshold,Wtot,W;




/*number of pixels*/
NPX=atoi(CLI[2]);

double *NEE,*GPP,*POO,*PAR,*LAI,*NPP,*RSO;
double *NEEU,*GPPU,*POOU,*PARU,*LAIU,*NPPU,*RSOU;
/*output fields*/
NEE=calloc(D.nodays*NPX,sizeof(double));
RSO=calloc(D.nodays*NPX,sizeof(double));
GPP=calloc(D.nodays*NPX,sizeof(double));
NPP=calloc(D.nodays*NPX,sizeof(double));
LAI=calloc(D.nodays*NPX,sizeof(double));
POO=calloc(D.nodays*NPX*npoo,sizeof(double));
PAR=calloc(npar*NPX,sizeof(double));
/*output field uncertainty*/
/*sqrt(total((a-mean(a)).^2/(numel(a)-1)))*/

NEEU=calloc(D.nodays*NPX,sizeof(double));
RSOU=calloc(D.nodays*NPX,sizeof(double));
GPPU=calloc(D.nodays*NPX,sizeof(double));
NPPU=calloc(D.nodays*NPX,sizeof(double));
LAIU=calloc(D.nodays*NPX,sizeof(double));
POOU=calloc(D.nodays*NPX*npoo,sizeof(double));
PARU=calloc(npar*NPX,sizeof(double));

/*************/

/*mean fluxes/pools w. uncertainty and overall change*/
double *NEEM,*GPPM,*POOM,*NPPM;
double *NEEMU,*GPPMU,*POOMU,*NPPMU;
double *meannee,*meangpp,*meannpp,*meanpoo;

/*output fields*/
NEEM=calloc(NPX,sizeof(double));
GPPM=calloc(NPX,sizeof(double));
NPPM=calloc(NPX,sizeof(double));
POOM=calloc(NPX*npoo,sizeof(double));
/*output field uncertainty*/
/*sqrt(total((a-mean(a)).^2/(numel(a)-1)))*/

NEEMU=calloc(NPX,sizeof(double));
GPPMU=calloc(NPX,sizeof(double));
NPPMU=calloc(NPX,sizeof(double));
POOMU=calloc(NPX*npoo,sizeof(double));

meannee=malloc(maxnp*sizeof(double));
meangpp=malloc(maxnp*sizeof(double));
meannpp=malloc(maxnp*sizeof(double));
meanpoo=malloc(maxnp*sizeof(double)*npoo);

/*pool changes for full time period*/
double *DPOO,*DPOOU,*dpoo;
DPOO=calloc(NPX*npoo,sizeof(double));
DPOOU=calloc(NPX*npoo,sizeof(double));
dpoo=malloc(maxnp*npoo*sizeof(double));


/************/


/*defining parameter structure*/
PARS PARS;
PARS.npar=npar;
PARS.pars=malloc(PARS.npar*maxnp*sizeof(double));
PARS.prob=malloc(maxnp*sizeof(double));
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
	CARDAMOM_READ_BINARY_DATA(metfile,&D);
	/*runnning all pars with met data*/
	/*looping through all parameters*/

	/*deriving max likelyhood - using this to determine if it is worth running model*/
	Pmax=PARS.prob[0];for (np=0;np<PARS.nsol;np++){if (PARS.prob[np]>Pmax){Pmax=PARS.prob[np];}}
	Pthreshold=Pmax-log(10000);Pthreshold=-1e6;
	Wtot=0;for (np=0;np<PARS.nsol;np++){if (PARS.prob[np]>Pthreshold){Wtot+=1;}}
	printf("Wtot = %f\n",Wtot);

	for (np=0;np<PARS.nsol;np++){
		if (PARS.prob[np]>Pthreshold){
		/*RUNNING DALEC*/
		DALEC_CDEA(D.MET, PARS.pars +np*PARS.npar, D.deltat,D.nodays, D.LAT, 
		M_LAI+ND*np, 
		M_NEE+ND*np, 	
		M_FLUXES+ND*np*nflu, 
		M_POOLS+(ND+1)*np*npoo);
		W=1;
		/*************store/write all data here***********/
	}}
	

	/*all dalec output is stored, loop 1 to work out ts and means here*/
	for (np=0;np<PARS.nsol;np++){
		if (PARS.prob[np]>Pthreshold){

		/*initialising temporary variables*/
		meannee[np]=0;meangpp[np]=0;meannpp[np]=0;
		for (p=0;p<npoo;p++){meanpoo[np*npoo+p]=0;dpoo[np*npoo+p]=0;}		


		/*fluxes*/
		for (m=0;m<D.nodays;m++){
		/*8day*/
		NEE[n*D.nodays+m]+=M_NEE[m+ND*np]*W/Wtot;
		RSO[n*D.nodays+m]+=M_FLUXES[m*nflu+13+(ND*np*nflu)]*W/Wtot;
		LAI[n*D.nodays+m]+=M_LAI[m+(ND*np)]*W/Wtot;
		GPP[n*D.nodays+m]+=M_FLUXES[m*nflu+(ND*np*nflu)]*W/Wtot;
		NPP[n*D.nodays+m]+=(M_FLUXES[m*nflu+(ND*np*nflu)]-M_FLUXES[m*nflu+2+(ND*np*nflu)])*W/Wtot;
		/*overall mean*/
		NEEM[n]+=M_NEE[m+ND*np]*W/Wtot/D.nodays;
                GPPM[n]+=M_FLUXES[m*nflu+(ND*np*nflu)]*W/Wtot/D.nodays;
                NPPM[n]+=(M_FLUXES[m*nflu+(ND*np*nflu)]-M_FLUXES[m*nflu+2+(ND*np*nflu)])*W/Wtot/D.nodays;
		
		/*mean for each solution - for mean uncertainty estimate later*/
		meannee[np]+=M_NEE[m+ND*np]/D.nodays;
		meangpp[np]+=M_FLUXES[m*nflu+ND*np*nflu]/D.nodays;
		meannpp[np]+=(M_FLUXES[m*nflu+(ND*np*nflu)]- M_FLUXES[m*nflu+2+(ND*np*nflu)])/D.nodays;
	        for (p=0;p<npoo;p++){meanpoo[np*npoo+p]+=M_POOLS[(m+1)*npoo+p+(ND+1)*np*npoo]/D.nodays;}


		/*pools*/
		for (p=0;p<npoo;p++){POO[p*D.nodays*NPX+n*D.nodays+m]+=M_POOLS[(m+1)*npoo+p+((ND+1)*npoo*np)]*W/Wtot;}
		for (p=0;p<npoo;p++){POOM[p*NPX+n]+=M_POOLS[(m+1)*npoo+p+((ND+1)*npoo*np)]*W/Wtot/D.nodays;}
		}
	
		/*calulating overall pool change	*/
		for (p=0;p<npoo;p++){dpoo[np*npoo+p]=M_POOLS[(ND)*npoo+p+(ND+1)*npoo*np]-M_POOLS[p+(ND+1)*npoo*np];
		DPOO[p*NPX+n]+=dpoo[np*npoo+p]*W/Wtot;}
		

		/*storing mean pars logarithmically*/
		/*Weigh correctly with parbrob*/
		for (npa=0;npa<PARS.npar;npa++){PAR[n*npar+npa]+=log(PARS.pars[PARS.npar*np+npa])*W/Wtot;}
	}}
	

	/*Loop 2: uncertainty calculations here*/

	for (np=0;np<PARS.nsol;np++){
		if (PARS.prob[np]>Pthreshold){
		/*flux uncertainty*/

		for (m=0;m<D.nodays;m++){
		/*8day uncertainty*/
		NEEU[n*D.nodays+m]+=pow(M_NEE[m+ND*np]-NEE[n*D.nodays+m],2)/(Wtot-1);
		RSOU[n*D.nodays+m]+=pow(M_FLUXES[m*nflu+13+(ND*np*nflu)]-RSO[n*D.nodays+m],2)/(Wtot-1);		
		LAIU[n*D.nodays+m]+=pow(M_LAI[m+(ND*np)]-LAI[n*D.nodays+m],2)/(Wtot-1);
		GPPU[n*D.nodays+m]+=pow(M_FLUXES[m*nflu+(ND*np*nflu)]-GPP[n*D.nodays+m],2)/(Wtot-1);
		NPPU[n*D.nodays+m]+=pow(M_FLUXES[m*nflu+(ND*np*nflu)]-M_FLUXES[m*nflu+2+(ND*np*nflu)]-NPP[n*D.nodays+m],2)/(Wtot-1);
		/*pools*/
		for (p=0;p<npoo;p++){
			POOU[p*D.nodays*NPX+n*D.nodays+m]+=pow(M_POOLS[(m+1)*npoo+p+((ND+1)*np*npoo)]-POO[p*D.nodays*NPX+n*D.nodays+m],2)/(Wtot-1);}
		
		}

		/*mean flux uncertainty*/
		NEEMU[n]+=pow(meannee[np]-NEEM[n],2)/(Wtot-1);		
		GPPMU[n]+=pow(meangpp[np]-GPPM[n],2)/(Wtot-1);		
		NPPMU[n]+=pow(meannpp[np]-NPPM[n],2)/(Wtot-1);		
		for (p=0;p<npoo;p++){POOMU[p*NPX+n]+=pow(meanpoo[np*npoo+p]-POOM[p*NPX+n],2)/(Wtot-1);}
		for (p=0;p<npoo;p++){DPOOU[p*NPX+n]+=pow(dpoo[np*npoo+p]-DPOO[p*NPX+n],2)/(Wtot-1);}



		/*storing mean pars logarithmically*/
		/*Weigh correctly with parbrob*/
		for (npa=0;npa<PARS.npar;npa++){PARU[n*npar+npa]+=pow(log(PARS.pars[PARS.npar*np+npa])-PAR[n*npar+npa],2)/(Wtot-1);}
	}}



      /*************DONE STORING D****************/

	}
/*done with if*/
}


/*calculating probabilities sum*/

/*transforming log averaged quantities*/
for (n=0;n<npoo*D.nodays*NPX;n++){POO[n]=POO[n];}
for (n=0;n<npar*NPX;n++){PAR[n]=exp(PAR[n]);}
for (n=0;n<D.nodays*NPX;n++){GPP[n]=GPP[n];}
for (n=0;n<D.nodays*NPX;n++){NPP[n]=NPP[n];}
for (n=0;n<D.nodays*NPX;n++){NEE[n]=NEE[n];}
for (n=0;n<D.nodays*NPX;n++){RSO[n]=RSO[n];}
for (n=0;n<D.nodays*NPX;n++){LAI[n]=LAI[n];}
/*corresponding uncertainties*/
for (n=0;n<npoo*D.nodays*NPX;n++){POOU[n]=exp(sqrt(POOU[n]));}
for (n=0;n<npar*NPX;n++){PARU[n]=exp(sqrt(PARU[n]));}
for (n=0;n<D.nodays*NPX;n++){GPPU[n]=sqrt(GPPU[n]);}
for (n=0;n<D.nodays*NPX;n++){NPPU[n]=sqrt(NPPU[n]);}
for (n=0;n<D.nodays*NPX;n++){NEEU[n]=sqrt(NEEU[n]);}
for (n=0;n<D.nodays*NPX;n++){RSOU[n]=sqrt(RSOU[n]);}
for (n=0;n<D.nodays*NPX;n++){LAIU[n]=sqrt(LAIU[n]);}

/*sqrt of mean flux and pool total variances*/
for (n=0;n<NPX;n++){NEEMU[n]=sqrt(NEEMU[n]);}
for (n=0;n<NPX;n++){GPPMU[n]=sqrt(GPPMU[n]);}
for (n=0;n<NPX;n++){NPPMU[n]=sqrt(NPPMU[n]);}
for (n=0;n<NPX*npoo;n++){POOMU[n]=sqrt(POOMU[n]);}
/*pool change*/
for (n=0;n<NPX*npoo;n++){DPOOU[n]=sqrt(DPOOU[n]);}




/*output files*/
FILE *neefd=fopen(neefile,"wb");
FILE *rsofd=fopen(rsofile,"wb");
FILE *laifd=fopen(laifile,"wb");
FILE *gppfd=fopen(gppfile,"wb");
FILE *nppfd=fopen(nppfile,"wb");
FILE *poofd=fopen(poofile,"wb");
FILE *parfd=fopen(parfile,"wb");


/*mean output files*/
FILE *neemfd=fopen(neemfile,"wb");
FILE *gppmfd=fopen(gppmfile,"wb");
FILE *nppmfd=fopen(nppmfile,"wb");
FILE *poomfd=fopen(poomfile,"wb");

/*pool change*/
FILE *dpoofd=fopen(dpoofile,"wb");

/**********/



/*writing all data*/
printf("D.nodays*NPX = %d\n", D.nodays*NPX);


fwrite(NEE,sizeof(double),D.nodays*NPX,neefd);
fwrite(RSO,sizeof(double),D.nodays*NPX,rsofd);
fwrite(LAI,sizeof(double),D.nodays*NPX,laifd);
fwrite(GPP,sizeof(double),D.nodays*NPX,gppfd);
fwrite(NPP,sizeof(double),D.nodays*NPX,nppfd);
fwrite(POO,sizeof(double),D.nodays*NPX*npoo,poofd);
fwrite(PAR,sizeof(double),npar*NPX,parfd);

/*writing mean pools*/
fwrite(NEEM,sizeof(double),NPX,neemfd);
fwrite(GPPM,sizeof(double),NPX,gppmfd);
fwrite(NPPM,sizeof(double),NPX,nppmfd);
fwrite(POOM,sizeof(double),NPX*npoo,poomfd);
/*writing dpoo*/
fwrite(DPOO,sizeof(double),NPX*npoo,dpoofd);


/*writing uncertainties in the same file*/
fwrite(NEEU,sizeof(double),D.nodays*NPX,neefd);fclose(neefd);
fwrite(RSOU,sizeof(double),D.nodays*NPX,rsofd);fclose(rsofd);
fwrite(LAIU,sizeof(double),D.nodays*NPX,laifd);fclose(laifd);
fwrite(GPPU,sizeof(double),D.nodays*NPX,gppfd);fclose(gppfd);
fwrite(NPPU,sizeof(double),D.nodays*NPX,nppfd);fclose(nppfd);
fwrite(POOU,sizeof(double),D.nodays*NPX*npoo,poofd);fclose(poofd);
fwrite(PARU,sizeof(double),npar*NPX,parfd);fclose(parfd);

/*writing mean field uncertainties*/
fwrite(NEEMU,sizeof(double),NPX,neemfd);fclose(neemfd);
fwrite(GPPMU,sizeof(double),NPX,gppmfd);fclose(gppmfd);
fwrite(NPPMU,sizeof(double),NPX,nppmfd);fclose(nppmfd);
fwrite(POOMU,sizeof(double),NPX*npoo,poomfd);fclose(poomfd);
/*writing dpoo unc.*/
fwrite(DPOOU,sizeof(double),NPX*npoo,dpoofd);fclose(dpoofd);


oksofar("All files written");

/*meanpoo[np*npoo+p]+=M_POOLS[(m+1)*npoo+p+(ND+1)*np*npoo]/D.nodays;*/

/*free all memory*/
free(PARS.pars);free(PARS.prob);
free(NEE);free(GPP);free(POO);free(PAR);free(LAI);free(NPP);free(RSO);
free(NEEU);free(GPPU);free(POOU);free(PARU);free(LAIU);free(NPPU);free(RSOU);


free(meannee);free(meangpp);free(meannpp);free(meanpoo);
/*freeing temp fpoo*/
free(dpoo);

free(NEEM);free(GPPM);free(POOM);free(NPPM);
free(DPOO);
free(NEEMU);free(GPPMU);free(POOMU);free(NPPMU);
free(DPOOU);


oksofar("about to free data struct");

free(M_LAI);
free(M_GPP);
free(M_NEE);
free(M_FLUXES);
free(M_POOLS);

FREE_DATA_STRUCT(D);


return 0;}




